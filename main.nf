#!/usr/bin/env nextflow

params.version = workflow.manifest.version

// Help message
helpMessage = """
===============================================================================
IKMB Virus pipeline | version ${params.version}
===============================================================================
Usage: nextflow run ikmb/virus-pipe --reads 'path/to/*_{1,2}_001.fastq.gz'
This example will perform assembly of viral reads, substracting any potential human reads using bloom filters
Required parameters:
--reads                        Input reads as a set of one or more PE Illumina reads
--email                        Email address to send reports to (enclosed in '')
Optional parameters:
--reference			A reference genome in fasta format to compare against
--iva_guided			Perform reference-guided instead of de-novo (default: false - experimental and may crash)
Output:
--outdir                       Local directory to which all output is written (default: results)
"""

params.help = false

// Show help when needed
if (params.help){
    log.info helpMessage
    exit 0
}


if (params.reference) {
	REF = file(params.reference)
	if (!REF.exists() ) {
		exit 1, "Reference file not found!"
	}
	Channel.fromPath(params.reference)
		.set { FastaIndex }

} else {
	REF = file("${baseDir}/assets/reference/NC_045512.2.fa")
	Channel.fromPath("${baseDir}/assets/reference/NC_045512.2.fa")
	.set { FastaIndex }
}

if (!REF.exists() ) {
	exit 1, "Could not find reference file..."
}

BLOOMFILTER = params.bloomfilter

OUTDIR = params.outdir

Channel.fromFilePairs(params.reads, flat: true)
	.ifEmpty { exit 1, "Did not find any reads matching your criteria!" }
	.set { reads_fastp }

process runFastp {

	scratch true

        input:
        set val(id), file(fastqR1),file(fastqR2) from reads_fastp

        output:
	set val(id),file(left),file(right) into (inputBioBloom , inputBwa)
        set file(html),file(json) into fastp_results

        script:

        left = fastqR1.getBaseName() + "_trimmed.fastq.gz"
        right = fastqR2.getBaseName() + "_trimmed.fastq.gz"
        json = fastqR1.getBaseName() + ".fastp.json"
        html = fastqR1.getBaseName() + ".fastp.html"

        """
                fastp --in1 $fastqR1 --in2 $fastqR2 --out1 $left --out2 $right --detect_adapter_for_pe -w ${task.cpus} -j $json -h $html --length_required 35
        """
	
}

process Bloomfilter {

	publishDir "${OUTDIR}/${id}/Bloomfilter", mode: 'copy'

	input:
	set id,file(left_reads),file(right_reads) from inputBioBloom

	output:
	set id,file(clean_reads) into inputIva
	set id,file(bloom) into BloomReport

	script:

	bloom = id + "_summary.tsv"
	clean_reads = id + ".filtered.fastq.gz"

	"""
		biobloomcategorizer -p $id --gz_output -d -n -e -s 0.01 -t ${task.cpus} -f "$BLOOMFILTER" $left_reads $right_reads | gzip > $clean_reads
	"""

}

process runIva {

        publishDir "${OUTDIR}/${id}/Assembly/Iva", mode: 'copy'

	input:
	set val(id),file(reads) from inputIva
	

	output:
	set val(id),file(contigs) into ( contigsIva, contigsAlign)

	script:

	outdir = "iva_" + id
	contigs = outdir + "/contigs.fasta"

	def options = ""
	if (params.iva_guided) {
		options = "--reference $REF"
	}

	"""
		iva --fr $reads $options $outdir
	"""

}

if (REF) {

	process alignRef {

		publishDir "${OUTDIR}/${id}/Assembly/Iva/Alignment", mode: 'copy'

		input:
		set val(id),file(contigs) from contigsAlign

		output:
		set val(id),file(msa) into contigsMSA

		script:
		ref_name = REF.getBaseName()
		merged_fa = id + "." + ref_name + ".fa"
		msa = merged_fa + ".aln"

		"""
			cat $REF $contigs >> $merged_fa
			muscle -in $merged_fa -out $msa -clwstrict
		"""

	}

	process makeBwaIndex {
		
		input:
		file(fasta) from FastaIndex

		output:
		set file(fasta),file(amb),file(ann),file(bwt),file(pac),file(sa) into BwaIndex
		
		script:
		base_name = fasta.getName()
		amb = base_name + ".amb"
		ann = base_name + ".ann"
		bwt = base_name + ".bwt"
		pac = base_name + ".pac"
		sa = base_name + ".sa"
		
		"""		
			bwa index $fasta
		"""
	}

	process runBwa {

		publishDir "${OUTDIR}/${sampleID}/BAM/raw", mode: 'copy'

		input:
		set val(sampleID),file(left),file(right) from inputBwa
		set file(fasta),file(amb),file(ann),file(bwt),file(pac),file(sa) from BwaIndex.collect()


		output:
		set val(sampleID),file(bam),file(bai) into bwaBam

		script:
		ref_name = REF.getBaseName()
		bam = sampleID + "." + ref_name + ".aligned.bam"
		bai = bam + ".bai"

		"""
			samtools dict $fasta > header.txt
			bwa mem -H header.txt -M -R "@RG\\tID:${sampleID}\\tPL:ILLUMINA\\tSM:${sampleID}\\tLB:${sampleID}\\tDS:${REF}\\tCN:CCGA" -t ${task.cpus} $fasta $left $right | samtools sort -O bam -m 2G -@ 4 - > $bam
			samtools index $bam
		"""
	}

	process runMD {
	
		publishDir "${OUTDIR}/${id}/BAM", mode: 'copy'
		
		input:
		set val(id),file(bam),file(bai) from bwaBam

		output:
		set val(id),file(bam_md),file(bai_md) into bamMD		

		script:
		bam_md = bam.getBaseName() + ".md.bam"
		bai_md = bam_md + ".bai"

		"""
			samtools sort -n $bam | samtools fixmate -m - fix.bam
			samtools sort -O BAM -o sorted.bam fix.bam
			samtools markdup sorted.bam $bam_md
			samtools index $bam_md
			rm fix.bam sorted.bam* 
		"""

	}

	process runFreebayes {

		publishDir "${OUTDIR}/${id}/Variants", mode: 'copy'

		input:
		set val(id),file(bam),file(bai) from bamMD

		output:
		set val(id),file(vcf) into fbVcf

		script:
		base_name = bam.getBaseName()
		vcf = base_name + ".vcf"

		"""
			freebayes --ploidy 1 -f $REF --genotype-qualities $bam > $vcf
		"""

	}

	process runFilterVcf {

		publishDir "${OUTDIR}/${id}/Variants", mode: 'copy'

		input:
		set val(id),file(vcf) from fbVcf

		output:
		file(vcf_filtered) into finalVcf

		script:
		vcf_filtered = vcf.getBaseName() + ".filtered.vcf"

		"""
			bcftools filter ${params.filter_options} $vcf > $vcf_filtered
		"""
	}

}

