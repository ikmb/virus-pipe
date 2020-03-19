#!/usr/bin/env nextflow

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
	REF = false
	FastaIndex = Channel.empty()
}

BLOOMFILTER = params.bloomfilter

OUTDIR = params.outdir

Channel.fromReadPairs(params.reads)
	.ifEmpty { exit 1, "Did not find any reads matching your criteria!" }
	.into { reads_fastp }

process runFastp {

	scratch true

        input:
        set val(id), file(fastqR1),file(fastqR2) from reads_fastp

        output:
	
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
	set id,file(left_reads),file(right_reads) into outputBloomfilter
	set id,file(bloom) into BloomReport

	script:

	bloom = id + "_summary.tsv"
	clean_reads = id + ".filtered.fastq.gz"

	"""
		biobloomcategorizer -p $id -gz_output -d -n -e -s 0.01 -t ${task.cpus} -f "$BLOOMFILTER" $left_reads $right_reads | gzip > $clean_reads
	"""

}

process runIva {

        publishDir "${OUTDIR}/${id}/Assembly/Iva", mode: 'copy'

	input:
	set val(id),file(reads) from inputIva
	

	output:
	set val(id),file(contigs) into contifsIva

	script:

	outdir = "iva_" + id
	contigs = outdir + "/contigs.fasta"

	def options = ""
	if (params.reference) {
		options = "--reference ${params.reference}"
	}

	"""
		iva --fr $reads $options $outdir
	"""

}

if (REF) {

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

		input:
		set val(id),file(right),file(left) from inputBwa
		set file(fasta),file(amb),file(ann),file(bwt),file(pac),file(sa) from BwaIndex.collect()


		output:
		set val(id),file(bam),file(bai) into bwaBam

		script:
		ref_name = REF.getBaseName()
		bam = id + "." ref_name + ".aligned.bam"
		bai = bam + ".bai"

		"""
			samtools dict $fasta > header.txt
			bwa mem -H header.txt -t ${task.cpus} $fasta | samtools sort  -O bam -m 2G -@ 4 - > $bam
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
			samtools markdup $bam $bam_md
			samtools index $bam_md
		"""

	}

	process runFreebayes {

		publishDir "${OUTDIR}/${id}/Variants", mode: 'copy'

		input:
		set val(id),file(bam),file(bai) from baMD

		output:
		set val(id),file(vcf) into fbVcf

		script:
		base_name = bam.getBaseName()
		vcf = base_name + ".vcf"

		"""
			freebayes -f $REF --genotype-qualities $bam > $vcf
		"""

	}

}

