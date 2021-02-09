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
--reads                         Input reads as a set of one or more PE Illumina reads (or:...)
--samples			Sample sheet with additional sample information (instead of --reads). See github for formatting hints.
--email                         Email address to send reports to (enclosed in '')
Optional parameters:
--clip				Remove x bases from both ends of the reads (default: 20)
--var_call_cov			Coverage of a site to be considered in variant calling (default: 20)
--var_call_frac			Fraction of reads required to support a variant call (default: 0.1)
--var_call_count		Number of reads required to support a variant call (default: 10)
--var_filter_mqm		Mean mapping quality for a variant to survive filtering (default: 40)
--var_filter_qual		Call quality for variant to survive filtering (default: 20)
--var_filter_sap		Read strand bias for variant to survive filtering (default: 100)
--primer_set                    Name of the primer set used (ARTIC-v3, Eden)
--primer_fasta			Primer sequences in FASTA format (overrides --primer_set)
--run_name			Specify a name for this analysis run
--kraken2_db			A kraken2-formatted database with virus species for taxonomic mapping
--assemble			Assemble genomes de-novo (in addition to the ref-flipping way)
--guided 			Guide assembly with known reference
--email				Specify email to send report to
Output:
--outdir                       Local directory to which all output is written (default: results)
"""

params.help = false

/*
*/

// Show help when needed
if (params.help){
    log.info helpMessage
    exit 0
}

def summary = [:]

REF = file("${baseDir}/assets/reference/NC_045512.2.fa")
REF_GFF = file("${baseDir}/assets/reference/NC_045512.2.gff")
REF_NAME = "NC_045512.2"

REF_WITH_HOST = file(params.ref_with_host)

host_genome = Channel.fromPath("${params.host_index}*")

// Minimum file size in bytes to attempt assembly from
size_limit = params.size_limit

/*
Primer sequences
*/

// Selection of amplicon primers
if (params.primer_fasta) {
	primers = Channel.fromPath(file(params.primer_fasta))
} else if (params.primer_set) {
	if (params.primer_sets[params.primer_set]) {
		primers = Channel.fromPath(params.primer_sets[params.primer_set].fasta)
} else {
		exit 1, "Did not recognize primer set name"
	}
} else {
	primers = Channel.empty()
}  

if (params.guided && !params.assemble) {
	log.info "Requested a guided assembly but did not actually enable assembly - will do that for you now!"
	params.assemble = true
}

if (params.filter && params.fast_filter) {
	log.info "Requested filter and fast_filter - will only use fast_filter..."
	params.filter = false
}

// set basic global options like location of database files
KRAKEN2_DB = params.kraken2_db

PATHOSCOPE_INDEX_DIR=file(params.pathoscope_index_dir)

OUTDIR = params.outdir

run_name = ( !params.run_name) ? "${workflow.sessionId}" : "${params.run_name}"

if (!params.run_name ) {
        log.info "No run name was specified, using ${run_name} instead"
}

if (params.samples && params.reads) {
	log.info "Specified both --reads and --samples! Will only consider --samples"
}

summary['Reference'] = REF
summary['Kraken2DB'] = params.kraken2_db
summary['PathoscopeDB'] = params.pathoscope_index_dir
summary['MappingReference'] = REF
summary['VariantCalling'] = [:]

summary['VariantCalling']['VarCallCov'] = params.var_call_cov
summary['VariantCalling']['VarCallFrac'] = params.var_call_frac
summary['VariantCalling']['VarCallCount'] = params.var_call_count
summary['VariantCalling']['VarFilterMqm'] = params.var_filter_mqm
summary['VariantCalling']['VarFilterSap'] = params.var_filter_sap
summary['VariantCalling']['VarFilterQual'] = params.var_filter_qual
summary['VariantCalling']['ConsensusMinCov'] = params.cns_min_cov
summary['VariantCalling']['ConsensusGtAdjust'] = params.cns_gt_adjust

// Header log info
log.info "IKMB ------------------------------------------------------------------------------"
log.info "db    db d888888b d8888b. db    db .d8888.        d8888b. d888888b d8888b. d88888b "
log.info "88    88   `88'   88  `8D 88    88 88'  YP        88  `8D   `88'   88  `8D 88'     "
log.info "Y8    8P    88    88oobY' 88    88 `8bo.          88oodD'    88    88oodD' 88ooooo "
log.info "`8b  d8'    88    88`8b   88    88   `Y8b. C8888D 88~~~      88    88~~~   88~~~~~ "
log.info ".`8bd8'    .88.   88 `88. 88b  d88 db   8D        88        .88.   88      88.     "
log.info "...YP    Y888888P 88   YD ~Y8888P' `8888Y'        88      Y888888P 88      Y88888P "
log.info "==================================================================================="
log.info "${workflow.manifest.description}		v${params.version}"
log.info "Nextflow Version:             	$workflow.nextflow.version"
log.info "Mapping reference:		${params.ref_with_host}"
log.info "Viral reference:             	${REF}"
log.info "Virus DB:			${params.kraken2_db}"
log.info "Assemble de-novo:		${params.assemble}"
if (params.assemble) {
	log.info "Assemble guided:		${params.guided}"
}
if (params.primer_fasta) {
	log.info "Primers for trimming:		${params.primer_fasta}"
} else if (params.primer_set) {
	log.info "Primers for trimming:		${params.primer_set}"
} else {
	log.info "No designated primers for trimming - only Illumina primers used"
}
log.info "Read clipping 3'/5'		${params.clip}bp"
log.info "Command Line:			$workflow.commandLine"
if (workflow.containerEngine) {
        log.info "Container engine:		${workflow.containerEngine}"
}
log.info "==================================================================================="

// ********************
// WORKFLOW STARTS HERE
// ********************

// Helper function for the sample sheet parsing to produce sane channel elements
def returnFile(it) {
    // Return file if it exists
    inputFile = file(it)
    if (!file(inputFile).exists()) exit 1, "Missing file in TSV file: ${inputFile}, see --help for more information"
    return inputFile
}

Channel.fromPath(REF)
	.into {  inputNormalize; Ref2Consensus; Ref2BwaIdx ; inputBloomMaker; Ref2Freebayes }

if (params.samples) {
	Channel.from(file(params.samples))
        .splitCsv(sep: ';', header: true)
	.map { row ->
                        def patient = row.IndivID
			def sample = row.SampleID
                        def left = returnFile( row.R1 )
			def right = returnFile( row.R2)
                        [ patient, sample, left, right ]
                }
       .set {  reads_fastp }
} else if (params.reads) {
        Channel.fromFilePairs(params.reads, flat: true)
        .map { triple -> tuple( triple[0],triple[0],triple[1],triple[2]) }
        .set { reads_fastp }
}

process get_pangolin_version {

	executor 'local'

	label 'pangolin'

	output:
	file(pangolin_version) into pango_version

	script:
	pangolin_version = "v_pangolin.txt"

	"""
		pangolin -v > $pangolin_version
	"""
}

process get_software_versions {

    label 'std'

    publishDir "${OUTDIR}/Summary/versions", mode: 'copy'

    input:
    file(pangolin_version) from pango_version

    output:
    file("v*.txt")
    file(yaml_file) into software_versions_yaml
    file(tab_file) into software_versions_report

    script:
    yaml_file = "software_versions_mqc.yaml"
    tab_file = "software_versions.tab"

    """
	    echo $workflow.manifest.version &> v_ikmb_virus_pipe.txt
	    echo $workflow.nextflow.version &> v_nextflow.txt
	    echo "Kraken2 2.0.8_beta" > v_kraken2.txt
	    freebayes --version &> v_freebayes.txt
	    fastp -v &> v_fastp.txt
	    samtools --version &> v_samtools.txt
	    bcftools --version &> v_bcftools.txt
	    multiqc --version &> v_multiqc.txt
	    bwa &> v_bwa.txt 2>&1 || true	    
	    bowtie2 --version &> v_bowtie2.txt
	    parse_versions.pl >  $yaml_file
	    parse_versions_tab.pl > $tab_file
    """
}

process trim_reads {

	label 'fastp'

	scratch true

        publishDir "${OUTDIR}/${sampleID}/RawReads", mode: 'copy'

        input:
        set val(patientID),val(sampleID), file(fastqR1),file(fastqR2) from reads_fastp
	file(primer_fa) from primers.collect().ifEmpty(false)

        output:
	set val(sampleID),file(left),file(right) into inputBioBloomTarget 
	set val(patientID),val(sampleID),file(left),file(right) into ( inputBwa, inputBowtieFilter, inputBioBloomHost )
        set file(html),file(json) into fastp_qc

        script:
	def options = ""
	if (primer_fa) {
		options = "--adapter_fasta " + primer_fa
	}

        left = fastqR1.getBaseName() + "_trimmed.fastq.gz"
        right = fastqR2.getBaseName() + "_trimmed.fastq.gz"
        json = sampleID + ".fastp.json"
        html = sampleID + ".fastp.html"

        """
                fastp --in1 $fastqR1 --in2 $fastqR2 --out1 $left --out2 $right $options -f $params.clip -t $params.clip --detect_adapter_for_pe -w ${task.cpus} -j $json -h $html --length_required 35
        """
	
}

// Get the fraction of Sars-CoV2 reads against this bloom filter
process make_bloomfilter {

	label 'std'

        publishDir "${OUTDIR}/Reference/Bloomfilter", mode: 'copy'

        input:
        file(reference) from inputBloomMaker

        output:
        set file(bf),file(txt) into refBloom
        file(reference)

        script:
        base_name = reference.getBaseName()
        bf = base_name + ".bf"
        txt = base_name + ".txt"

        """
        	biobloommaker -p $base_name $reference
        """

}

// *****************************
// Filter reads against the human genome
// *****************************

if (params.filter) {

        // ****************
        // Using Bowtie2
        // ****************

        process remove_host_reads_bt {

                label 'std'

                publishDir "${OUTDIR}/${id}/CleanReads", mode: 'copy'

                input:
                set val(patientID),val(id),file(left_reads),file(right_reads) from inputBowtieFilter
                path(bwt_files) from host_genome.collect()

                output:
                set val(patientID),val(id),file(left_clean),file(right_clean) into ( inputPathoscopeMap, inputKraken, inputSpades, inputFailed  )

                script:
                left_clean = id + ".clean.R1.fastq.gz"
                right_clean = id + ".clean.R2.fastq.gz"
                unpaired_clean = id + ".clean.unpaired.fastq.gz"
                bowtie_log = id + ".txt"

                """
                        bowtie2 -x genome -1 $left_reads -2 $right_reads -S /dev/null --no-unal -p ${task.cpus} --un-gz $unpaired_clean  --un-conc-gz ${id}.clean.R%.fastq.gz 2> $bowtie_log
                """

        }


} else {
	inputBioBloomHost.into { inputKraken; inputPathoscopeMap; inputSpades; inputFailed  }
}

// ************************
// Get all the reads for a target species
// ************************
process select_covid_reads {

	label 'std'

	publishDir "${OUTDIR}/${id}/Bloomfilter/Target", mode: 'copy'

	input:
	set file(bf),file(txt) from refBloom.collect()
        set id,file(left_reads),file(right_reads) from inputBioBloomTarget

        output:
        set id,file(clean_reads)
        file(bloom) into BloomReportTarget

        script:
        analysis = id + ".Target"
        bloom = analysis + "_summary.tsv"
        clean_reads = analysis + ".filtered.fastq.gz"

        """
                biobloomcategorizer -p $analysis --gz_output -d -e -s 0.01 -t ${task.cpus} -f "$bf" $left_reads $right_reads | gzip > $clean_reads
        """
}

// *************************
// Get a list of non-host species in the sample
// *************************
process kraken2_search {

	label 'kraken'

	publishDir "${OUTDIR}/${id}/Taxonomy/", mode: 'copy'

	when:
	params.taxonomy

	input:
	set val(patientID),val(id),file(left),file(right) from inputKraken

	output:
	set val(patientID),val(id),file(report) into (KrakenReport, Kraken2Report)
	file(kraken_log)

	script:
	report = id + ".kraken2_report.txt"
	kraken_log = id + ".kraken2.log"

	"""
		kraken2 --db $KRAKEN2_DB --threads ${task.cpus} --output $kraken_log --report $report $left $right 
		if [ ! -f $report ]; then
			touch $report
		fi
	"""
}

process kraken2yaml {


	input:
	file(reports) from KrakenReport.collect()

	output:
	file(report_yaml) into KrakenYaml

	script:
	
	report_yaml = "kraken_report_mqc.yaml"

	"""
		kraken_covid2yaml.pl --outfile $report_yaml
	"""

}

process denovo_assemble_virus {

	label 'std'

	publishDir "${OUTDIR}/${id}/Denovo_Assembly", mode: 'copy'

	when:
	params.assemble

	input:
	set val(patientID),val(id),file(left),file(right) from inputSpades.filter{ p,i,l,r -> r.size() > 500000 }

	output:
	set val(patientID),val(id),file(assembly) optional true into contigsSpades
	
	script:
	assembly = id + ".spades.fasta"
	def options = ""
	if (params.guided) {
		options = "--untrusted-contigs ${REF}"
	}

	"""
		coronaspades.py -1 $left -2 $right $options -t ${task.cpus} -m ${task.memory.toGiga()} -o spades
		cp spades/scaffolds.fasta $assembly

	"""
}

process fail_sample {


	publishDir "${OUTDIR}/Failed", mode: 'copy'

	input:
	set val(patientID),val(id),file(left),file(right) from inputFailed.filter{ p,i,l,r -> r.size() < 500000 }

	output:
	file(failed_sample)

	script:
	failed_sample = id + ".failed.txt"

	"""
		echo "Failing $id due to low read count" > $failed_sample
	"""

}

process denovo_assembly_scaffold {

	label 'std'

        publishDir "${OUTDIR}/${id}/Denovo_Assembly", mode: 'copy'
	
	when:
	params.assemble

	input:
	set val(patientID),val(id),file(assembly) from contigsSpades

	output:
	set val(patientID),val(id),file(pseudo) into contigsScaffolds

	script:

	pseudo = id + ".spades.scaffolds.fasta"

	"""
		ragtag.py scaffold $REF $assembly -f 800 --remove-small -C -t ${task.cpus} -o ragtag

		cp ragtag/ragtag.scaffolds.fasta $pseudo
	"""
}

// **************************
// Run assembly QC with QUAST
// **************************

contigsScaffolds.into {assemblies; assemblies_qc }

process assembly_qc {

	label 'quast'

	publishDir "${OUTDIR}/${id}/Assembly/Spades/QC", mode: 'copy'

	input:
	set val(patientID),val(id),path(assembly) from assemblies_qc

	output:
	path(quast_dir) into QuastReport
	set val(patientID),val(id),file("${quast_dir}/report.tsv") into Quast2Report
	file("quast_results")

	script:
	named_assembly = id + ".fasta"
	quast_dir = id + "_quast"

	"""
		cp $assembly $named_assembly
		quast $named_assembly -r $REF -g $REF_GFF
		mkdir -p $quast_dir
		
		sed 's/_L/-L/' quast_results/latest/report.tsv | sed 's/COV_/COV-/' > ${quast_dir}/report.tsv
		
	"""
}
/*
*/

if (params.pathoscope) {
	// **********************
	// Run pathoscope MAP
	// **********************
	process pathoscope_map {

		label 'pathoscope'

		input:
		set val(patientID),val(id),file(left_reads),file(right_reads) from inputPathoscopeMap

		output:
		set val(patientID),val(id),file(pathoscope_sam) into inputPathoscopeId

		script:
		pathoscope_sam = id + ".sam"

		"""
        		pathoscope MAP -1 $left_reads -2 $right_reads -indexDir $PATHOSCOPE_INDEX_DIR -filterIndexPrefixes hg19_rRNA \
	        	-targetIndexPrefix A-Lbacteria.fa,M-Zbacteria.fa,virus.fa -outAlign $pathoscope_sam -expTag $id -numThreads ${task.cpus}
		"""

	}

	// **********************
	// Run pathoscope ID
	// **********************
	process pathoscope_id {

		label 'pathoscope'

		publishDir "${OUTDIR}/${id}/Pathoscope", mode: 'copy'

		input:
		set val(patientID),val(id),file(samfile) from inputPathoscopeId.filter{ i,r -> r.size() > 1000000 }

		output:
		set val(patientID),val(id),file(pathoscope_tsv) into outputPathoscopeId

		script:

		//pathoscope_sam = "updated_" + samfile
		pathoscope_tsv = id + "-sam-report.tsv"

		"""
        		pathoscope ID -alignFile $samfile -fileType sam -expTag $id
		"""

	}

} // end pathoscope

// **************************
// RKI-compliant: align reads against reference genome and flip variant bases
// **************************

process align_make_index_bwa {

	label 'std'

	publishDir "${OUTDIR}/Reference/BWA", mode: 'copy'

	input:
	file(fasta) from Ref2BwaIdx

	output:
	set file(fasta),file(bwa_amb),file(bwa_ann),file(bwa_btw),file(bwa_pac),file(bwa_sa) into BwaIndex

	script:
	base_name = fasta.getName()
	bwa_amb = base_name + ".amb"
	bwa_ann = base_name + ".ann"
	bwa_btw = base_name + ".bwt"
	bwa_pac = base_name + ".pac"
	bwa_sa = base_name + ".sa"
	
	"""
		bwa index $fasta
	"""
}

process align_viral_reads_bwa {

	label 'std'

	input:
	set val(patientID),val(sampleID),file(left),file(right) from inputBwa
	set file(fasta),file(bwa_amb),file(bwa_ann),file(bwa_btw),file(bwa_pac),file(bwa_sa) from BwaIndex.collect()

	output:
	set val(patientID),val(sampleID),file(bam),file(bai) into bwaBam

	script:
	ref_name = fasta.getBaseName()
        bam = sampleID + "." + ref_name + ".aligned.bam"
        bai = bam + ".bai"	

	"""
		samtools dict -a $ref_name -o assembly.dict -s Sars-CoV2 $fasta
		bwa mem -H assembly.dict -M -R "@RG\\tID:${sampleID}\\tPL:ILLUMINA\\tSM:${sampleID}\\tLB:${sampleID}\\tDS:${fasta}\\tCN:CCGA" \
                        -t ${task.cpus} ${fasta} $left $right \
                        | samtools sort -O bam -o $bam -
		samtools index $bam
	"""
}

// ***********************
// Mark duplicate reads
// ***********************
process mark_dups {

       	label 'std'

	publishDir "${OUTDIR}/${id}/BAM", mode: 'copy'

	input:
	set val(patientID),val(id),file(bam),file(bai) from bwaBam

	output:
	set val(patientID),val(id),file(bam_md_virus),file(bai_md_virus) into ( inputBamStats, inputBamCoverage, bam2mask)
	set val(patientID),val(id),file(bam_md_virus),file(bai_md_virus) into bamMD

	script:
	bam_md = id + ".dedup.bam"
	bai_md = bam_md + ".bai"
	bam_md_virus = id + ".dedup.viral_mapping.bam"
	bai_md_virus = bam_md_virus + ".bai"

	"""
		samtools sort -m 4G -t 2 -n $bam | samtools fixmate -m - fix.bam
		samtools sort -m 4G -t 2 -O BAM fix.bam | samtools markdup - tmp.md.bam
		samtools rmdup tmp.md.bam $bam_md
		samtools index $bam_md

		samtools view -bh -o $bam_md_virus $bam_md $REF_NAME
		samtools index $bam_md_virus
		rm fix.bam 
	"""

}

// ************************
// Coverage statistics
// ************************

process coverage_stats {
	
	label 'std'

	publishDir "${OUTDIR}/${id}/BAM", mode: 'copy'

	input:
	set val(patientID),val(id),file(bam),file(bai) from inputBamStats.filter{ p,i,b,d -> b.size() > 10000 }

	output:
	file(global_dist) into BamStats
	set val(patientID),val(id),file(sam_coverage) into BamCoverage
	set val(patientID),val(id),file(report),file(global_dist) into coverage_report

	script:
	global_dist = id + ".mosdepth.global.dist.txt"
	sam_coverage = id + ".coverage.samtools.txt"
	report = id + ".jpg"
	base_name = id 
	
	"""
		mosdepth -t ${task.cpus} $id $bam
		samtools depth -d 200  $bam > $sam_coverage
		bam2coverage_plot.R $sam_coverage $base_name
	"""
	
}

process align_stats {

        label 'std'

        publishDir "${OUTDIR}/${sampleID}/BAM", mode: 'copy'

        input:
        set val(patientID),val(sampleID),file(bam),file(bai) from inputBamCoverage

	output:
	file(align_stats) into BamAlignStats
	set val(patientID),val(sampleID),file(align_stats) into Samtools2Report

	script:
	align_stats = sampleID + ".txt"

	"""
		samtools stats  $bam > $align_stats
	"""

}

// ************************
// Variant calling
// ************************
process call_variants {
	
       	label 'std'

	publishDir "${OUTDIR}/${id}/Variants/Raw", mode: 'copy'

	input:
	set val(patientID),val(id),file(bam),file(bai) from bamMD
	file(reference) from Ref2Freebayes.collect()

	output:
	set val(patientID),val(id),file(vcf) into fbVcf

	script:
	base_name = bam.getBaseName()
	vcf = base_name + ".vcf"

	"""
		freebayes --genotype-qualities \
			--min-coverage $params.var_call_cov \
			--haplotype-length -1 \
			--min-alternate-fraction $params.var_call_frac \
			--min-alternate-count $params.var_call_count \
			--pooled-continuous \
			-f $reference $bam > $vcf
	"""

}

// ********************
// Filter variant calls
// ********************
process filter_vcf {
	
       	label 'std'

	//publishDir "${OUTDIR}/${id}/Variants", mode: 'copy'

	input:
	set val(patientID),val(id),file(vcf) from fbVcf

	output:
	set val(patientID),val(id),file(vcf_filtered) into VcfNormalize

	script:
	vcf_filtered = vcf.getBaseName() + ".filtered.vcf"

	"""
		bcftools filter -e 'INFO/MQM < ${params.var_filter_mqm} | INFO/SAP > ${params.var_filter_sap} | QUAL < ${params.var_filter_qual}'  $vcf > $vcf_filtered
	"""
}

process normalize_and_adjust_vcf {

        publishDir "${OUTDIR}/${id}/Variants", mode: 'copy'

	label 'std'

	input:
	set val(patientID),val(id),file(vcf) from VcfNormalize
        file(ref_genome) from inputNormalize.collect()

	output:
        set val(patientID),val(id),file(vcf_filtered) into (finalVcf,Vcf2Report,VcfPredict)
	set val(patientID),val(id),file(vcf_filtered_gz),file(vcf_filtered_tbi) into (Vcf2Mask,Vcf2Consensus)

	script:

	vcf_filtered = vcf.getBaseName() + ".normalized.vcf"
	vcf_filtered_gz = vcf_filtered + ".gz"
	vcf_filtered_tbi = vcf_filtered_gz + ".tbi"

	"""
		vt normalize -o tmp.vcf -r $ref_genome $vcf

		adjust_gt_rki.py -o $vcf_filtered --vf $params.cns_gt_adjust $vcf

		bgzip -c $vcf_filtered > $vcf_filtered_gz
		tabix $vcf_filtered_gz
		
	"""

}

// *************************
// Make consensus assembly
// *************************

inputMasking = bam2mask.join(Vcf2Mask, by: [0,1] )

process create_cov_mask {

	publishDir "${OUTDIR}/${sampleID}/BAM", mode: 'copy'	

	label 'bedtools'
	
	input:
	set val(patientID),val(sampleID),file(bam),file(bai),file(vcf),file(tbi) from inputMasking

	output:
	set val(patientID),val(sampleID),file(mask) into ConsensusMask

	script:

	mask = sampleID + ".coverage_mask.bed"

	"""

		bedtools genomecov -bga -ibam $bam | awk '\$4 < ${params.cns_min_cov}' | bedtools merge > tmp.bed
		bedtools intersect -v -a tmp.bed -b $vcf > $mask
	"""

}

MakeConsensus = Vcf2Consensus.join(ConsensusMask, by: [0,1])

process consensus_assembly {

	//publishDir "${OUTDIR}/RKI_Assemblies", mode: 'copy'

	label 'std'

	input:
	set val(patientID),val(id),file(vcf),file(tbi),file(bed) from MakeConsensus
	file(ref_assembly) from Ref2Consensus.collect()
	
	output:
	set val(patientID),val(id),file(consensus) into Consensus2Header

	script:

	consensus = id + ".consensus_assembly.fa"
	
	"""
		bcftools consensus -I \
			-o $consensus \
			-f $ref_assembly \
			-m $bed \
			--sample $id \
			$vcf 
	"""
}

// Replace the fasta header to include sample name
process consensus_header {

	publishDir "${OUTDIR}/RKI_Assemblies", mode: 'copy'

	label 'std'

	input:
	set val(patientID),val(id),file(consensus) from Consensus2Header

	output:
	set val(patientID),val(id),file(consensus_reheader) into (assemblies_pangolin, consensus2qc)
	file(consensus_masked_reheader)

	script:

	consensus_reheader = id + ".iupac_consensus.fasta"
	consensus_masked_reheader = id + ".masked_consensus.fasta"

	header = id + "_iupac_consensus_v" + params.version 
	masked_header = id + "_masked_consensus_v" + params.version

	"""
		echo '>$header' > $consensus_reheader
		tail -n +2 $consensus | fold -w 80 >> $consensus_reheader
		echo  >> $consensus_reheader

		echo '>$masked_header' > $consensus_masked_reheader
		
		tail -n +2 $consensus | tr "RYSWKMBDHVN" "N" | fold -w 80 >> $consensus_masked_reheader
		echo >> $consensus_masked_reheader
	"""

}

// QC of the final IUPAC assembly
process consensus_qc {

	label 'gaas'

	publishDir "${OUTDIR}/${id}/QC", mode: 'copy'

	input:
	set val(patientID),val(id),file(consensus_reheader)  from consensus2qc

	output:
	set val(patientID),val(id),file(stats) into ConsensusStats

	script:
	stats = id + "_assembly_report.txt"

	"""
		gaas_fasta_statistics.pl -f $consensus_reheader -o stats > $stats
	"""
}

// **********************
// Determine Pangolin lineage
// **********************
process assembly_pangolin {

        label 'pangolin'

        publishDir "${OUTDIR}/${id}/Pangolin", mode: 'copy'

        input:
        set val(patientID),val(id),path(assembly) from assemblies_pangolin

        output:
        set val(patientID),val(id),path(report) into (pangolin_report, Pangolin2Report)

        script:

        report = id + ".pangolin.csv"

        """
                pangolin --outfile $report $assembly
        """
}

// Convert Pangolin results to YAML format
process pangolin2yaml {

        label 'std'

        publishDir "${OUTDIR}/Pangolin", mode: 'copy'

        input:
        file(reports) from pangolin_report.collect()

        output:
        file(report) into PangolinYaml
	file(xls)
	file(csv)

        script:

        report = "pangolin_report_mqc.yaml"
	xls = "pangolin." + run_name + ".xlsx"
	csv = "pangolin." + run_name + ".csv"

        """
                pangolin2yaml.pl > $report
		pangolin2xls.pl --outfile $xls > $csv
        """

}

process vcf_stats {

	label 'std'

	publishDir "${OUTDIR}/${id}/Variants", mode: 'copy'

	input:
	set val(patientID),val(id),file(vcf) from finalVcf

	output:
	file(stats) into VcfStats
	set val(patientID),val(id),file(stats) into VcfReport

	script:
	stats = vcf.getBaseName() + ".stats"

	"""
		bcftools stats $vcf > $stats
	"""

}

process effect_prediction {

	label 'std'

	publishDir "${OUTDIR}/${id}/Variants", mode: 'copy'

	input:
	set val(patientID),val(id),file(vcf) from VcfPredict

	output:
	set val(patientID),val(id),file(effects) into EffectPrediction

	script:

	effects = id + ".snpeff." + REF_NAME + ".vcf"

	"""
		snpEff -v $REF_NAME -onlyProtein -no-upstream -no-downstream -canon $vcf > $effects
	"""

}

// **********************
// Write a per-patient report
// **********************

GroupedReports = Kraken2Report.join(Pangolin2Report,by: [0,1]).join(Samtools2Report,by: [0,1]).join(ConsensusStats,by: [0,1]).join(EffectPrediction,by: [0,1]).join(coverage_report,by: [0,1])

process final_report {

	label 'std'

	publishDir "${OUTDIR}/Reports", mode: 'copy'

	input:
	set val(patientID),val(id),file(kraken),file(pangolin),file(samtools),file(fasta_qc),file(variants),file(coverage_plot),file(mosdepth) from GroupedReports
	file(version_yaml) from software_versions_report

	output:
	file(patient_report) 
	file(patient_report_json)

	script:

	patient_report = id + "_report.pdf"
	patient_report_json = id + "_report.json"

	"""
		cp $baseDir/assets/ikmb_bfx_logo.jpg . 
		covid_report.pl --patient $patientID --kraken $kraken --software $version_yaml --pangolin $pangolin --depth $mosdepth --bam_stats $samtools --assembly_stats $fasta_qc --vcf $variants --plot $coverage_plot --outfile $patient_report > $patient_report_json
	"""

}

// **********************
// Compile quality metrics into report
// **********************
process MultiQC {

	label 'std'

	publishDir "${OUTDIR}/MultiQC", mode: 'copy'

	input:
	file('*') from fastp_qc.collect().ifEmpty('')
	file('*') from BamStats.collect().ifEmpty('')
	file('*') from BloomReportTarget.collect().ifEmpty('')
	file('*') from KrakenYaml.ifEmpty('')
	file('*') from BamAlignStats.collect()
	file('*') from QuastReport.collect().ifEmpty('')
	file('*') from PangolinYaml.ifEmpty('')
	file('*') from VcfStats.collect()
	file('*') from software_versions_yaml.collect()

	output:
	file(report) into multiqc_report

	script:
	def options = ""
	if (params.assemble) {
		options = "*/report.tsv"
	}
	report = run_name + "_multiqc.html"
	"""
		cp $baseDir/assets/multiqc_config.yaml multiqc_config.yaml
		cp $params.logo .

		multiqc -n $report *.* $options
	"""

}

workflow.onComplete {
  log.info "========================================="
  log.info "Duration:           $workflow.duration"
  log.info "========================================="

  def email_fields = [:]
  email_fields['version'] = workflow.manifest.version
  email_fields['session'] = workflow.sessionId
  email_fields['runName'] = run_name
  email_fields['Reads'] = params.reads
  email_fields['success'] = workflow.success
  email_fields['dateStarted'] = workflow.start
  email_fields['dateComplete'] = workflow.complete
  email_fields['duration'] = workflow.duration
  email_fields['exitStatus'] = workflow.exitStatus
  email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
  email_fields['errorReport'] = (workflow.errorReport ?: 'None')
  email_fields['commandLine'] = workflow.commandLine
  email_fields['projectDir'] = workflow.projectDir
  email_fields['script_file'] = workflow.scriptFile
  email_fields['launchDir'] = workflow.launchDir
  email_fields['user'] = workflow.userName
  email_fields['Pipeline script hash ID'] = workflow.scriptId
  email_fields['Reference'] = REF
  email_fields['manifest'] = workflow.manifest
  email_fields['summary'] = summary

  email_info = ""
  for (s in email_fields) {
        email_info += "\n${s.key}: ${s.value}"
  }

  def output_d = new File( "${params.outdir}/pipeline_info/" )
  if( !output_d.exists() ) {
      output_d.mkdirs()
  }

  def output_tf = new File( output_d, "pipeline_report.txt" )
  output_tf.withWriter { w -> w << email_info }

 // make txt template
  def engine = new groovy.text.GStringTemplateEngine()

  def tf = new File("$baseDir/assets/email_template.txt")
  def txt_template = engine.createTemplate(tf).make(email_fields)
  def email_txt = txt_template.toString()

  // make email template
  def hf = new File("$baseDir/assets/email_template.html")
  def html_template = engine.createTemplate(hf).make(email_fields)
  def email_html = html_template.toString()

  def subject = "Virus Pipe run finished ($run_name)."

  if (params.email) {

        def mqc_report = null
        try {
                if (workflow.success && !params.skip_multiqc) {
                        mqc_report = multiqc_report.getVal()
                        if (mqc_report.getClass() == ArrayList){
                                log.warn "[IKMB VirusPipe] Found multiple reports from process 'multiqc', will use only one"
                                mqc_report = mqc_report[0]
                        }
                }
        } catch (all) {
                log.warn "[IKMB VirusPipe] Could not attach MultiQC report to summary email"
        }

        def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.maxMultiqcEmailFileSize.toBytes() ]
        def sf = new File("$baseDir/assets/sendmail_template.txt")
        def sendmail_template = engine.createTemplate(sf).make(smail_fields)
        def sendmail_html = sendmail_template.toString()

        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
        }

  }

}



