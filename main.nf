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
--primer_set			Name of the primer set used (ARTIC-v3, Eden)
--clip				Remove x bases from both ends of the reads (default: 6)
--primer_fasta			Primer sequences in FASTA format (overrides --primer_set)
--run_name			Specify a name for this analysis run
--kraken2_db			A kraken2-formatted database with virus species for taxonomic mapping
--assemble			Assemble genomes de-novo
--guided 			Guide assembly with known reference
--email				Specify email to send report to
--primers			Primers used for amplification of genome fragments (default: Eden)
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

/*
Primer sequences
*/

pb_barcodes = file("${baseDir}/assets/pacbio/Sequel_16_Barcodes_v3.fasta")

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
host_genome = Channel.fromPath("${params.host_index}*")

// set basic global options like location of database files
BLOOMFILTER_HOST = params.bloomfilter_host

KRAKEN2_DB = params.kraken2_db

PATHOSCOPE_INDEX_DIR=file(params.pathoscope_index_dir)

OUTDIR = params.outdir

run_name = ( !params.run_name) ? "${workflow.sessionId}" : "${params.run_name}"

if (!params.run_name ) {
        log.info "No run name was specified, using ${run_name} instead"
}

summary['Reference'] = REF
summary['Kraken2DB'] = params.kraken2_db
summary['PathoscopeDB'] = params.pathoscope_index_dir
summary['HostBloomFilter'] = params.bloomfilter_host

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
log.info "Host DB:			${params.bloomfilter_host}"
log.info "Virus DB:			${params.kraken2_db}"
log.info "Assemble de-novo:		${params.assemble}"
if (params.assemble) {
	log.info "Assemble guided:		${params.guided}"
}
if (params.primers) {
	log.info "Primers for trimming:		${primers}"
} else if (params.primer_set) {
	log.info "Primers for trimming:		${params.primer_set}"
} else {
	log.info "No designated primers for trimming - only Illumina primers used"
}
log.info "Read clipping 3'/5'		${params.clip}"
log.info "Command Line:			$workflow.commandLine"
if (workflow.containerEngine) {
        log.info "Container engine:		${workflow.containerEngine}"
}
log.info "==================================================================================="


// ********************
// WORKFLOW STARTS HERE
// ********************

Channel.fromPath(REF)
	.into { inputBloomMaker; inputNormalize  }

if (params.reads) {
	Channel.fromFilePairs(params.reads, flat: true)
	.set { reads_fastp }
} else {
	reads_fastp = Channel.empty()
}

// **********************
// ILLUMINA WORKFLOW
// **********************

process trim_reads {

	label 'fastp'

	scratch true

        publishDir "${OUTDIR}/${sampleID}/RawReads", mode: 'copy'

        input:
        set val(sampleID), file(fastqR1),file(fastqR2) from reads_fastp
	file(primer_fa) from primers.collect().ifEmpty(false)

        output:
	set val(sampleID),file(left),file(right) into (inputBowtie, inputBioBloomHost, inputBioBloomTarget, inputBowtieFilter )
        set file(html),file(json) into fastp_qc

        script:
	def options = ""
	if (primer_fa) {
		options = "--adapter_fasta " + primer_fa
	}

        left = fastqR1.getBaseName() + "_trimmed.fastq.gz"
        right = fastqR2.getBaseName() + "_trimmed.fastq.gz"
        json = fastqR1.getBaseName() + ".fastp.json"
        html = fastqR1.getBaseName() + ".fastp.html"

        """
                fastp --in1 $fastqR1 --in2 $fastqR2 --out1 $left --out2 $right $options -f $params.clip -t $params.clip --detect_adapter_for_pe -w ${task.cpus} -j $json -h $html --length_required 35
        """
	
}

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

if (params.fast_filter) {

	// **********************
	// Using Bloom filter
	// **********************

	process remove_host_reads_bloom {

        	label 'std'

		publishDir "${OUTDIR}/${id}/Bloomfilter/Host", mode: 'copy'

		input:
		set val(id),file(left_reads),file(right_reads) from inputBioBloomHost

		output:
		set id,file(clean_reads) into inputReformat
		file(bloom) into BloomReportHost

		script:
		analysis = id + ".Host" 
		bloom = analysis + "_summary.tsv"
		clean_reads = analysis + ".filtered.fastq.gz"

		"""
			biobloomcategorizer -p $analysis --gz_output -d -n -e -s 0.01 -t ${task.cpus} -f "$BLOOMFILTER_HOST" $left_reads $right_reads | gzip > $clean_reads
		"""

	}

	// *************************
	// Take the interlaved non-host reads and produce sane PE data 
	// *************************
	process deinterleave_reads {

		label 'std'

	        publishDir "${OUTDIR}/${id}//CleanReads", mode: 'copy'

        	input:
	        set val(id),file(reads) from inputReformat

        	output:
	        set val(id),file(left),file(right) into ( inputPathoscopeMap, inputKraken, inputSpades )

        	script:
	        left = id + "_R1_001.bloom_non_host.fastq.gz"
        	right = id + "_R2_001.bloom_non_host.fastq.gz"
	
        	"""
	                reformat.sh in=$reads out1=$left out2=$right addslash int
        	"""

	}

} else if (params.filter) {

        // ****************
        // Using Bowtie2
        // ****************

        process remove_host_reads_bt {

                label 'std'

                publishDir "${OUTDIR}/${id}/CleanReads", mode: 'copy'

                input:
                set val(id),file(left_reads),file(right_reads) from inputBowtieFilter
                path(bwt_files) from host_genome.collect()

                output:
                set val(id),file(left_clean),file(right_clean) into ( inputPathoscopeMap, inputKraken, inputSpades )

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
	inputBioBloomHost.into { inputKraken; inputPathoscopeMap; inputSpades  }
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
	set val(id),file(left),file(right) from inputKraken

	output:
	set val(id),file(report) into (KrakenReport, Kraken2Report)
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

process assemble_virus {

	label 'std'

	publishDir "${OUTDIR}/${id}/Assembly/Spades", mode: 'copy'

	when:
	params.assemble

	input:
	set val(id),file(left),file(right) from inputSpades

	output:
	set val(id),file(scaffolds) into contigsSpades
	
	script:
	scaffolds = "spades/scaffolds.fasta"

	def options = ""
	if (params.guided) {
		options = "--trusted-contigs ${REF}"
	}

	"""
		coronaspades.py -1 $left -2 $right $options -t ${task.cpus} -m ${task.memory.toGiga()} -o spades
	"""
}

// **************************
// Run assembly QC with QUAST
// **************************

contigsSpades.into {assemblies; assemblies_qc; assemblies_pangolin }

process assembly_pangolin {

	label 'pangolin'

	publishDir "${OUTDIR}/${id}/Pangolin", mode: 'copy'

	input:
	set val(id),path(assembly) from assemblies_pangolin

	output:
	set val(id),path(report) into (pangolin_report, Pangolin2Report)

	script:

	report = id + ".pangolin.csv"

	"""
		pangolin --outfile $report $assembly 
	"""
}

process pangolin2yaml {

	label 'std'

	publishDir "${OUTDIR}/Pangolin", mode: 'copy'

	input:
	file(reports) from pangolin_report.collect()

	output:
	file(report) into PangolinYaml

	script:

	report = "pangolin_report_mqc.yaml"

	"""
		pangolin2yaml.pl > $report
	"""

}

process assembly_qc {

	label 'quast'

	publishDir "${OUTDIR}/${id}/Assembly/Spades/QC", mode: 'copy'

	input:
	set val(id),path(assembly) from assemblies_qc

	output:
	path(quast_dir) into QuastReport
	set val(id),file("${quast_dir}/report.tsv") into Quast2Report
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
		set id,file(left_reads),file(right_reads) from inputPathoscopeMap

		output:
		set id,file(pathoscope_sam) into inputPathoscopeId

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
		set id,file(samfile) from inputPathoscopeId.filter{ i,r -> r.size() > 1000000 }

		output:
		set id,file(pathoscope_tsv) into outputPathoscopeId

		script:

		//pathoscope_sam = "updated_" + samfile
		pathoscope_tsv = id + "-sam-report.tsv"

		"""
        		pathoscope ID -alignFile $samfile -fileType sam -expTag $id
		"""

	}

} // end pathoscope

// **************************
// align reads against reference genome
// **************************
process align_viral_reads {

       	label 'std'

	publishDir "${OUTDIR}/${sampleID}/BAM/raw", mode: 'copy'

	input:
	set val(sampleID),file(left),file(right) from inputBowtie

	output:
	set val(sampleID),file(bam),file(bai) into bowtieBam

	script:
	ref_name = REF_WITH_HOST.getBaseName()
	bam = sampleID + "." + ref_name + ".aligned.bam"
	bai = bam + ".bai"

	"""
		bowtie2 --rg-id ${sampleID} --rg PL:Illumina --rg SM:${sampleID} --rg CN:CCGA -x $REF_WITH_HOST -p ${task.cpus} --no-unal --sensitive -1 $left -2 $right | samtools sort - | samtools view -bh -o $bam - 
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
	set val(id),file(bam),file(bai) from bowtieBam

	output:
	set val(id),file(bam_md_virus),file(bai_md_virus) into ( bamMD, inputBamStats, inputBamCoverage )
	set val(id),file(bam_md),file(bai_md) into HostBam

	script:
	bam_md = id + ".dedup.bam"
	bai_md = bam_md + ".bai"
	bam_md_virus = id + ".dedup.viral_mapping.bam"
	bai_md_virus = bam_md_virus + ".bai"

	"""
		samtools sort -m 4G -t 2 -n $bam | samtools fixmate -m - fix.bam
		samtools sort -m 4G -t 2 -O BAM fix.bam | samtools markdup - $bam_md
		samtools index $bam_md

		samtools view -bh -o $bam_md_virus $bam_md $REF_NAME
		samtools index $bam_md_virus
		rm fix.bam 
	"""

}

// ************************
// Coverage statistics
// ************************
process bam2bed {

	label 'bedtools'

	publishDir "${OUTDIR}/${id}/BAM", mode: 'copy'

	input:
	set val(id),file(bam_md),file(bai_md) from HostBam.filter { b -> b.size() > 8000 }

	output:
	file(bed) 

	script:
	
	bed = bam_md.getBaseName() + ".Host.bed"

	"""
		bedtools bamtobed -i $bam_md | grep ^chr > $bed
	"""

}

process coverage_stats {
	
	label 'std'

	publishDir "${OUTDIR}/${id}/BAM", mode: 'copy'

	input:
	set val(id),file(bam),file(bai) from inputBamStats.filter{ i,b,d -> b.size() > 10000 }

	output:
	file(global_dist) into BamStats
	set val(id),file(sam_coverage) into BamCoverage
	file(report)

	script:
	global_dist = id + ".mosdepth.global.dist.txt"
	sam_coverage = id + ".coverage.samtools.txt"
	report = id + ".coverage.pdf"
	
	"""
		mosdepth -t ${task.cpus} $id $bam
		samtools depth -d 200  $bam > $sam_coverage
		bam2coverage_plot.R $sam_coverage ${params.cov_lim} $report
	"""
	
}

process align_stats {

        label 'std'

        publishDir "${OUTDIR}/${sampleID}/BAM", mode: 'copy'

        input:
        set val(sampleID),file(bam),file(bai) from inputBamCoverage

	output:
	file(align_stats) into BamAlignStats
	set val(sampleID),file(align_stats) into Samtools2Report

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

	publishDir "${OUTDIR}/${id}/Variants", mode: 'copy'
	input:
	set val(id),file(bam),file(bai) from bamMD

	output:
	set val(id),file(vcf) into fbVcf

	script:
	base_name = bam.getBaseName()
	vcf = base_name + ".vcf"

	"""
		freebayes --genotype-qualities --min-coverage 10 -V --ploidy 1 -f $REF --genotype-qualities $bam > $vcf
	"""

}

// ********************
// Filter variant calls
// ********************
process filter_vcf {
	
       	label 'std'

	//publishDir "${OUTDIR}/${id}/Variants", mode: 'copy'

	input:
	set val(id),file(vcf) from fbVcf

	output:
	set val(id),file(vcf_filtered) into VcfNormalize

	script:
	vcf_filtered = vcf.getBaseName() + ".filtered.vcf"

	"""
		bcftools filter ${params.filter_options} $vcf > $vcf_filtered
	"""
}

process normalize_vcf {

	label 'std'

	input:
	set val(id),file(vcf) from VcfNormalize
        file(ref_genome) from inputNormalize.collect()

	output:
        set val(id),file(vcf_filtered) into (finalVcf,Vcf2Report,VcfPredict)

	script:

	vcf_filtered = vcf.getBaseName() + ".normalized.vcf"

	"""
		vt normalize -o $vcf_filtered -r $ref_genome $vcf
	"""

}

process vcf_stats {

	label 'std'

	publishDir "${OUTDIR}/${id}/Variants", mode: 'copy'

	input:
	set val(id),file(vcf) from finalVcf

	output:
	file(stats) into VcfStats
	set val(id),file(stats) into VcfReport

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
	set val(id),file(vcf) from VcfPredict

	output:
	set val(id),file(effects) into EffectPrediction

	script:

	effects = id + ".snpeff." + REF_NAME + ".vcf"

	"""
		snpEff -v $REF_NAME -onlyProtein -no-upstream -no-downstream -canon $vcf > $effects
	"""

}

// **********************
// Write a per-patient report
// **********************

GroupedReports = Kraken2Report.join(Pangolin2Report).join(Samtools2Report).join(Quast2Report).join(EffectPrediction)


process final_report {

	label 'std'

	publishDir "${OUTDIR}/${id}/Report", mode: 'copy'

	input:
	set val(id),file(kraken),file(pangolin),file(samtools),file(quast),file(variants) from GroupedReports

	output:
	file(patient_report) 

	script:

	patient_report = id + "_report.pdf"

	"""
		cp $baseDir/assets/ikmb_bfx_logo.jpg . 
		covid_report.pl --kraken $kraken --pangolin $pangolin --bam_stats $samtools --assembly_stats $quast --vcf $variants --outfile $patient_report
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
	//file('*') from BloomReportHost.collect().ifEmpty('')

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



