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
--pacbio 		       Pacbio subread movie file in BAM format
--email                        Email address to send reports to (enclosed in '')
Optional parameters:
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

REF = file("${baseDir}/assets/reference/MN908947.3.fasta")
REF_GFF = file("${baseDir}/assets/reference/MN908947.3.gff")
REF_NAME = "MN908947.3"

REF_WITH_HOST = file(params.ref_with_host)

/*
Primer sequences
*/

pb_barcodes = file("${baseDir}/assets/pacbio/Sequel_16_Barcodes_v3.fasta")
if (params.primers) {
        primers = file(params.primers)
} else {
        primers = file("${baseDir}/assets/primers/Eden_Sydney.fasta")
}

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
log.info "${workflow.manifest.description}	v${params.version}"
log.info "Nextflow Version:             $workflow.nextflow.version"
log.info "Viral reference:             	${REF}"
log.info "Host DB:			${params.bloomfilter_host}"
log.info "Virus DB:			${params.kraken2_db}"
log.info "Assemble de-novo:		${params.assemble}"
if (params.assemble) {
	log.info "Align assembly:			${params.align}"
}
log.info "Primer to trimming:		${primers}"
log.info "Command Line:			$workflow.commandLine"
if (workflow.containerEngine) {
        log.info "Container engine:		${workflow.containerEngine}"
}
log.info "==================================================================================="


// ********************
// WORKFLOW STARTS HERE
// ********************

Channel.fromPath(REF)
	.set { inputBloomMaker }

if (params.reads) {
	Channel.fromFilePairs(params.reads, flat: true)
	.set { reads_fastp }
} else {
	reads_fastp = Channel.empty()
}

if (params.pacbio) {
	Channel.fromPath(params.pacbio)
	.set { pb_reads }
} else {
	pb_reads = Channel.empty()
}

// ***************************
// PACBIO WORFKLOW
// ***************************

process runLima {

	input:
	path pb_read from pb_reads

	output:
	path "*.bam*" into pb_reads_demux
	path "*pbi" into pb_reads_demux_index
	script:
	
	prefix = pb_read.getBaseName() + ".demux"
	"""
		lima --num-threads ${task.cpus} --split-bam-named --same $pb_read $pb_barcodes $prefix 
	"""

}

// remove empty barcodes
pb_reads_demux.join(pb_reads_demux_index).filter { b,i -> b.size() > 200000 }.into { pb_reads_demux_valid_fasta; pb_reads_demux_valid_css }

process bam2fasta {

        input:
        set path(bam),path(pbi) from pb_reads_demux_valid_fasta

        output:
        path fasta into pb_reads_clean

        script:
        base_name = bam.getBaseName()
        fasta = base_name + ".fasta.gz"

        """
                bam2fasta -o $base_name $bam
        """

}

process runFlye {

        publishDir "${params.outdir}/flye/${lib}", mode: 'copy'

        input:
        path fasta from pb_reads_clean

        output:
        set val(lib),path(assembly) optional true into flye_assemblies
        path assembly_stats optional true

        script:
        lib = fasta.getBaseName()
        assembly = "outdir/assembly.fasta"
        assembly_stats = "outdir/assembly_info.txt"

        """
                 flye --pacbio-raw $fasta --genome-size 30k --threads ${task.cpus} --out-dir outdir 2>&1 || true
        """
}

process runCss {

	input:
	set path(bam),path(pbi) from pb_reads_demux_valid_css

	output:
	path(css) into css_reads_bam
	path(fasta) into css_reads_fasta

	script:
	base_name = bam.getBaseName()
	css = base_name + ".css.bam"
	fasta = css.getBaseName() + ".fasta.gz"

	"""
		css $bam $css
		bam2fasta -o $base_name $css
	"""
}

process alignPacbio {

	input:
	path fa_reads from css_reads_fasta

	output:
	set val(lib_name),path(bam),path(bai) into (pacbioBam_fb,pacbioBam)
	
	script:
	lib_name = fasta.getBaseName()
	bam = lib_name + ".align.bam"
	bai = bam + ".bai"

	"""
		minimap2 -t ${task.cpus} -ax asm20 $REF_WITH_HOST $fa_reads | samtoos -bh - | samtools sort -o $bam -
		samtools index $bam
	"""

}

process callVariantsPb {

	publishDir "${params.outdir}/pacbio", mode: 'copy'

	input:
	set val(lib_name),path(bam),path(bai) from pacbioBam_fb

	output:
	path vcf into pacbio_variants

	script:

	vcf = bam.getBaseName() . ".vcf"

	"""
		samtools index $bam
		freebayes --ploidy 1 -f $REF --genotype-qualities $bam > $vcf
	"""

}

// **********************
// ILLUMINA WORKFLOW
// **********************
process runFastp {

	label 'std'

	scratch true

        input:
        set val(sampleID), file(fastqR1),file(fastqR2) from reads_fastp

        output:
	set val(sampleID),file(left),file(right) into inputBowtie, inputBioBloomHost, inputBioBloomTarget, inputSpades
        set file(html),file(json) into fastp_qc

        script:

        left = fastqR1.getBaseName() + "_trimmed.fastq.gz"
        right = fastqR2.getBaseName() + "_trimmed.fastq.gz"
        json = fastqR1.getBaseName() + ".fastp.json"
        html = fastqR1.getBaseName() + ".fastp.html"


        """
                fastp --in1 $fastqR1 --in2 $fastqR2 --out1 $left --out2 $right --adapter_fasta $primers --detect_adapter_for_pe -w ${task.cpus} -j $json -h $html --length_required 35
        """
	
}

process runBloomMakerTarget {

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


// **********************
// Filter reads against a host contamination (e.g. human) - only reads not matching the host will survive
// **********************
process BloomfilterHost {

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
process runDeinterlave {

	label 'std'

        publishDir "${OUTDIR}/${id}/Bloomfilter/Host", mode: 'copy'

        input:
        set val(id),file(reads) from inputReformat

        output:
        set val(id),file(left),file(right) into ( inputPathoscopeMap, inputKraken )

        script:
        left = id + "_R1_001.bloom_non_host.fastq.gz"
        right = id + "_R2_001.bloom_non_host.fastq.gz"

        """
                reformat.sh in=$reads out1=$left out2=$right addslash int
        """

}

// ************************
// Get all the reads for a target species
// ************************
process BloomfilterTarget {

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
process runKraken2 {

	label 'kraken'

	publishDir "${OUTDIR}/${id}/Taxonomy/", mode: 'copy'

	when:
	params.taxonomy

	input:
	set val(id),file(left),file(right) from inputKraken

	output:
	set val(id),file(report) into KrakenReport
	file(kraken_log)

	script:
	report = id + ".kraken2_report.txt"
	kraken_log = id + ".kraken2.log"

	"""
		kraken2 --db $KRAKEN2_DB --threads ${task.cpus} --output $kraken_log --report $report $left $right 
	"""
}

process Kraken2Yaml {


	input:
	file(reports) from KrakenReport.collect()

	output:
	file(report_yaml) into KrakenYaml

	script:
	
	report_yaml = "kraken_report_mqc.yaml"
	"""
		kraken2yaml.pl --outfile $report_yaml
	"""

}

process runSpades {

	label 'spades'

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
		spades.py --meta --1 $left --2 $right $options -t ${task.cpus} -m ${task.memory.toGiga()} -o spades
	"""
}

// **************************
// Run assembly QC with QUAST
// **************************

contigsSpades.concat(flye_assemblies).into {assemblies; assemblies_qc }

process runQuast {

	label 'quast'

	publishDir "${OUTDIR}/${id}/Assembly/Spades/QC", mode: 'copy'

	input:
	set val(id),path(assembly) from assemblies_qc

	output:
	path(quast_dir) into QuastReport
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


// **********************
// Run pathoscope MAP
// **********************
process runPathoscopeMap {

	label 'pathoscope'

	when:
	params.taxonomy

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
process runPathoscopeId {

	label 'pathoscope'

	publishDir "${OUTDIR}/${id}/Pathoscope", mode: 'copy'

	when:
	params.taxonomy

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

// **************************
// align reads against reference genome
// **************************
process runBowtie {

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

all_bams = bowtieBam.concat(pacbioBam)

// ***********************
// Mark duplicate reads
// ***********************
process runMD {

       	label 'std'

	publishDir "${OUTDIR}/${id}/BAM", mode: 'copy'

	input:
	set val(id),file(bam),file(bai) from all_bams

	output:
	set val(id),file(bam_md_virus),file(bai_md_virus) into ( bamMD, inputBamStats, inputBamCoverage )
	set val(id),file(bam_md),file(bai_md)

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
process runCoverageStats {
	
	label 'std'

	publishDir "${OUTDIR}/${id}/BAM", mode: 'copy'

	input:
	set val(id),file(bam),file(bai) from inputBamStats.filter{ i,b,d -> b.size() > 10000 }

	output:
	file(global_dist) into BamStats
	file(report)

	script:
	global_dist = id + ".mosdepth.global.dist.txt"
	sam_coverage = id + ".coverage.samtools.txt"
	report = id + ".coverage.pdf"
	
	"""

		mosdepth -t ${task.cpus} $id $bam
		samtools depth -d 200  $bam > $sam_coverage
		bam2coverage_plot.R $sam_coverage $report
		
	"""
	
}

process runAlignStats {

        label 'std'

        publishDir "${OUTDIR}/${sampleID}/BAM", mode: 'copy'

        input:
        set val(sampleID),file(bam),file(bai) from inputBamCoverage

	output:
	file(align_stats) into BamAlignStats

	script:
	align_stats = sampleID + ".txt"

	"""
		samtools stats  $bam > $align_stats
	"""

}

// ************************
// Variant calling
// ************************
process runFreebayes {
	
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
		freebayes --ploidy 1 -f $REF --genotype-qualities $bam > $vcf
	"""

}

// ********************
// Filter variant calls
// ********************
process runFilterVcf {
	
       	label 'std'

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

// **********************
// Compile quality metrics into report
// **********************
process runMultiQC {

	label 'std'

	publishDir "${OUTDIR}/MultiQC", mode: 'copy'

	input:
	file('*') from fastp_qc.collect().ifEmpty('')
	file('*') from BamStats.collect().ifEmpty('')
	file('*') from BloomReportTarget.collect().ifEmpty('')
	file('*') from KrakenYaml.ifEmpty('')
	file('*') from BamAlignStats.collect()
	file('*') from QuastReport.collect().ifEmpty('')
	//file('*') from BloomReportHost.collect().ifEmpty('')

	output:
	file("multiqc_report.html") into multiqc_report

	script:
	def options = ""
	if (params.assemble) {
		options = "*/report.tsv"
	}
	"""
		cp $baseDir/assets/multiqc_config.yaml multiqc_config.yaml
		cp $params.logo .

		multiqc *.* $options
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



