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
--run_name			Specify a name for this analysis run
--reference			A reference genome in fasta format to compare against (default: NC_045512.2.fa)
--iva_guided			Perform reference-guided instead of de-novo (default: false - experimental and may crash)
--assemble			Assemble genomes de-novo
--email				Specify email to send report to
Output:
--outdir                       Local directory to which all output is written (default: results)
"""

params.help = false

// Show help when needed
if (params.help){
    log.info helpMessage
    exit 0
}

def summary = [:]

// Enable reference-based analyses like variant calling
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

// set basic global options like location of database files
BLOOMFILTER_HOST = params.bloomfilter_host

KRAKEN2_DB = params.kraken2_db

PATHOSCOPE_INDEX_DIR=file(params.pathoscope_index_dir)

OUTDIR = params.outdir

run_name = ( params.run_name == false) ? "${workflow.sessionId}" : "${params.run_name}"

if (params.run_name == false) {
        log.info "No run name was specified, using ${run_name} instead"
}

summary['Reference'] = REF
summary['Kraken2DB'] = params.kraken2_db
summary['PathoscopeDB'] = params.pathoscope_index_dir
summary['HostBloomFilter'] = params.bloomfilter_host


// Header log info
log.info "========================================="
log.info "${workflow.manifest.description} v${params.version}"
log.info "Nextflow Version:             $workflow.nextflow.version"
log.info "Reference:             		${REF}"
log.info "Command Line:                 $workflow.commandLine"
if (workflow.containerEngine) {
        log.info "Container engine:             ${workflow.containerEngine}"
}
log.info "========================================="


// ********************
// WORKFLOW STARTS HERE
// ********************

Channel.fromPath(REF)
	.set { inputBloomMaker }

Channel.fromFilePairs(params.reads, flat: true)
	.ifEmpty { exit 1, "Did not find any reads matching your criteria!" }
	.set { reads_fastp }

// **********************
// Perform read-trimming
// **********************
process runFastp {

	label 'std'

	scratch true

        input:
        set val(id), file(fastqR1),file(fastqR2) from reads_fastp

        output:
	set val(id),file(left),file(right) into (inputBioBloomHost , inputBioBloomTarget, inputBwa )
        set file(html),file(json) into fastp_qc

        script:

        left = fastqR1.getBaseName() + "_trimmed.fastq.gz"
        right = fastqR2.getBaseName() + "_trimmed.fastq.gz"
        json = fastqR1.getBaseName() + ".fastp.json"
        html = fastqR1.getBaseName() + ".fastp.html"

        """
                fastp --in1 $fastqR1 --in2 $fastqR2 --out1 $left --out2 $right --detect_adapter_for_pe -w ${task.cpus} -j $json -h $html --length_required 35
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
	set id,file(left_reads),file(right_reads) from inputBioBloomHost

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
        set id,file(clean_reads) into inputIva
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

	script:
	report = id + ".kraken2_report.txt"

	"""
		kraken2 --db $KRAKEN2_DB --threads ${task.cpus} --report $report $left $right 
	"""
}


// **********************
// Assemble reads, optionally with a reference
// **********************
process runIva {

        label 'std'

        publishDir "${OUTDIR}/${id}/Assembly/Iva", mode: 'copy'

	when:
	params.assemble

	input:
	set val(id),file(reads) from inputIva.filter{ i,r -> r.size() > 500000 }

	output:
	set val(id),file(contigs) into ( contigsIva, contigsAlign)

	script:

	outdir = "iva_" + id
	contigs = outdir + "/contigs.fasta"

	def options = ""
	if (params.iva_guided) {
		options = "--reference $REF"
	}

	println reads.size()
	"""
		iva --fr $reads -t ${task.cpus} $options $outdir
	"""

}

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


// **********************
// Align assembly against a reference genome
// **********************
process alignRef {

       	label 'std'

	publishDir "${OUTDIR}/${id}/Assembly/Iva/Alignment", mode: 'copy'
	
	when:
	params.align

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

// *********************
// Create mapping index from reference genome
// *********************
process makeBwaIndex {

       	label 'std'

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

// **************************
// align reads against reference genome
// **************************
process runBwa {

       	label 'std'

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
		bwa mem -H header.txt -M -R "@RG\\tID:${sampleID}\\tPL:ILLUMINA\\tSM:${sampleID}\\tLB:${sampleID}\\tDS:${REF}\\tCN:CCGA" -t ${task.cpus} $fasta $left $right | samtools sort -O bam -m 2G -@ 4 - > mapped.bam
		samtools view -b -o $bam -F 4 mapped.bam
		samtools index $bam
		rm mapped.bam
	"""
}

// ***********************
// Mark duplicate reads
// ***********************
process runMD {

       	label 'std'

	publishDir "${OUTDIR}/${id}/BAM", mode: 'copy'

	input:
	set val(id),file(bam),file(bai) from bwaBam

	output:
	set val(id),file(bam_md),file(bai_md) into ( bamMD, inputBamStats )

	script:
	bam_md = id + ".dedup.bam"
	bai_md = bam_md + ".bai"

	"""
		samtools sort -m 4G -t 2 -n $bam | samtools fixmate -m - fix.bam
		samtools sort -m 4G -t 2 -O BAM fix.bam | samtools markdup - $bam_md
		samtools index $bam_md
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
	set val(id),file(bam),file(bai) from inputBamStats

	output:
	set val(id),file(coverage_stats) into BamStats

	script:
	coverage_stats = id + ".wgs_coverage.txt"

	"""
		picard -Xmx${task.memory.toGiga()}G SortSam \
		I=$bam O=/dev/stdout SORT_ORDER=coordinate | picard -Xmx${task.memory.toGiga()}G CollectWgsMetrics \
                I=/dev/stdin \
       	        O=$coverage_stats \
               	R=$REF \
		MINIMUM_MAPPING_QUALITY=5
 
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
	//file('*') from BloomReportHost.collect().ifEmpty('')

	output:
	file("multiqc_report.html") into multiqc_report

	script:

	"""
		cp $baseDir/assets/multiqc_config.yaml multiqc_config.yaml
		cp $params.logo .

		multiqc *
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



