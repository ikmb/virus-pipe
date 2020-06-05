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
--reference			A reference genome in fasta format to compare against (default: MN908947.3.fasta)
--gff				A GFF annotation for the custom reference
--kraken2_db			A kraken2-formatted database with virus species for taxonomic mapping
--assemble			Assemble genomes de-novo
--guided 			Guide assembly with known reference
--email				Specify email to send report to
--primers			Primers used for amplification of genome fragments
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

// Enable custom reference-based analyses like variant calling
if (params.reference) {
	REF = file(params.reference)
	if (!REF.exists() ) {
		exit 1, "Reference file not found!"
	}
	Channel.fromPath(params.reference)
		.ifEmpty { exit 1; "Could not find custom reference file (--reference)" }
		.set { FastaIndex }

	if (params.assemble && !params.gff) {
		exit 1, "Must provide a matching GFF file to your custom reference (--gff)"
	}
	if (params.gff) {
		REF_GFF = file(params.gff)
	}

} else {
	REF = file("${baseDir}/assets/reference/MN908947.3.fasta")
	REF_GFF = file("${baseDir}/assets/reference/MN908947.3.gff")
        Channel.fromPath("${baseDir}/assets/reference/MN908947.3.fasta")
       .set { FastaIndex }
}

if (!REF.exists() ) {
	exit 1, "Could not find reference file..."
}


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
	.into { inputBloomMaker; fasta_index }

if (params.reads) {
	Channel.fromFilePairs(params.reads, flat: true)
	.set { reads_fastp }
} else {
	reads_fastp = Channel.empty()
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
                fastp --in1 $fastqR1 --in2 $fastqR2 --out1 $left --out2 $right --detect_adapter_for_pe -w ${task.cpus} -j $json -h $html --length_required 35
        """
	
}

// *********************
// Create mapping index from reference genome
// *********************
process makeBowtieIndex {

       	label 'std'

	input:
	file(fasta) from FastaIndex

	output:
	set file(fasta),file(i1),file(i2),file(i3),file(i4),file(r1),file(r2) into BowtieIndex	

	script:
	base_name = fasta.getBaseName()

	i1 = base_name + ".1.bt2"
	i2 = base_name + ".2.bt2"
	i3 = base_name + ".3.bt2"
	i4 = base_name + ".4.bt2"
	r1 = base_name + ".rev.1.bt2"
	r2 = base_name + ".rev.2.bt2"
		
	"""		
		bowtie2-build $fasta $base_name
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
	set file(fasta),file(i1),file(i2),file(i3),file(i4),file(r1),file(r2) from BowtieIndex.collect()

	output:
	set val(sampleID),file(bam),file(bai) into bowtieBam

	script:
	ref_name = fasta.getBaseName()
	bam = sampleID + "." + ref_name + ".aligned.bam"
	bai = bam + ".bai"

	"""
		bowtie2 --rg-id ${sampleID} --rg PL:Illumina --rg SM:${sampleID} --rg CN:CCGA -x $ref_name -p ${task.cpus} --no-unal --sensitive -1 $left -2 $right | samtools sort - | samtools view -bh -o $bam - 
		samtools index $bam
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
        set val(id),file(left),file(right) from inputSpades.filter{ i,r -> r.size() > 500000 }

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

