# Usage information

This pipeline can process Illumina short reads from sequenced viral cDNA, Metatranscriptomes or similar to:

* clean/trim reads
* remove reads mapping the a host genome (default: human)
* extract reads matching a reference genome (default: COVID19)
* - perform taxonomic identification against a database of viral reference genome
* - perform taxonomic identification against a database of clinical pathogens (bacterial and viral, 2015)
* map viral reads against a reference genome and perform variant calling (default: COVID19)
* assembly viral genome
* - align de-novo assembly against reference genome

## Basic execution

Assuming an config file exists (default = Medcluster of the IKMB), the following commend will start the pipeline (nextflow and singularity must be available in the PATH):

`nextflow run ikmb/virus-pipe --reads '/path/to/*_R{1,2}_001.fastq.gz'`

For more details on available options, see below.

## Available options and settings

### `--reads` 
Path to a folder with PE Illumina reads for analysis (e.g. --reads /path/to/*_R{1,2}_001.fastq.gz)

### `--clip` (default: 6)
Remove n bases from both 3' and 5' to account for fragmented amplicon primers that cannot be detected otherwise.

### `--primer_set` (default: ARTIC-v3)
Defines which set of PCR primers was used if this is sequenced from amplicons using Illumina (ARTIC-v3 or Eden). 

### `--primer_fasta` (default: false)
Provide a set of primer sequences in FASTA format (overrides --primer_set option)

### `--assembly` (default: true)
Assemble the reads using Spades. 

### `--guided` (default: true)
Use the built-in Covid19 reference to guide assembly (this will likely inflate assembly metrics)

### `--taxonomy` (default: false)
Run the Pathoscope pipeline to check for other metagenomic signals in the data.

### `--filter` (default: true)
If this option is set, the trimmed reads will be run against a bloom filter to remove likely human reads prior to taxnomic assignment.

### `--fast_filter` (default: false)
Use a rapid bloom filter to remove host reads. This is significantly faster than an exhaustive alignment (--filter), but is less testedt. Use only when processing high sample numbers and/or very deep samples. 

## Expert options (safe to ignore!)

The below settings are typically set in a site-specific config file and not manipulated from the command line. 

### `--kraken2_db` 
The location of a Kraken2 formatted DB with viral sequences

### `--pathoscope_index_dir`
The location of a folder containing index files for Pathoscope (last available reference data from 2015 assumed)

### `--ref_with_host`
The base name of a bowtie2 formatted index containing a human genome sequence together with the Covid19 genome


