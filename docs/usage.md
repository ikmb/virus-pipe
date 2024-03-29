# Usage information

This pipeline can process Illumina short reads from sequenced viral cDNA, Metatranscriptomes or similar to:

* clean/trim reads
* align reads against a defined reference (virus or virus+human)
* variant calling against viral reference + normalization + filtering
* create RKI-compliant consensus assembly (ref-flipping)
* lineage-typing using Pangolin

Optional/non-essential additional steps:

* remove reads mapping the a host genome (default: human)
* extract reads matching a reference genome (default: COVID19)
* - perform taxonomic identification against a database of viral reference genome
* assemble viral genome de-novo
* - scaffold de-novo contigs against reference
* - align de-novo assembly against reference genome

## Basic execution

Assuming an config file exists (default = Medcluster of the IKMB), the following commend will start the pipeline (nextflow and singularity must be available in the PATH):

`nextflow run ikmb/virus-pipe --reads '/path/to/*_R{1,2}_001.fastq.gz'`

For more details on available options, see below.

## Available options and settings

### `--folder` 
Path to a folder containing Illumina PE reads. Will automatically group files into libraries and lanes. Assumes the CCGA (standard Illumina´) naming scheme: `*_L0*_R{1,2}_001.fastq.gz`. Mutually exclusive with `--reads` and `--samples`.

### `--reads` 
Regexp pointing to a list of PE Illumina reads for analysis (e.g. --reads /path/to/*_L0*_R{1,2}_001.fastq.gz). Must be enclosed in single-quotes. Mutually exclusive with `--folder` and `--samples`.

### `--samples`
Path to a CSV formatted sample sheet as an alternative to --reads. Expects the following columns:

```
IndivID;SampleID;R1;R2
21Ord1339;21Ord1339-L1;/path/to/reads_R1_001.fastq.gz;/path/to/reads_R2_001.fastq.gz
```

A script is included with this code base to produce such a file from a folder of fastQ files

Mutually exclusive with `--reads` and `--folder`. 

```
ruby /path/to/samplesheet_from_folder.rb -f /path/to/folder > Samples.csv
```

You can edit this file to replace the library-derived labels for IndivID and/or SampleID with an order number or patient ID. 

### `--run_name`
Provide a usefull name to this analysis run (could be the LIMS project ID)

### `--outdir` (default: results)
A folder name to which all outputs are written. This can also be a path (relative or absolute) - the target location with be created if necessary. 

### `--clip` (default: 20)
Remove n bases from both 3' and 5' of each read to account for fragmented amplicon primers that cannot be detected otherwise.

### `--primer_set` (default: ARTIC-v3)
Defines which set of PCR primers was used if this is sequenced from amplicons using Illumina (ARTIC-v3 or Eden). This option is likely not useful if the library prep fragments/nicks the adapters and/or reads are generated from fragmented amplicons. If this is the case, the option `--clip` should be used to simply trim off bases from the ends of all reads to remove potential PCR primer contamination. 

### `--primer_fasta` (default: false)
Provide a set of primer sequences in FASTA format (overrides --primer_set option)

### `--var_call_cov` (default: 10)
Minimum coverage at a site required for analysis

### `--var_call_count` (default: 10)
Minimum number of reads required to support a SNP call

### `---var_call_frac` (default: 0.1)
Minimum fraction of reads required to call a SNP

### `--var_filter_mqm` (default: 40)
Minimum mean mapping quality required for a variant to pass

### `--var_filter_sap` (default: 2000)
Maximum strand bias probability to include a variant call (note: this value may not be useful with amplicon data - hence the high default!)

### `--var_filer_qual` (default: 10)
Mimum call quality to include a variant

### `--assemble` (default: false)
Assemble the reads using Spades. This option is always on by default.

### `--guided` (default: false)
Use the built-in Covid19 reference to guide assembly (this will likely inflate assembly metrics). This option is always on by default. 

### `--filter` (default: true)
If this option is set, the trimmed reads will be mapped and filtered against the human genome using Bowtie2 prior to taxonomic assignment.

## Expert options (safe to ignore!)

The below settings are typically set in a site-specific config file and not manipulated from the command line. 

### `--metadata`
Enables an additional process that retrieves sample-specific meta data to produce a RKI compatible meta data sheet. This uses an in-house script and only works on the IKMB infrastructure. It is enabled by default on systems that support this option. 

### `--kraken2_db` 
The location of a Kraken2 formatted DB with viral sequences

If you want to set up your own database, please see: https://github.com/DerrickWood/kraken2/wiki/Manual#custom-databases

### `--host_index`
A Bowtie2 formatted index of the human genome (e.g. from iGenomes) - is used to remove host reads prior to taxonomic assignment and (optional) de-novo assembly.

### `--ref_with_host`
The base name of a BWA formatted index containing a human genome sequence together with the Sars-CoV2 [genome](../assets/reference/NC_045512.2.fa). This is used to have decoys for mapping and variant calling. If this is not set, the built-in viral reference genome is used only.  
This should be fine, we have not yet seen a clear advantage of using a combined reference when working off amplicon data (this might be different for metagenomic data). 

### `--pango_data`
A local version of the pangolin database; can overwrite whatever version is shipped with the pangolin container. 
