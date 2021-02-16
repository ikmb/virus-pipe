# Output structure

This pipeline produces several outputs, roughly structured into "general pipeline infos", "per sample infos" and "joint infos". 

## Sample and run metrics

### `results/MultiQC`

A multi-sample summary report including information from all relevant processing steps. This is usually the first output to look at to gauge the overall quality of the data. 

### `results/pipeline_info`

This folder contains run metrics, like core hours used, the processing time for each step and a graph of the overall workflow. 

### `results/Reports`

Per-sample reports to be used for individual validation - contains information on assembly coverage, annotated SNPs and effect prediction. Available as human-readable PDF as well as JSON.

### `results/Pangolin`

The joint Pangolin lineage predictions.

### `results/RKI_Assemblis`

Per-sample genome assemblies for submission to the RKI. Each assembly is available as IUPAC-aware version and one, in which all variable bases are masked. Only the former is typically used. 

### `results/MySample`

Each sample has a number of specific outputs:

* BAM
* Bloomfilter
* CleanReads
* Pangolin  
* QC 
* RawReads  
* Taxonomy 
* Variants

The outputs should be largely self-explanatory. 
