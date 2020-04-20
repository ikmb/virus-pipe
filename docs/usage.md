# Usage information

This pipeline can process Illumina short reads from sequenced viral cDNA, Metatranscriptomes or similar to:

* clean/trim reads
* remove reads mapping the a host genome (default: human)
* extract reads matching a reference genome (default: COVID19)
* - perform taxonomic identification against a database of viral reference genome
* - perform taxonomic identification against a database of clinical pathogens (bacterial and viral, 2015)
* map viral reads against a reference genome and perform variant calling (default: COVID19)
* assembly viral genome (optionally: guided against reference genome)
* - align de-novo assembly against reference genome

## Available options and settings


