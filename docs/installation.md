# Installation instructions

This pipeline is configured to run on the IKMB MedCluster. If you are using this system, nothing needs to be done. 

IF you want to try and run this pipeline on your own system, here are some basic pointers. Feel free to open an issue so 
I know what there is interest to have this pipeline a bit more portable. 

All the instructions below assume that you are familiar with Nextflow. 

## Nextflow config file

Your site-specific config file should, in addition to the cluster parameters, contain the following information:

```
params.kraken2_db = "/work_ifs/ikmb_repository/databases/Kraken2/2020-03_viruses"
params.host_index = "/work_ifs/ikmb_repository/references/iGenomes/references/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/genome"
```

The option kraken2_db should point to a Kraken2 Virus database - make sure that this included the Sars-CoV2 genome. When in doubt, you can pass the genome sequence after the build finishes. 
Instructions can be found [here](https://github.com/DerrickWood/kraken2/wiki/Manual#custom-databases)

The option host_index points to a Bowtie2 index of the human genome. We use this to extract reads that do NOT map to human for taxonomic profiling with Kraken as well as the optional de-novo assembly of the viral genome. 
If you are using e.g. iGenomes, you can easily just use that (or make your own index). 


```
params.ref_with_host="/work_ifs/ikmb_repository/databases/custom_indices/Homo_sapiens_GRCh38_no_alts_with_virus.fa"
```

This setting is purely optional and points to the base of a BWA index containing both the human genome as well as the Sars-CoV2 Referenz (from this code base under assets/references). We use this as a kind of decoy for mapping and variant calling.
But it is very likely that you can achieve the same results just using the built-in reference only. 

