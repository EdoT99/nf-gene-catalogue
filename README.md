## This repo contains a NF (NExtflow) based pipeline for building a Gene catalogue for Nickel & Iron binidnign enzymes from metagenomcs samples
###


### Set up
Install Nextflow
```bash
conda create -n nextflow -c conda-forge -c bioconda nextflow
conda activate nextflow
```
Make sure you have docker or docker engine on you machine
```bash

```

### Testing pipeline with a mock dataset
```bash
nextflow run EdoT99/nf-gene-catalogue -r v1.0.0 -profile test,docker
```