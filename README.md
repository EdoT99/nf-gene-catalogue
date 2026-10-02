## Is a Nextflow based pipeline for building a Gene catalogue of enzymes from metagenomcs samples
###
The pipeline pools together all predicted orf (prodigal output), it searches a series of HMM profiles passed by the user for the proteins of interest, it dereplicates the HITs and build a Gene-catalogue of unique representative sequences from several environments

### Set up
Install Nextflow
```bash
conda create -n nextflow -c conda-forge -c bioconda nextflow
conda activate nextflow
```
Make sure you have docker or docker engine on you machine
```bash
docker --version
whcih docker
```
if not, follow the installation link for your specifc OS: [https://docs.docker.com/engine/install/](docker_installation)

### Testing pipeline with a mock dataset
```bash
nextflow run EdoT99/nf-gene-catalogue -r v1.0.0 -profile test,docker
```

### Tips on execution
You can edit and choose your preferred configurations using the nextflow.config file, the profiles enable to seelct multipl configuraionts on machine and type of
Find more info about how to set up a config at: [https://training.nextflow.io/latest/hello_nextflow/transcripts/06_hello_config/](nextflow_confi)
#### End-to-end
```bash
nextflow run main.nf -profile docker
```
#### Intermediate-to-end
```bash
nextflow run main.nf -profile docker --pooled_faa /path/to/pooled_orf.file
```
