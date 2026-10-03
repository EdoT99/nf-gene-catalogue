# nf-gene-catalogue

A Nextflow pipeline that builds a **non-redundant gene catalogue of target enzymes** from metagenomic samples.

Starting from paired-end reads, it trims and assembles each sample, predicts genes, pools the proteins and contigs of all samples, screens the proteins against user-provided HMM profiles, extracts the annotated proteins together with their gene sequences, and dereplicates them into a non-redundant set.

If you already have assemblies and predicted proteins (e.g. from a previous Geomosaic or MEGAHIT/Prodigal run), you can skip the read-based steps and start directly from them.

## Workflow

```
 reads samplesheet (sample, r1, r2)
        │
        ▼
 PREPROCESSING ──► fastp: adapter / quality trimming
        │
        ▼
 ASSEMBLY ──► MEGAHIT  or  metaSPAdes  (--assembler)
        │        └──► length filter + contigs renamed contig_1, contig_2, ...
        ▼
 ORF_PREDICTION ──► Prodigal (metagenomic mode)
        │
        │   ◄── or start here with an ORF samplesheet (sample, orfs, contigs)
        ▼
 POOL_ORFS / POOL_ASSEMBLIES ──► one protein and one contig file, IDs prefixed with the sample name
        │
        ├──► HMMSEARCH ──► one search per HMM profile
        │        │
        │        ▼
        │   FILTER_PROTEINS ──► annotated ORFs passing --filter_evalue
        │        │
        │        ▼
        │   EXTRACT_GENES ──► nucleotide sequence of each annotated gene
        │
        └──► DEREPLICATION ──► non-redundant gene catalogue (MMseqs2)
```

### Tools

| Step | Tool |
|---|---|
| Read trimming | [fastp](https://github.com/OpenGene/fastp) |
| Assembly | [MEGAHIT](https://github.com/voutcn/megahit) or [metaSPAdes](https://github.com/ablab/spades) |
| Contig filtering, sequence extraction | [SeqKit](https://bioinf.shenwei.me/seqkit/) |
| Gene prediction | [Prodigal](https://github.com/hyattpd/Prodigal) |
| Profile search | [HMMER](http://hmmer.org/) |
| Clustering / dereplication | [MMseqs2](https://github.com/soedinglab/MMseqs2) |

All tools are pulled automatically as containers (Docker / Singularity) or conda environments, one per step. You do not need to install them yourself.

## Requirements

**1. Nextflow** (≥ 25.10, requires Java 17+)

```bash
conda create -n nextflow -c conda-forge -c bioconda nextflow
conda activate nextflow
nextflow -version
```

**2. A software engine**: Docker on a workstation, Singularity/Apptainer or conda on HPC.

```bash
docker --version
docker run hello-world   # must work without sudo
```

If Docker is not installed, follow the [official installation guide](https://docs.docker.com/engine/install/). On Linux, add your user to the `docker` group so Docker runs without `sudo`:

```bash
sudo usermod -aG docker $USER   # then log out and back in
```

## Quick start: test run

Run the pipeline on the bundled mock dataset to check your setup:

```bash
nextflow run EdoT99/nf-gene-catalogue -r main -profile test,docker
```

## Input

The pipeline has two entry points. Use **one** of them.

### A. From reads (`--input_sample_table`)

A CSV file with one row per sample and its paired-end FASTQ files:

```csv
sample,r1,r2
SAMPLE_A,/data/SAMPLE_A_R1.fastq.gz,/data/SAMPLE_A_R2.fastq.gz
SAMPLE_B,/data/SAMPLE_B_R1.fastq.gz,/data/SAMPLE_B_R2.fastq.gz
```

| Column | Content |
|---|---|
| `sample` | Unique sample name, no spaces |
| `r1` | Forward reads (`.fastq.gz`) |
| `r2` | Reverse reads (`.fastq.gz`) |

The pipeline checks the header at launch and stops if a required column is missing.

### B. From existing assemblies and proteins (`--input_orf_table`)

If assembly and gene prediction were already done, skip trimming, assembly and Prodigal:

```csv
sample,orfs,contigs
SAMPLE_A,/data/SAMPLE_A/prodigal/orf_predicted.faa,/data/SAMPLE_A/megahit/geomosaic_contigs.fasta
SAMPLE_B,/data/SAMPLE_B/prodigal/orf_predicted.faa,/data/SAMPLE_B/megahit/geomosaic_contigs.fasta
```

| Column | Content |
|---|---|
| `sample` | Unique sample name, no spaces |
| `orfs` | Protein FASTA predicted by **Prodigal** (`.faa`) |
| `contigs` | The assembly the proteins were predicted from (plain or `.gz`) |

> **Note:** gene sequences are extracted using the coordinates stored in Prodigal's protein headers (`>ID # start # end # strand # ...`). The `orfs` files must therefore be unmodified Prodigal output, predicted from the `contigs` file in the same row.

### Contig and protein names

After assembly, contigs are filtered by length and renamed `contig_1`, `contig_2`, … (the same convention as [Geomosaic](https://github.com/giovannellilab/Geomosaic)); a mapping to the original assembler headers is saved for each sample. When pooling, the sample name is prepended to every contig **and** protein ID:

```
contig:   >contig_12      →  >SAMPLE_A_contig_12
protein:  >contig_12_3    →  >SAMPLE_A_contig_12_3   (gene 3 on that contig)
```

Every hit can therefore be traced back to its sample and contig.

### HMM profiles (`--hmm_db`)

A folder containing one or more profile files ending in `.hmm`. Defaults to the profiles shipped in `assets/hmm_db/`. The pipeline stops with an error if no profiles are found.

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--input_sample_table` | – | Reads samplesheet (`sample,r1,r2`). Entry point A. |
| `--input_orf_table` | – | Assembly/protein samplesheet (`sample,orfs,contigs`). Entry point B. |
| `--assembler` | `megahit` | `megahit` or `metaspades`. metaSPAdes needs considerably more memory. |
| `--hmm_db` | `assets/hmm_db` | Folder with `.hmm` profiles. |
| `--hmm_evalue` | `1e-5` | E-value reporting threshold for `hmmsearch`. |
| `--filter_evalue` | `1e-10` | Stricter E-value used to keep hits when filtering. |
| `--batch` | `batch01` | Name used for the gene catalogue outputs. |
| `--pooled_faa` | – | Existing pooled protein FASTA. Skips protein pooling. |
| `--pooled_contigs` | – | Existing pooled contig FASTA. Skips contig pooling. |

Parameters can be given on the command line or collected in a YAML file:

```yaml
# params.yaml
input_sample_table: samplesheet.csv
assembler: megahit
hmm_db: /data/hmm_profiles
filter_evalue: 1e-10
batch: batch01
```

```bash
nextflow run EdoT99/nf-gene-catalogue -r main -profile docker -params-file params.yaml
```

Command-line parameters override the params file, which overrides config files.

### Tool options

Extra options for individual tools are set in a config file through `ext.args`, without changing the pipeline:

```nextflow
// my_options.config
process {
    withName: 'MEGAHIT'               { ext.args = '--presets meta-sensitive' }
    withName: 'FILTER_RENAME_CONTIGS' { ext.args = '-m 1500' }          // minimum contig length (default 1000 bp)
    withName: 'DEREPLICATION'         { ext.args = '--min-seq-id 0.9 -c 0.8 --cov-mode 1' }
}
```

```bash
nextflow run EdoT99/nf-gene-catalogue -r main -profile docker -params-file params.yaml -c my_options.config
```

## Usage

### From reads

```bash
nextflow run EdoT99/nf-gene-catalogue -r main -profile docker \
    --input_sample_table samplesheet.csv \
    --assembler megahit \
    --hmm_db /path/to/hmm_profiles
```

### From existing assemblies and proteins

```bash
nextflow run EdoT99/nf-gene-catalogue -r main -profile docker \
    --input_orf_table orf_samplesheet.csv \
    --hmm_db /path/to/hmm_profiles
```

### From existing pooled files

```bash
nextflow run EdoT99/nf-gene-catalogue -r main -profile docker \
    --pooled_faa /path/to/pooled_orf_predicted.faa \
    --pooled_contigs /path/to/pooled_contigs.fa \
    --hmm_db /path/to/hmm_profiles
```

> **Important:** pooled files must use the same naming as the pipeline (see [Contig and protein names](#contig-and-protein-names)): removing the last `_<number>` from a protein ID must give the exact ID of its contig. If a contig cannot be found, `EXTRACT_GENES` stops with an error reporting how many genes were expected and how many were extracted. The safest option is to reuse the `pooled_orfs/` and `pooled_contigs/` files of a previous run.

### Resuming

Add `-resume` to reuse steps that already finished. For example, to try a different filtering threshold without rerunning assembly or HMMER:

```bash
nextflow run EdoT99/nf-gene-catalogue -r main -profile docker -params-file params.yaml \
    --filter_evalue 1e-20 -resume
```

## Profiles

Combine one profile for **where** to run with one for **which software** to use, separated by commas without spaces (e.g. `-profile univ_hpc,singularity`).

| Profile | Type | Description |
|---|---|---|
| `docker` | software | Run every tool in a Docker container |
| `singularity` | software | Run every tool in a Singularity/Apptainer container (HPC) |
| `conda` | software | Build a conda environment for every tool |
| `univ_hpc` | where | Submit every task as a Slurm job |
| `test` | data | Small mock dataset to check the installation |

Without a "where" profile, tasks run on the local machine.

## Running on an HPC

With the Slurm profile, Nextflow itself only coordinates: each task is submitted as its own Slurm job with the CPUs, memory and time set for that step, and runs in its own container or conda environment.

**1. Keep Nextflow running for the whole pipeline**, either in `tmux`/`screen` on the login node or as a small Slurm job:

```bash
#!/bin/bash
#SBATCH --job-name=nf-head
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=7-00:00:00

nextflow run EdoT99/nf-gene-catalogue -r main -profile univ_hpc,singularity \
    -params-file params.yaml -c my_cluster.config -resume
```

**2. Use a shared cache for images / environments**, so they are downloaded once and visible to all compute nodes. Add to `~/.bashrc` on the cluster (use a location on shared project or work storage, not `/tmp`):

```bash
export NXF_SINGULARITY_CACHEDIR=/path/to/shared/singularity_cache
export NXF_CONDA_CACHEDIR=/path/to/shared/conda_cache        # if you use conda
```

**3. Put cluster-specific settings in your own config** and pass it with `-c`:

```nextflow
// my_cluster.config
process.queue          = 'standard'            // your partition
process.clusterOptions = '--account=mygroup'   // if your cluster requires an account
process.resourceLimits = [ cpus: 36, memory: 750.GB, time: 7.d ]   // largest node / partition limit
```

Run the `work/` directory and results on shared storage that all nodes can read.

## Output

Results are written to `results/` (change it with `-output-dir`):

```
results/
├── preprocessing/fastp/<sample>/      trimmed reads, fastp HTML/JSON reports
├── assembly/<assembler>/<sample>/
│   ├── <sample>.final_contigs.fa          filtered contigs, renamed contig_N
│   ├── <sample>.contig_mapping.tsv        original header ⇥ new name
│   └── <sample>.<assembler>.log           assembler log
├── prodigal/<sample>/                 predicted proteins (.faa), genes (.fna), GFF
├── pooled_orfs/                       pooled_orf_predicted.faa
├── pooled_contigs/                    pooled_contigs.fa
├── hmmsearch/                         <profile>.tblout / .domtblout / .out, one set per profile
├── hmm_filtered/
│   ├── hmm_hits.tsv                   all ORF–profile hits passing --filter_evalue
│   ├── hmm_best_hits.tsv              best (lowest E-value) profile per ORF
│   ├── hmm_hits.faa                   protein sequences of the annotated ORFs
│   ├── hmm_hits.fna                   gene (nucleotide) sequences, same IDs
│   └── hmm_hits_coords.tsv            contig, start, end, strand, ORF ID (1-based)
└── gene_catalog/<batch>/
    ├── <batch>_rep_seq.fasta          non-redundant representative sequences
    ├── <batch>_cluster.tsv            representative ⇥ member mapping
    ├── <batch>_all_seqs.fasta         all sequences grouped by cluster
    └── <batch>_kept_ids.txt           IDs of the representatives
```

## Citations

If you use this pipeline, please cite the tools it relies on:

- **fastp**: Chen S, Zhou Y, Chen Y, Gu J. *fastp: an ultra-fast all-in-one FASTQ preprocessor.* Bioinformatics (2018).
- **MEGAHIT**: Li D, Liu CM, Luo R, Sadakane K, Lam TW. *MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph.* Bioinformatics (2015).
- **metaSPAdes**: Nurk S, Meleshko D, Korobeynikov A, Pevzner PA. *metaSPAdes: a new versatile metagenomic assembler.* Genome Research (2017).
- **Prodigal**: Hyatt D, et al. *Prodigal: prokaryotic gene recognition and translation initiation site identification.* BMC Bioinformatics (2010).
- **HMMER**: Eddy SR. *Accelerated profile HMM searches.* PLoS Computational Biology (2011).
- **MMseqs2**: Steinegger M, Söding J. *MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets.* Nature Biotechnology (2017).
- **SeqKit**: Shen W, et al. *SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation.* PLoS ONE (2016).
- **Nextflow**: Di Tommaso P, et al. *Nextflow enables reproducible computational workflows.* Nature Biotechnology (2017).

Contig filtering and renaming follow the conventions of [Geomosaic](https://github.com/giovannellilab/Geomosaic).

## License

Released under the [MIT License](LICENSE).
