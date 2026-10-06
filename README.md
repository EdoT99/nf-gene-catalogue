# nf-gene-catalogue

A Nextflow pipeline that builds a **non-redundant gene catalogue of target enzymes** from metagenomic samples.

Each sample is processed **independently and in parallel**: reads are trimmed and assembled, genes are predicted, and the proteins are screened against user-provided HMM profiles. Only the proteins that pass the scoring filters, together with their gene sequences, are then pooled across samples and dereplicated into a non-redundant catalogue.

If you already have assemblies and Prodigal predictions (e.g. from a previous Geomosaic run), you can skip the read-based steps and start directly from them.

## Workflow

```
 reads samplesheet (sample, r1, r2) Table A
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
        │   ◄── or start here with a protein/contig samplesheet (sample, orfs, contigs) Table B
        ▼
 HMMSEARCH ──► all HMM profiles, fixed database size (-Z) for comparable E-values
        │
        ▼
 FILTER_HITS ──► E-value + bit score, best profile per ORF
        │
        ▼
 EXTRACT_GENES ──► nucleotide sequence of each hit, from the sample's own contigs
        │
        ▼
 POOL_HITS / POOL_GENES ──► hits of all samples, IDs prefixed with the sample name
        │
        ▼
 DEREPLICATION ──► non-redundant gene catalogue (MMseqs2)

```

| # | Step | Tool | Output | Runs |
|---|---|---|---|---|
| 1 | Trimming | fastp | trimmed R1 / R2, report | per sample |
| 2 | Assembly | MEGAHIT or metaSPAdes | contigs | per sample |
| 3 | Contig filtering & renaming | SeqKit | filtered contigs `contig_N` + name mapping | per sample |
| 4 | Gene prediction | Prodigal | proteins (`contig_N_M`) | per sample |
| 5 | Profile search | HMMER `hmmsearch` | hits table | per sample |
| 6 | Hit filtering | — | proteins passing E-value and bit score, best profile per ORF | per sample |
| 7 | Gene extraction | SeqKit | nucleotide sequences + coordinates of the hits | per sample |
| 8 | Pooling | — | one protein and one gene file, IDs prefixed with the sample name | once |
| 9 | Dereplication | MMseqs2 | non-redundant gene catalogue | once |

Starting from existing Prodigal output (entry point B) skips steps 1–4.

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

The pipeline has two entry points. Use **exactly one** of them; the pipeline stops if neither or both are given.

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

After assembly, contigs are filtered by length and renamed `contig_1`, `contig_2`, … (the same convention as [Geomosaic](https://github.com/giovannellilab/Geomosaic)); a mapping to the original assembler headers is saved for each sample. Prodigal names each protein after its contig plus a gene number (`contig_12_3` = gene 3 on `contig_12`).

During per-sample processing IDs stay as they are. When the hits are pooled, the sample name is prepended so IDs from different samples never collide:

```
protein:  >contig_12_3   →  >SAMPLE_A_contig_12_3
gene:     >contig_12_3   →  >SAMPLE_A_contig_12_3
```

Every catalogue entry can therefore be traced back to its sample, contig and gene.

### HMM profiles (`--hmm_db`)

A folder containing one or more profile files ending in `.hmm`. Defaults to the profiles shipped in `assets/hmm_db/`. All profiles are searched together in one `hmmsearch` per sample. The pipeline stops with an error if no profiles are found.

## Hit scoring

Each sample's proteins are filtered in three steps:

1. **E-value** ≤ `--filter_evalue` (full-sequence E-value)
2. **Bit score** ≥ `--min_bitscore`
3. **Best profile per ORF:** if an ORF matches several profiles, only the highest-scoring one is kept.

> **Why a fixed database size (`--hmm_z`)?** An HMMER E-value depends on how many sequences were searched. Because each sample is searched separately, a larger sample would otherwise get larger E-values than a smaller one for the very same hit. The pipeline therefore gives every search the same database size (`hmmsearch -Z`): by default the total number of ORFs over all samples, so E-values behave as if all samples had been searched together. Bit scores do not depend on database size and are directly comparable between samples.

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--input_sample_table` | – | Reads samplesheet (`sample,r1,r2`). Entry point A. |
| `--input_orf_table` | – | Protein/contig samplesheet (`sample,orfs,contigs`). Entry point B. |
| `--assembler` | `megahit` | `megahit` or `metaspades`. metaSPAdes needs considerably more memory. |
| `--hmm_db` | `assets/hmm_db` | Folder with `.hmm` profiles. |
| `--hmm_evalue` | `1e-5` | Loose reporting threshold inside `hmmsearch`. |
| `--hmm_z` | total ORFs | Database size used for E-values in every sample. Set it to skip counting ORFs on large datasets. |
| `--filter_evalue` | `1e-10` | E-value cutoff for keeping a hit. |
| `--min_bitscore` | `30` | Bit-score cutoff for keeping a hit. |
| `--batch` | `batch01` | Name used for the gene catalogue outputs. |

Parameters can be given on the command line or collected in a YAML file:

```yaml
# params.yaml
input_sample_table: samplesheet.csv
assembler: megahit
hmm_db: /data/hmm_profiles
filter_evalue: 1e-10
min_bitscore: 30
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

### Resuming

Add `-resume` to reuse steps that already finished. Because every sample is processed separately, this also works when samples are added: only the new ones are assembled and annotated.

Changing a scoring threshold reruns only the filtering and the steps after it, not assembly or HMMER:

```bash
nextflow run EdoT99/nf-gene-catalogue -r main -profile docker -params-file params.yaml \
    --min_bitscore 50 -resume
```

> Adding samples changes the default `--hmm_z` (total ORFs), which reruns `HMMSEARCH` for every sample. Set `--hmm_z` to a fixed value if you plan to add samples over time.

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

With the Slurm profile, Nextflow itself only coordinates: each task (one step of one sample) is submitted as its own Slurm job with the CPUs, memory and time set for that step, and runs in its own container or conda environment.

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
├── preprocessing/fastp/<sample>/          trimmed reads, fastp HTML/JSON reports          ┐
├── assembly/<assembler>/<sample>/                                                         │
│   ├── <sample>.final_contigs.fa              filtered contigs, renamed contig_N          │
│   ├── <sample>.contig_mapping.tsv            original header ⇥ new name                  │  one folder
│   └── <sample>.<assembler>.log               assembler log                               │  per sample
├── hmm_annotation/<sample>/                                                               │
│   ├── <sample>.tblout                        raw hmmsearch hits                          │
│   ├── <sample>.hmm_hits.tsv                  hits passing E-value + bit score, best per ORF
│   ├── <sample>.hits.faa                      protein sequences of the hits               │
│   ├── <sample>.hits.fna                      gene (nucleotide) sequences, same IDs       │
│   └── <sample>.hits_coords.tsv               contig, start, end, strand, ORF ID          ┘
├── hmm_annotation/all_samples_hmm_hits.tsv    filtered hits of all samples in one table
├── pooled_hits/
│   ├── proteins/                              hit proteins of all samples, IDs SAMPLE_contig_N_M
│   └── genes/                                 hit genes of all samples, same IDs
└── gene_catalog/<batch>/
    ├── <batch>_rep_seq.fasta                  non-redundant representative sequences
    ├── <batch>_cluster.tsv                    representative ⇥ member mapping
    ├── <batch>_all_seqs.fasta                 all sequences grouped by cluster
    └── <batch>_kept_ids.txt                   IDs of the representatives
```

`preprocessing/` and `assembly/` are only produced when starting from reads (entry point A).

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
