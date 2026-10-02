# nf-gene-catalogue

A Nextflow pipeline that builds a **non-redundant gene catalogue of target enzymes** from metagenomic samples.

It pools the predicted proteins (ORFs) of many samples, screens them against user-provided HMM profiles, extracts the annotated proteins, and dereplicates them into a non-redundant set.

## Workflow

```
samplesheet (sample, orfs)
        │
        ▼
   POOL_ORFS ──────────────► pooled ORFs, headers prefixed with sample name
        │
        ├──► HMMSEARCH ────► one search per HMM profile (HMMER)
        │        │
        │        ▼
        │   FILTER_PROTEINS ► hits passing the E-value threshold + their sequences (SeqKit)
        │
        └──► MMSEQ ─────────► non-redundant gene catalogue (MMseqs2)
```

Third-party tools used:

| Tool | Purpose |
|---|---|
| [HMMER](http://hmmer.org/) | screening proteins against HMM profiles |
| [MMseqs2](https://github.com/soedinglab/MMseqs2) | clustering / dereplication |
| [SeqKit](https://bioinf.shenwei.me/seqkit/) | filtering and extracting sequences |

All tools are pulled automatically as containers (or conda environments). You do not need to install them yourself.

## Requirements

**1. Nextflow** (≥ 25.10, requires Java 17+)

```bash
conda create -n nextflow -c conda-forge -c bioconda nextflow
conda activate nextflow
nextflow -version
```

**2. A container engine**: Docker on a workstation, or Singularity/Apptainer on HPC.

```bash
docker --version
docker run hello-world   # must work without sudo
```

If Docker is not installed, follow the [official installation guide](https://docs.docker.com/engine/install/) for your OS. On Linux, add your user to the `docker` group so Docker runs without `sudo`:

```bash
sudo usermod -aG docker $USER   # then log out and back in
```

## Quick start: test run

Run the pipeline on the bundled mock dataset to check your setup:

```bash
nextflow run EdoT99/nf-gene-catalogue -r v1.0.0 -profile test,docker
```

## Input

### Samplesheet (`--input`)

A CSV file with one row per sample, pointing to its predicted proteins (e.g. Prodigal `.faa` output):

```csv
sample,orfs
SAMPLE_A,/data/SAMPLE_A/prodigal/orf_predicted.faa
SAMPLE_B,/data/SAMPLE_B/prodigal/orf_predicted.faa
```

Sample names must be unique and must not contain spaces. They are prepended to every protein ID (`>k141_1_1` → `>SAMPLE_A_k141_1_1`), so hits can always be traced back to their sample.

### HMM profiles (`--hmm_db`)

A folder containing one or more profile files ending in `.hmm`. The pipeline stops with an error if no profiles are found.

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--input` | – | Samplesheet CSV (`sample,orfs`). Not needed if `--pooled_faa` is given. |
| `--hmm_db` | – | **Required.** Folder with `.hmm` profiles. |
| `--pooled_faa` | – | Existing pooled protein FASTA. Skips the pooling step. |
| `--hmm_evalue` | `1e-5` | E-value reporting threshold for `hmmsearch`. |
| `--filter_evalue` | `1e-10` | Stricter E-value used to keep hits when filtering. |
| `--batch` | `batch01` | Name used for the gene catalogue outputs. |

Parameters can be given on the command line or collected in a YAML file:

```yaml
# params.yaml
input: samplesheet.csv
hmm_db: /data/hmm_profiles
filter_evalue: 1e-10
batch: batch01
```

```bash
nextflow run EdoT99/nf-gene-catalogue -r v1.0.0 -profile docker -params-file params.yaml
```

Command-line parameters override the params file, which overrides the config files.

## Usage

### Full run

```bash
nextflow run EdoT99/nf-gene-catalogue -r v1.0.0 -profile docker \
    --input samplesheet.csv \
    --hmm_db /path/to/hmm_profiles
```

### Starting from an existing pooled file

If you already have a pooled protein FASTA, skip the pooling step with `--pooled_faa`:

```bash
nextflow run EdoT99/nf-gene-catalogue -r v1.0.0 -profile docker \
    --pooled_faa /path/to/pooled_orf_predicted.faa \
    --hmm_db /path/to/hmm_profiles
```

### Resuming

Add `-resume` to reuse steps that already finished. For example, to try a different filtering threshold without rerunning HMMER:

```bash
nextflow run EdoT99/nf-gene-catalogue -r v1.0.0 -profile docker -params-file params.yaml \
    --filter_evalue 1e-20 -resume
```

## Profiles

Select one or more with `-profile`, separated by commas without spaces (e.g. `-profile test,docker`).

| Profile | Description |
|---|---|
| `docker` | Run every tool in a Docker container |
| `singularity` | Run every tool in a Singularity/Apptainer container (HPC) |
| `conda` | Build a conda environment for every tool |
| `test` | Small mock dataset to check the installation |

Machine-specific settings (CPUs, memory, Slurm queue, account) can go in your own config file, loaded with `-c`:

```bash
nextflow run EdoT99/nf-gene-catalogue -r v1.0.0 -profile singularity -c my_cluster.config -params-file params.yaml
```

See the Nextflow [configuration training](https://training.nextflow.io/latest/hello_nextflow/06_hello_config/) for how to write one.

## Output

Results are written to `results/` (change it with `-output-dir`):

```
results/
├── pooled_orfs/          pooled_orf_predicted.faa
├── hmmsearch/            <profile>.tblout / .domtblout / .out, one set per profile
├── hmm_filtered/
│   ├── hmm_hits.tsv          all ORF–profile hits passing --filter_evalue
│   ├── hmm_best_hits.tsv     best (lowest E-value) profile per ORF
│   └── hmm_hits.faa          sequences of the annotated ORFs
└── gene_catalog/<batch>/
    ├── <batch>_rep_seq.fasta    non-redundant representative sequences
    ├── <batch>_cluster.tsv      representative ⇥ member mapping
    ├── <batch>_all_seqs.fasta   all sequences grouped by cluster
    └── <batch>_kept_ids.txt     IDs of the representatives
```

## Citations

If you use this pipeline, please cite the tools it relies on:

- **HMMER**: Eddy SR. *Accelerated profile HMM searches.* PLoS Comput Biol (2011).
- **MMseqs2**: Steinegger M, Söding J. *MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets.* Nat Biotechnol (2017).
- **SeqKit**: Shen W, et al. *SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation.* PLoS ONE (2016).
- **Nextflow**: Di Tommaso P, et al. *Nextflow enables reproducible computational workflows.* Nat Biotechnol (2017).

## License

Released under the [MIT License](LICENSE).
