# czid-workflows - portable CZ ID production pipeline logic

#### Infectious Disease Sequencing Platform
CZ ID is a hypothesis-free global software platform that helps scientists identify pathogens in metagenomic sequencing
data.

- **Discover** - Identify the pathogen landscape
- **Detect** - Monitor and review potential outbreaks
- **Decipher** - Find potential infecting organisms in large datasets

CZ ID is a collaborative open project of [Chan Zuckerberg Initiative](https://www.chanzuckerberg.com/) and
[Chan Zuckerberg Biohub](https://czbiohub.org).

## Workflows
Currently we have 5 main workflows. The details of each pipeline are in a README in each of the workflow folders. 

* [short-read-mngs](workflows/short-read-mngs/README.md) 
* [consensus-genome](workflows/consensus-genome/README.md)
* [phylotree-ng](workflows/phylotree-ng/README.md)
* [long-read-mngs](workflows/long-read-mngs/README.md)
* [amr](workflows/amr/README.md)

## Running these workflows
This repository contains [WDL](https://openwdl.org/) workflows that the [CZ ID](https://czid.org) platform uses in
production. See [Running WDL workflows locally](https://github.com/chanzuckerberg/czid-workflows/wiki/Running-WDL-workflows-locally)
to get started with them.

### Running short and long-read-mngs on Skyline HPC

Execution of the **short-read-mngs** (`Illumina`) and **long-read-mngs** (`ONT`) CZID workflows on Skyline HPC is orchestrated through [cluster/run-czid](cluster/run-czid) which invokes Slurm submissions for all samples defined in a samplesheet. 

#### HPC Configuration
* Uses miniwdl with Slurm and Singularity configuration:  [cluster/miniwdl.cfg](cluster/miniwdl.cfg)
* Container images are available on this fork's GHCR as czid-short-read-mngs and czid-long-read-mngs: [images](https://github.com/orgs/OpenOmics/packages?repo_name=czid-workflows)

---
#### Prepare a samplesheet

Prepare a CSV samplesheet with your input data that looks as follows. Both short and long-read-mngs can be run in the same samplesheet:

[example_samplesheet.csv](cluster/example_samplesheet.csv):

```csv
sample_name,seq_input,seq_platform,seq_fastq_R1,seq_fastq_R2,host_genome
sample_PE_1,DNA,Illumina,/path/to/fastqs/sample1_R1.fastq.gz,/path/to/fastqs/sample1_R2.fastq.gz,human
sample_SE_2,DNA,Illumina,/path/to/fastqs/sample2_R1.fastq.gz,,human
sample_3,DNA,ONT,/path/to/fastqs/sample3.fastq.gz,,human
sample_4,RNA,ONT,/path/to/fastqs/sample4.fastq.gz,,mosquito
```

Required columns must contain a value for every sample. Optional columns may be omitted or left blank.

| Column | Required? | Description |
| --- | --- | --- |
| `sample_name` | Yes | Unique sample name |
| `seq_input` | Yes | `DNA` or `RNA` (case-insensitive) |
| `seq_platform` | Yes | `Illumina` for short-read-mngs or `ONT` for long-read-mngs (case-insensitive) |
| `seq_fastq_R1` | Yes | Absolute FASTQ file path, or a path relative to `--fastq-dir` |
| `seq_fastq_R2` | For paired-end Illumina | R2 FASTQ file path, using the same path rules as R1. Leave blank for single-end reads or ONT |
| `host_genome` | No | Defaults to `human` when blank or omitted. Other options being added are `mosquito`, `bat`, and `macaca_fascicularis` |


Alternatively, if all fastq files are in a shared directory, `seq_fastq_R1` and `seq_fastq_R2` can use filenames relative to the shared FASTQ directory which is passed with `--fastq-dir`. See the 2nd example in next section:

`samplesheet.csv`

```csv
sample_name,seq_input,seq_platform,seq_fastq_R1,seq_fastq_R2,host_genome
sample_1,DNA,Illumina,sample1_R1.fastq.gz,sample1_R2.fastq.gz,human
sample_2,RNA,ONT,sample2.fastq.gz,,mosquito
```
---
#### Submit jobs

[cluster/run-czid](cluster/run-czid) takes at least 2 inputs:
1. CSV Samplesheet file
2. Output directory passed with `--output-dir`

For HPC use, enclose the execution command within a sbatch Slurm script.

```bash
/data/openomics/prod/czid-workflows/cluster/run-czid \
    example_samplesheet.csv \
    --output-dir /path/to/output/dir
```

If the samplesheet contains filenames relative to the shared FASTQ directory, supply `--fastq-dir /absolute/path/to/fastqs`:

```bash
/data/openomics/prod/czid-workflows/cluster/run-czid \
    samplesheet.csv \
    --fastq-dir /absolute/path/to/fastqs \
    --output-dir /path/to/output/dir
```

Each sample's directory under `--output-dir` contains `submit.sbatch` with the miniwdl workflow invocation, Slurm `.out`/`.err` logs, and run directories with workflow results.

---
#### Main pipeline reference inputs

- The orchestrator selects the workflow and host-filtering references.

- Reference files are sourced from CZID's public S3 buckets. The ncbi-indexes-prod version is 2024-02-06.
  - [short-read config YAML](cluster/skyline_short-read-mngs_config_2024-02-06.yml)
  - [long-read config YAML](cluster/skyline_long-read-mngs_config_2024-02-06.yml)

- **Data provenance for taxonomy annotation**
  - Both workflows use the **2026-11-01** NCBI taxdump alongwith TaxonKit for annotating the results report with name, category, 7-rank lineage. This matches the taxids with the 2026-11-01 taxonomy database from: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/.
  - The reference sequence databases used for analysis are the 2024-02-06 reference files from CZID's buckets.

#### Main pipeline outputs

| Output | Short-read mNGS | Long-read mNGS |
| --- | --- | --- |
| Taxon report | `short_read_report.csv`: NT/NR counts, reads per million (RPM), alignment metrics, and contig metrics | `long_read_report.csv`: NT/NR counts, bases per million (BPM), alignment metrics, and contig metrics |

For data provenance with 2026-11-01 taxonomy annotation, the workflows additionally output metadata containing taxid resolution differences if present.
| Output | Purpose |
| --- | --- |
| `<report>.csv.taxonomy.csv` | Per-row audit of original and resolved taxids, resolution status, and rank/genus changes. |
| `<report>.csv.taxonomy.json` | Annotation method, TaxonKit version, and reference taxonomy date |


### System Requirements 
The system requirements differs for each workflow and depending on the database being used. For example running the short-read-mngs workflow with the full NT and NR databases would require an instance with >1TB of disk space and >100GB of memory. Running other workflows (e.g. consensus-genome, amr) requires much less space. 

If using an ARM Mac try following [these](RunningCZIDWorkflowsOnARMMacs.md) setup instructions
### Software requirements
* docker with buildx support (version >= 19.03)
* python3 
* virtualenv
* requirements-dev.txt - to automatically install this run `make python-dependencies`

### Quick Setup
To get setup, first set the workflow you want to run with 

```export WORKFLOW=<workflow-name>``` e.g.

```export WORKFLOW=amr```

You can see available workflows with `make ls`

Either `build` or `pull` the workflow docker image with 

```make pull  ## The faster option``` or 

```make build ## The slower option, but necessary if you're modifying the docker container```

### Running a workflow
Run a workflow with 

```make run```

Which simply runs the ```miniwdl run path_to_wdl.wdl``` command with some defaults

Each workflow has a number of required and optional inputs, and all require at least an input file (usually a fastq). Default inputs are set from the `workflows/<workflow-name>/test/local_test.yml` file. These may or may not be accurate for every analysis. You can override these defaults and add your own with:

```make run INPUT='-i your_file_here.yml'```

If you're happy with the defaults, you can add arguments to the `miniwdl` command using 

```make run EXTRA_INPUTS='input_fastq=/path/to/input.fastq' ```

### An example
Lets say I want to run a consensus-genome workflow. I would run the following:

```
export WORKFLOW=consensus-genome
make pull # pull the latest docker container from github packages
make python-dependencies # create a .venv and install the requirements-dev.txt dependencies 

make run EXTRA_INPUTS='fastqs_0=workflows/consensus-genome/test/sample_sars-cov-2_paired_r1.fastq.gz \
        fastqs_1=workflows/consensus-genome/test/sample_sars-cov-2_paired_r2.fastq.gz \
        technology="Illumina" \
        sample="my_sample_name"'
```
## CI/CD

We use GitHub Actions for CI/CD. Lint and unit tests run on GitHub from jobs in `.github/workflows/wdl-ci.yml`
(triggered on every commit).

## Contributing

This project adheres to the Contributor Covenant code of conduct. By participating, you are expected to uphold this code. Please report unacceptable behavior to opensource@chanzuckerberg.com.

## Security

Please disclose security issues responsibly by contacting security@chanzuckerberg.com.
