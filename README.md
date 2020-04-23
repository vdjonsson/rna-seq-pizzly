# Snakemake workflow: RNA-Seq Pizzly

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.12.3-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/rna-seq-pizzly.svg?branch=master)](https://travis-ci.org/snakemake-workflows/rna-seq-pizzly)

This workflow performs fusion detection and quantification by pseudoalignment using [kallisto](https://pachterlab.github.io/kallisto) and [pizzly](https://github.com/pmelsted/pizzly). 
The fusions detection algortihm was described in: 

Melsted, P. et al. [Fusion detection and quantification by pseudoalignment.](https://www.biorxiv.org/content/10.1101/166322v1) bioRxiv 166322 (2017) doi:10.1101/166322.

Workflow structure is modeled after [snakemake workflow template](https://github.com/snakemake-workflows/cookiecutter-snakemake-workflow) and [kallisto sleuth workflow](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth). 

## Outline 
1) Perform kallisto fusion on each sample separately. 
2) Perform pizzly on each sample separately. 
3) Combine pizzly fusion sequences across multiple samples. For repeated fusion sequences that come from different transcript pairs but the same gene pair, only first occurring transcript pair is kept, as downstream analysis focuses on gene-level changes and each unique transcript should be defined once. 
4) The filtered combined set of fusions are added to the original transcriptome to build a new kallisto index. 
5) Perform kallisto quant on each sample using new transcript-fusion index. 

## Authors

* Rachel Ng (@racng) - designed and authored workflow
* Vanessa Jonsson (@vdjonsson) - led analysis design & [Jonsson Lab](https://www.vdjonsson.com)

## Requirements:
- snakemake-≥5.12.3
- [kallisto](https://pachterlab.github.io/kallisto/download)
- [pizzly](https://github.com/pmelsted/pizzly)
- Python packages: json, pickle, pandas 
- R packages: rhdf5, tidyverse

## Usage

### Step 1: Install workflow

If you simply want to use this workflow, download and extract the [latest release](https://github.com/racng/rna-seq-pizzly/releases) or clone this repository:
```
git clone https://github.com/racng/rna-seq-pizzly.git
```
If you intend to modify and further extend this workflow or want to work under version control, fork this repository as outlined in [Advanced](#advanced). The latter way is recommended.

In any case, if you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository and, if available, its DOI (see above).

### Step 2: Configure workflow

Under the `rna-seq-pizzly/`, create a project folder. Copy the template configuration files `config.yaml` from `config/` to your project folder. 
```
mkdir project_name
cp config/config.yaml project_name/
```
Configure the workflow according to your needs via editing the configuration files.
Make sure to edit `workdir`, `refdir`, and `reference` and have the desired reference transcriptome fasta and gtf files downloaded into a folder in `refdir`.

Create `samples.tsv` in the project folder to include the names, FASTQ paths, and other metadata for samples that you want to include in the project. 
The sample sheet should have the following headers:
- sample: sample names
- fqdir: path to directory containing Read 1 and Read 2 FASTQ files
- fq1: filename of Read 1 FASTQ file
- fq2: filename of Read 2 FASTQ file
- Any additional metadata variables to be used in downstream differential expression analysis. (condition, gender, treatment, etc.)

### Step 3: Execute workflow

Navigate to the project folder containing `config.yaml` and `samples.tsv`. 
```
cd project_name
```
#### Using env modules:

Make sure env-modules names are correct. 
Test your configuration by performing a dry-run via

    snakemake -s ../workflow/Snakefile --use-envmodules -n 

Execute the workflow locally using `$N` cores via

    snakemake -s ../workflow/Snakefile --use-envmodules --cores $N

Execute on slurm cluster environment using `$n` maximum number of jobs:

    snakemake -s ../workflow/Snakefile --use-envmodules --profile ../profile/slurm -j $n --latency-wait 120 

#### Using conda:

Test your configuration by performing a dry-run via

    snakemake -s ../workflow/Snakefile --use-conda -n

Execute the workflow locally using `$N` cores via

    snakemake -s ../workflow/Snakefile --use-conda --cores $N

Execute on slurm cluster environment using `$n` maximum number of jobs:

    snakemake -s ../workflow/Snakefile --use-conda --profile ../profile/slurm -j $n --latency-wait 120 

See [Snakemake profiles](https://github.com/Snakemake-Profiles/slurm) for other cluster deployment. 

If you not only want to fix the software stack but also the underlying OS, `--use-singularity` in combination with any of the modes above.

See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

To run snakemake in the background and pipe the terminal output to file, wrap the above commands 

    nohup [INSERT SNAKEMAKE COMMAND] > out.text &

Check background job and its jobid with `jobs`. Bring back to foreground with `fg [jobid]`. 
### Step 4: Investigate results

After successful execution, you can create a self-contained interactive HTML report with all results via:

    snakemake -s ../workflow/Snakefile --report report.html

This report can, e.g., be forwarded to your collaborators.

### Advanced

The following recipe provides established best practices for running and extending this workflow in a reproducible way.

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the repo to a personal or lab account.
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the fork to the desired working directory for the concrete project/run on your machine.
3. [Create a new branch](https://git-scm.com/docs/gittutorial#_managing_branches) (the project-branch) within the clone and switch to it. The branch will contain any project-specific modifications (e.g. to configuration, but also to code).
4. Modify the config, and any necessary sheets (and probably the workflow) as needed.
5. Commit any changes and push the project-branch to your fork on github.
6. Run the analysis.
7. Optional: Merge back any valuable and generalizable changes to the [upstream repo](https://github.com/snakemake-workflows/rna-seq-pizzly) via a [**pull request**](https://help.github.com/en/articles/creating-a-pull-request). This would be **greatly appreciated**.
8. Optional: Push results (plots/tables) to the remote branch on your fork.
9. Optional: Create a self-contained workflow archive for publication along with the paper (snakemake --archive).
10. Optional: Delete the local clone/workdir to free space.

## Acknowledgment

Köster, Johannes and Rahmann, Sven. “Snakemake - A scalable bioinformatics workflow engine”. Bioinformatics 2012.

Creators of Snakemake software, workflow, templates, and profiles. 
- Johannes Köster (https://koesterlab.github.io) - snakemake, rna-seq-kallisto-sleuth
- David Lähnemann - rna-seq-kallisto-sleuth
- Per Unneberg (@percyfal) - slurm profile

## Testing

Under development
