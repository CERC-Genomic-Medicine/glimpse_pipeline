# Genotype imputation using low depth sequencing data

## 1. Description
Pipeline for genotype imputation from low-depth sequencing data using [GLIMPSE1](https://github.com/odelaneau/GLIMPSE). This pipeline uses [GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller) for estimating genotype probabilities at reference sites before imputation. The pipeline was developed using [Nextflow](https://www.nextflow.io/) and was tested on [SLURM](https://slurm.schedmd.com/documentation.html) job scheduler.

> [!IMPORTANT]
> As recommended by the GLIMPSE1 authors, this pipeline imputes indels while ignoring indels' genotype likelihoods computed from the sequencing data. If you are interested in more accurate indel imputation, then consider switching to GLIMPSE2 (but keep in mind that GLIMPSE2 doesn't support a joint imputation).

## 2. Prerequisites
The following software is required:
- apptainer
- Nextflow
- bcftools

## 3. Installation
To run this pipline you will need:
1. Download and install [GLIMPSE1](https://github.com/odelaneau/GLIMPSE):
   ```
   git clone https://github.com/odelaneau/GLIMPSE.git
   cd GLIMPSE/
   git checkout glimpse1
   cd ..
   ```
   Then, follow the documentation at [GLIMPSE webpage](https://odelaneau.github.io/GLIMPSE/glimpse1/installation.html).

2. Build GATK singularity container (see available versions at [broadinstitute/gatk](https://hub.docker.com/r/broadinstitute/gatk): `apprainer build gatk_VERSION.sif docker://broadinstitute/gatk:VERSION`
3. Clone this repo: `git clone https://github.com/CERC-Genomic-Medicine/glimpse_pipeline.git`

## 4. Execution
1. Modify `nextflow.config` configuration file:
* `params.reference_vcfs` -- path to VCF/BCF files with phased reference panel genotypes. Each VCF/BCF file must have the corresponding tbi/csi index. The genotypes must be split by chromosome and contain prefixes: `chr1`, `chr2`, ..., `chr22`, `chrX_nonpar`, `chrX_par1`, `chrX_par2`.
* `params.reference_sites_vcfs` -- path to sites-only VCF files of the reference panel. Each VCF file must have the corresponding tbi/csi index. The sites must be split by chromosome in the same way as the genotypes.
* `params.study_bams` -- path to BAM/CRAM files. One BAM/CRAM file per study participant. Each BAM/CRAM file must have the corresponding bai/crai index.
* `params.study_sample_ploidy` -- path to the text file listing sample ploidy in chromosome X. No header; space-delimited; one sample per line; two columns - sample name and number of chrX copies (1 or 2). You can estimate the ploidy of chromosome X from chromosome X sequencing depth. You can ignore this argument if you are not interested in imputing chromosome X.
* `params.referenceDir` -- path to the folder with the reference genome *.fa file.
* `params.referenceGenome` -- name of the reference genome *.fa file (e.g. hs37d5.fa).
* `params.gatkContainer` -- path to the GATK singularity image file (.sif).
* `params.window_size` -- imputation window size in base-pairs. This is a parameter to the `GLIMPSE_chunk` executable. See [GLIMPSE1 documentation](https://odelaneau.github.io/GLIMPSE/glimpse1/tutorial_b38.html#run_chunk) for more details.
* `params.buffer_size` -- imputation window buffer size in base-pairs. This is a parameter to the `GLIMPSE_chunk` executable. See [GLIMPSE1 documentation](https://odelaneau.github.io/GLIMPSE/glimpse1/tutorial_b38.html#run_chunk) for more details.
* `params.chunk_exec` -- path to the `GLIMPSE_chunk` executable.
* `params.phase_exec` -- path to the `GLIMPSE_phase` executable.
* `params.ligate_exec` -- path to the `GLIMPSE_ligate` executable.
* `params.glimpse_maps` -- path to the GLIMPSE's genetic maps folder with the corresponding human genome build version.
* `process.*` and `executor.*` -- set this arguments according to your compute cluster configuration.

2. Run pipleine. Example of interactive SLURM job:
```
salloc --time=12:00:00 --ntasks=1 --mem-per-cpu=16G
module load nextflow
module load apptainer
module load bcftools
nextflow run Imputation.nf
```
