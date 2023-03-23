# Genotype imputation using low depth sequencing data

## 1. Description
Pipeline for genotype imputation from low depth sequencing data using [GLIMPSE1](https://github.com/odelaneau/GLIMPSE). This pipeline uses [GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller) for estimating genotype probabilities at reference sites prior to imputation. The pipeline was developed using [Nextflow](https://www.nextflow.io/) and was tested on [SLURM](https://slurm.schedmd.com/documentation.html) job scheduler.

## 2. Prerequisites
The following software is required:
- Singularity
- Nextflow
- bcftools

## 3. Installation
To run this pipline you will need:
1. Download and install [GLIMPSE1](https://github.com/odelaneau/GLIMPSE)
2. Build GATK singularity container: `singularity build gatk_VERSION.sif docker://broadinstitute/gatk:VERSION`
3. Clone this repo: `git clone https://github.com/CERC-Genomic-Medicine/glimpse_pipeline.git`

## 4. Execution
1. Modify `nextflow.config` configuration file:
* `params.reference_vcfs` -- path to VCF/BCF files with phased reference panel genotypes. Each VCF/BCF file must have the corresponding tbi/csi index.
* `params.reference_sites_vcfs` -- path to sites-only VCF/BCF files of the reference panel. Each VCF/BCF file must have the corresponding tbi/csi index.
* `params.study_bams` -- path to BAM/CRAM files. One BAM/CRAM file per study participant. Each BAM/CRAM file must have the corresponding bai/crai index.
* `params.referenceDir` -- path to the folder with the reference genome *.fa file.
* `params.referenceGenome` -- name of the reference genome *.fa file (e.g. hs37d5.fa).
* `params.gatkContainer` -- path to the GATK singularity image file (.sif).
* `params.window_size` -- imputation window size in base-pairs. This is a parameter to the `GLIMPSE_chunk` executable. See [GLIMPSE1 documentation](https://odelaneau.github.io/GLIMPSE/glimpse1/tutorial_b38.html#run_chunk) for more details.
* `params.buffer_size` -- imputation window buffer size in base-pairs. This is a parameter to the `GLIMPSE_chunk` executable. See [GLIMPSE1 documentation](https://odelaneau.github.io/GLIMPSE/glimpse1/tutorial_b38.html#run_chunk) for more details.
* `params.chunk_exec` -- path to the `GLIMPSE_chunk` executable.
* `params.phase_exec` -- path to the `GLIMPSE_phase` executable.
* `params.ligate_exec` -- path to the `GLIMPSE_ligate` executable.
* `params.glimpse_maps` -- path to the GLIMPSE's genetic maps folder with the corresponding human genome build version (unarchived).
* `process.*` and `executor.*` -- set this arguments according to your compute cluster configuration.

2. Run pipleine. Example of interactive SLURM job:
```
salloc --time=12:00:00 --ntasks=1 --mem-per-cpu=16G
module load nextflow
module load singularity
module load bcftools
nextflow run Pipeline.nf
```
