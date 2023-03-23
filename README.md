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
