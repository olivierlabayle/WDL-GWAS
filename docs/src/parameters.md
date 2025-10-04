# Workflow Parameters

## Input Arguments

- `covariates_file`: A CSV file containing the set of covariates and phenotypes. Missing values can be either NAN, NA, NULL or an empty character.
- `genotypes`: A PLINK BED fileset in GRCH38.
- `imputed_genotypes`: A set of PLINK PGEN filesets, one for each chromosome in GRCH38.

## General Options

- `phenotypes`: The set of phenotypes for which each GWAS will be run.
- `groupby` (default: []): A set of variables used to stratify individuals for which a GWAS will be run independently. If empty, the full dataset is used.
- `filterby` (default: []):
- `covariates` (default: ["AGE", "SEX", "AGE_x_AGE", "AGE_x_SEX"]): A set of covariates used to adjust for confounding or increase power in the association testing step. Product of variables can be defined using the `_x_` syntax, for example: ["AGE", "SEX", "AGE_x_SEX", "AGE_x_AGE"].
- `min_cases_controls` (default: 10): Minimum number of phenotype cases/controls within a group to proceed to GWAS.

## PCA Options

- `ip_values` (default: "1000 50 0.05"): Values used to create independent genotypes for PCA (see [here](https://www.cog-genomics.org/plink/2.0/ld)).
- `npcs` (default 10): Number of principal components to use to account for population structure.
- `approx_pca` (default: true): Whether to use an approximation to the PCA algorithm (see [here](https://www.cog-genomics.org/plink/2.0/strat)).
- `maf` (default: 0.01): Minor allele frequency threshold. This is used to filter variants for PCA and for plotting GWAS results. The association testing step is still performed across all variants.

## GWAS Options

- `mac` (default: 10): Minor allele count used to filter variants (also used for PCA).
- `regenie_cv_folds` (default: 5): Number of folds for Regenie step 1. Can also be `loocv`.
- `regenie_bsize` (default 1000): Regenie block size.

For Regenie's options see the [online documentation](https://rgcgithub.github.io/regenie/options/).

## Finemapping Options

## Miscellaneous Options

- `docker_image` (default: olivierlabayle/wdl-gwas:main): The docker image used to run the workflow.

