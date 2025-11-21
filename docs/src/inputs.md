# Workflow's Inputs

!!! note
    Due to some cromwell parsing issues, all options and arguments must be provided as strings unless stated otherwise. That is, in your input file you should write `npcs="10"`, not `npcs=10`.

## Required Arguments

These arguments must be provided or the workflow will not run.

- `covariates_file`: A CSV file containing the set of covariates and phenotypes. Missing values can be either NAN, NA, NULL or an empty character.
- `genotypes`: A PLINK BED fileset in GRCh38 containing variants for all chromosomes. This is only used for PCA, so you can either used the typed variants or an appropriate subset of the imputed genotypes.
- `imputed_genotypes`: A set of PLINK2 PGEN filesets, one for each chromosome in GRCh38.
- `phenotypes`: A list of phenotypes for which each GWAS will be run independently.

!!! warning
    For `genotypes` and `imputed_genotypes` it is required that each variant `ID` is set and unique.

## General Options

These options describe the general behaviour of the workflow.

- `groupby` (default: []): A list of variables used to stratify individuals. If empty, the full dataset is used.
- `filterby` (default: []): A list of conditions on covariates present in the dataset to filter individuals. At the moment this can either be in the format:
  - `COLUMN=value`: which can be used to filter for individuals with categorical traits, e.g., `SUPERPOPULATION=AFR` if the `covariates_file` has a `SUPERPOPULATION` column.
  - `COLUMN>value`: which can be used to filter individuals with continuous variables, e.g., `AGE>70`. Any of `<`, `>`, `=`, `<=`, `>=` is a valid operator.
  - A combination of the above, e.g., `filterby=[SUPERPOPULATION=AFR, AGE>70]`.
- `covariates` (default: `[AGE, SEX, AGE_x_AGE, AGE_x_SEX]`): A list of covariates used to adjust for confounding or increase power in the association testing step. Product of variables can be defined using the `_x_` syntax, for example: `AGE_x_SEX`.
- `min_cases_controls` (default: 10): For binary traits only, the minimum number of cases/controls within a group to proceed to the GWAS.

## PCA Options

PCA is performed using [plink2's PCA](https://www.cog-genomics.org/plink/2.0/strat).

- `npcs` (default 10): Number of principal components to use to account for population structure.
- `ip_values` (default: `1000 50 0.05`): A string of values used to create independent genotypes for PCA (see [here](https://www.cog-genomics.org/plink/2.0/ld)).
- `approx_pca` (default: true): Whether to use an approximation to the PCA algorithm (see [here](https://www.cog-genomics.org/plink/2.0/strat)). Turning this to `false` if used only for small datasets, for instance during testing.
- `loco_pca` (default: false): Whether principal components should be computed in a LOCO fashion to be used as part of the GWAS covariates (cannot be used with `gwas_software=saige`).


## GWAS Options

We offer two different options to run GWAS, [SAIGE](https://saigegit.github.io/SAIGE-doc/) (default) and [REGENIE](https://rgcgithub.github.io/regenie/options/). Please refer to their respective documentation for the below options.

- `gwas_software` (default: `saige`): One of `regenie` or `saige`.
- `regenie_cv_folds` (default: `loocv`): Number of folds for Regenie step 1. Any integer is valid.
- `regenie_bsize` (default 1000): Regenie block size.

## Fine-Mapping Options

Finemapping proceeds in two stages. First clumps are formed using [plink2 LD-based result clumping](https://www.cog-genomics.org/plink/2.0/postproc) and lead variants are identified. Then a window is formed around the lead variant to be further finemapped with [SuSiE](https://stephenslab.github.io/susieR/).

- `min_sig_clump_size` (default: `10`): Defines the minimum number of variants within a clump for a locus to be considered for finemapping.
`lead_pvalue` (default: `5e-8`): A clump's lead variant must have at least this p-value.
- `p2_pvalue` (default: `5e-5`): Other variants in the clump must have at least this p-value.
- `r2_threshold` (default: `0.1`): The LD between the lead variant and other variants in the clump must also as strong as `r2_threshold`.
- `clump_kb` (default: `250`): The maximum distance between the lead variant and a clumped variant.
- `n_causal` (default: `10`): SuSiE's `L` variable, the expected number of causal variants in a locus.
- `finemap_strategy` (default: `rss`): Either `rss` fine-mapping will use the summary stats from the GWAS step. Any other option (e.g., `all`) leads to individual level data-based fine-mapping.
- `susie_max_iter` (default `1000`): Maximum number of iterations used by the SuSiE algorithm.

## Meta-Analysis Options

Meta analysis is performed using [METAL](https://github.com/statgen/METAL).

- `meta_analysis` (default: `length(groupby) > 0`): The only parameter that must be provided as a boolean (not a string). It defaults to true if groups are provided and false otherwise.
- `meta_exclude` (default: ["ADMIXED"]): Whether some groups should be excluded from the meta-analysis stage. Note that this will result in lower power. We default to exclusing admixed individuals if the grouping column contains superpopulations.
- `meta_method` (default: "STDERR"). The METAL's meta analysis method, please refer to their documentation for more information.

## Miscellaneous Options

- `maf` (default: `0.01`): Minor allele frequency threshold. This is used to filter variants for PCA and for plotting GWAS results. The association testing step is always performed across all variants. Post fitering for downstream analysis is left to the user.
- `mac` (default: `10`): Minor allele count used to filter variants in REGENIE and PCA QC.

## Developpement Options

In principle you only need to change these if you are trying to contribute to the workflow.

- `docker_image` (default: `olivierlabayle/wdl-gwas:main`): The docker image used to run the workflow. You shouldn't need to change that if you are not contributing to the workflow.
- `julia_use_sysimage` (default: `true`): Whether to use the Julia's system image built within the container with [PackageCompiler](https://julialang.github.io/PackageCompiler.jl/stable/). During development it is better to have this to false to reflect the current state of the code.
- `julia_threads` (default: `auto`): Julia's number of [threads](https://docs.julialang.org/en/v1/manual/multi-threading/). Some external tools seem to crash the parent process when multiple threading is used locally. Typically set to false in the development stage.

