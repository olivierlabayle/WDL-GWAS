# Workflow Parameters


- `docker_image` (default: olivierlabayle/wdl-gwas:main): The docker image used to run the workflow.
- `covariates_file`: A CSV file containing the set of covariates and phenotypes. Missing values should be represented by empty fields.
- `genotypes`: PLINK BED files.
- `imputed_genotypes`: PLINK PGEN files, split by chromosome.
- `groupby` (default: []): A set of variables used to stratify individuals for which a GWAS will be run independently. If empty, the full dataset is used.
- `covariates` (default: ["AGE", "SEX", "AGE_x_AGE", "AGE_x_SEX"]): A set of covariates used to adjust for confounding or increase power in the association testing step. Product of variables can be defined using the `_x_` syntax, for example: ["AGE", "SEX", "AGE_x_SEX", "AGE_x_AGE"].
- `phenotypes`: A set of phenotypes for which a GWAS will be run independently.
- `min_cases_controls` (default: 10): Minimum number of phenotype cases/controls within a group to proceed to GWAS.
- `high_ld_regions`: File containing high LD regions to be excluded when performing LD pruning for PCA. The file is stored in `assets/exclude_b38.txt` and needs to be uploaded to the RAP.
- `ip_values` (default: "1000 50 0.05"): Values used to create independent genotypes for PCA (see [here](https://www.cog-genomics.org/plink/2.0/ld)).
- `npcs` (default 10): Number of principal components to use to account for population structure.
- `approx_pca` (default: true): Whether to use an approximation to the PCA algorithm (see [here](https://www.cog-genomics.org/plink/2.0/strat)).
- `maf` (default: 0.01): Minor allele frequency threshold. This is used to filter variants for PCA and for plotting GWAS results. The association testing step is still performed across all variants.
- `mac` (default: 10): Minor allele count used to filter variants for PCA and entering REGENIE's variants.
- `regenie_cv_folds` (default: 5): Number of folds for Regenie step 1. Can also be `loocv`.
- `regenie_bsize` (default 1000): Regenie block size.

For Regenie's options see the [online documentation](https://rgcgithub.github.io/regenie/options/).
