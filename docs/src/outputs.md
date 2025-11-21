# Workflow's Outputs

The workflow generates a (potentially large) number of results files that can be a bit difficult to navigate at first.

In what follows, `GROUP_NAME` represents the value of a group if the `groupby` option was selected. If no groups are provided, the workflow sets it to `all`. Similarly `PHENOTYPE` is one of the provided `phenotypes` input.

## GWAS Summary Statistics

The GWAS summary statistics are stored in files with names `GROUP_NAME.PHENOTYPE.gwas.tsv`. One for each group and phenotype.

Each output file stores variant-level information and has the following columns:

- CHROM: The chromosome ID.
- POS: The variant position.
- ID: The variant ID.
- ALLELE_0: The baseline allele. This should be the REF allele provided in the input PGEN files.
- ALLELE_1: The effect allele, this should be the ALT allele provided in the input PGEN files.
- ALLELE_1_FREQ: The ALLELE_1 frequency.
- BETA: The effect size.
- SE: The standard error.
- LOG10P: The negative log10(p-value).
- N: The number of samples used for this variant.

Furthermore:

- **SAIGE** outputs the extra following: ALLELE_1_COUNT, MISSING_RATE, T_STAT, VAR.
- **REGENIE** outputs theextra following: TEST, CHISQ, EXTRA.

## Finemapping Summary Statistics

Similarly to GWAS results, finemapping summary statistics are stored in files with names `GROUP_NAME.PHENOTYPE.finemapping.tsv`. One for each group and phenotype.

The columns of those files are:

- CHROM: The chromosome ID.
- POS: The variant position.
- ID: The variant ID.
- REF: The reference allele from the input PGEN file
- ALT: The alternate allele from the input PGEN file
- LOCUS_ID: The GWAS lead variant ID identifying the locus.	
- PIP: The posterior inclusion probability.
- CS: The credible set ID.	
- UNPHASED_R2: The R2 with the GWAS lead variant ID.
- SUSIE_CONVERGED: Whether the SuSiE algorithm converged.

## Plots

Finally, the workflow outputs a series of plots:

- The GWAS plots, formatted as `GROUP_NAME.PHENOTYPE.PLOT_TYPE.png`. One for each group, phenotype and plot type (currently QQ, Manhattan).
- Optional locus plots, formatted as `GROUP_NAME.PHENOTYPE.LEAD_VARIANT_ID.locuszoom.png`. One for each group, phenotype and locus passing the fine-mapping thresholds.

## Meta Analysis

If the `groupby` and `meta_analysis` options are specified, the workflow will also produce meta-analysis outputs. These are similar to the above specified outputs with `GROUP_NAME=META_ANALYSIS`.

The GWAS output fill then contain the following additional columns:

- DIRECTION: Summary of effect direction for each study, with one '+' or '-' per study
- HET_ISQ: I^2 statistic which measures heterogeneity on scale of 0-100%
- HET_CHISQ: Chi-squared statistic in simple test of heterogeneity
- HET_DF: Degrees of freedom for heterogeneity statistic
- LOG10P_HET: log10 of p-value for heterogeneity statistic

## Results Location

All results are found in the output folder whose lcoation and organisation depend on the execution engine.

### Local Execution with Cromwell

The outputs will be stored in the output folder defined by the `final_workflow_outputs_dir` variable specified in the cromwell's [option file](https://cromwell.readthedocs.io/en/latest/wf_options/Overview/). The output folder is typically organised in subfolders corresponding to the workflow's tasks.

For example, for a grouped analysis with meta-analysis you should have the following subfolders:

- `call-gwas_group_plots`: Containing plots for each analysed group.
- `call-gwas_meta_plots`: Containing the plots for the meta-analysis.
- `call-merge_gwas_group_chr_results`: Containing GWAS results for each group.
- `call-meta_analyse`: Containing meta-analysed GWAS results.
- `call-merge_fp_group_chr_results`: Containing finemapping results for each group.
- `call-merge_fp_meta_chr_results`: Containing finemapping results from the meta-analysis stage using summary statistics.

### UK Biobank RAP

In that case, the output folder is defined by the `destination` argument provided to `dx run`. Results will be found as a flat list in the output folder and all intermediate results are stored in the `intermediate` subfolder.