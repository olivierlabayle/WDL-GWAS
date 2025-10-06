# Workflow's Outputs

The workflow generates a (potentially large) number of results files that can be a bit difficult to navigate at first.

The results of the main workflow are:

- The GWAS summary statistics, formatted as `GROUP_NAME.PHENOTYPE.gwas.tsv`. One for each group and phenotype.
- The GWAS plots, formatted as `GROUP_NAME.PHENOTYPE.PLOT_TYPE.png`. One for each group, phenotype and plot type (currently QQ, Manhattan).
- The finemapping summary statistics, formatted as `GROUP_NAME.PHENOTYPE.finemapping.tsv`. One for each group and phenotype.
- Optional locus plots, formatted as `GROUP_NAME.PHENOTYPE.LEAD_VARIANT_ID.locuszoom.png`. One for each group, phenotype and locus passing the fine-mapping thresholds.

If no groups are provided, the groupname is set to `GROUP_NAME=all`. Otherwise, and if the meta-analysis stage was requested, the following outputs will be created.

- Meta-Analysed GWAS results, formatted as `META_ANALYSIS.PHENOTYPE.gwas.tsv`. One for each phenotype.
- Meta-Analysed GWAS plots, formatted as `META_ANALYSIS.PHENOTYPE.PLOT_TYPE.png`. One for each phenotype and plot type (currently QQ, Manhattan).
- Meta-Analysis based finemapping results, formatted as `META_ANALYSIS.PHENOTYPE.finemapping.tsv`. One for each phenotype.
- Optional meta-analysis based locus plots, formatted as `META_ANALYSIS.PHENOTYPE.LEAD_VARIANT_ID.locuszoom.png`. One for each phenotype and locus passing the finemapping thresholds.

All results are found in the output folder whose lcoation and organisation depend on the execution engine.

## Local Execution with Cromwell

The outputs will be stored in the output folder defined by the `final_workflow_outputs_dir` variable specified in the cromwell's [option file](https://cromwell.readthedocs.io/en/latest/wf_options/Overview/). The output folder is typically organised in subfolders corresponding to the workflow's tasks.

For example, for a grouped analysis with meta-analysis you should have the following subfolders:

- `call-gwas_group_plots`: Containing plots for each analysed group.
- `call-gwas_meta_plots`: Containing the plots for the meta-analysis.
- `call-merge_gwas_group_chr_results`: Containing GWAS results for each group.
- `call-meta_analyse`: Containing meta-analysed GWAS results.
- `call-merge_fp_group_chr_results`: Containing finemapping results for each group.
- `call-merge_fp_meta_chr_results`: Containing finemapping results from the meta-analysis stage using summary statistics.

## UK Biobank RAP

In that case, the output folder is defined by the `destination` argument provided to `dx run`. Results will be found as a flat list in the output folder and all intermediate results are stored in the `intermediate` subfolder.