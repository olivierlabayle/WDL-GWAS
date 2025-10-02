# Running WDL-GWAS on the UK Biobank RAP

First you need to compile the WDL workflow and upload it to the RAP, this can be done with the following:

```bash
java -jar $DX_COMPILER_PATH compile rap_workflows/gwas/workflow.wdl \
-f -project $PROJECT_ID \
-reorg \
-folder /workflows/gwas \
-inputs rap_workflows/gwas/inputs.json
```

where the `DX_COMPILER_PATH` and `PROJECT_ID` have to be set appropriately. The compiler might output some warnings like `missing input for non-optional parameter` but you can ignore these.

!!! warning "Updating workflows"
    At this point in time it seems like compiling multiple times the same workflow does not replace the old files. You will need to manually erase them from the RAP.

Then, you can run the workflow with the following command

```bash
dx run -y \
-f rap_workflows/gwas/inputs.dx.json \
--priority high \
--destination /gwas_outputs_covid_meta/ \
/workflows/gwas/gwas
```

## Outputs

Outputs will be stored in the `/gwas_outputs/` folder. There will be many files as I haven't figured out how to organise them yet. You can filter them to access, as a few examples:

- plots: searching for `png`
- group/phenotype results: searching for `regenie.results.${group_name}.${phenotype}.tsv`