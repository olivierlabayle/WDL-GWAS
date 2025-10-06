# Running WDL-GWAS on the UK Biobank RAP

!!! note
    A more thorough example is coming soon!
    
## Describing the Inputs

As per the [Running WDL-GWAS Locally](@ref) example, the workflow's inputs are provided via a JSON file. However this time, the file paths need to point to teir location in your UK Biobank RAP project. There are two main ways you can fill those in:

- Via object identifiers: `dx://RAP_PROJECT_ID:FILE_ID`
- Via a file path: `dx://RAP_PROJECT_ID:/path/to/file`

An ecample inputs JSON file could look like:

```json
{
    "gwas.meta_exclude": ["AMR"],
    "gwas.covariates_file": "dx://RAP_PROJECT_ID:/path/to/covariates/file",
    "gwas.groupby": ["SUPERPOPULATION"],
    "gwas.covariates": ["AGE", "SEX", "AGE_x_AGE"],
    "gwas.phenotypes": ["SEVERE_PNEUMONIA", "SEVERE_COVID_19"],
    "gwas.genotypes": {
        "chr": "all",
        "bed": "dx://RAP_PROJECT_ID:/path/to/bed/file",
        "bim": "dx://RAP_PROJECT_ID:/path/to/bim/file",
        "fam": "dx://RAP_PROJECT_ID:/path/to/fam/file"
    },
    "gwas.imputed_genotypes": [
        {
            "chr": "1",
            "pgen": "dx://RAP_PROJECT_ID:PGEN_CHR_1_FILE_ID",
            "psam": "dx://RAP_PROJECT_ID:PSAM_CHR_1_FILE_ID",
            "pvar": "dx://RAP_PROJECT_ID:PVAR_CHR_1_FILE_ID"
        },
        {
            "chr": "2",
            "pgen": "dx://RAP_PROJECT_ID:PGEN_CHR_2_FILE_ID",
            "psam": "dx://RAP_PROJECT_ID:PSAM_CHR_2_FILE_ID",
            "pvar": "dx://RAP_PROJECT_ID:PVAR_CHR_2_FILE_ID"
        },
        {
            "chr": "3",
            "pgen": "dx://RAP_PROJECT_ID:PGEN_CHR_3_FILE_ID",
            "psam": "dx://RAP_PROJECT_ID:PSAM_CHR_3_FILE_ID",
            "pvar": "dx://RAP_PROJECT_ID:PVAR_CHR_3_FILE_ID"
        },
    ],
}
```

## Running the Workflow

First you need to compile the workflow and upload it to the RAP, this can be done with [dxCompiler](https://github.com/dnanexus/dxCompiler/blob/develop/doc/ExpertOptions.md):

```bash
java -jar $DX_COMPILER_PATH compile workflow.wdl \
-f -project $RAP_PROJECT_ID \
-reorg \
-folder /workflows/gwas \
-inputs inputs.json
```

where:

- `DX_COMPILER_PATH` is set as per the installation's instruction.
- `RAP_PROJECT_ID` is your UK Biobank RAP project ID.
- `inputs.json` is your workflow's input file. 

The compiler might output some warnings like `missing input for non-optional parameter` but you can ignore these.

The command will do a few things:

1. It will compile and upload the workflow to the UK Biobank RAP 
2. It will create a corresponding `inputs.dx.json` file locally which can be used to run the workflow 

Then, to run the workflow, run:

```bash
dx run -y \
-f inputs.dx.json \
--priority high \
--destination /workflow_outputs/ \
/workflows/gwas/gwas
```

Of course you can change the destination `/workflow_outputs` to your need. You should be able to visualize the run execution on your UK Biobank RAP account under the `MONITOR` tab.

## Outputs

Outputs will be stored in the `/workflow_outputs//` folder. Since the RAP preserves intermediate files, this will be quite dense.

- plots: searching for `png`
- group/phenotype results: searching for `regenie.results.${group_name}.${phenotype}.tsv`