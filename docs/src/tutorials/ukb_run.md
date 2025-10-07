# Running WDL-GWAS on the UK Biobank RAP

!!! note
    This page is under development!

In this example, we provide a tutorial to reproduce the two GWAS studies performed in the WDL-GWAS paper: BMI and colorectal cancer. This involves a bit more than running the workflow since we need to prepare the data for it.

## Create the Covariates File

First we need to extract the relevant fields from the UK Biobank's main dataset, we can do this using `dx extract_dataset`:

```bash
dx extract_dataset DATASET_RECORD_ID --fields-file docs/src/assets/paper_fields.txt --output raw_covariates.csv
```

However, the extracted dataset cannot be directly used with the WDL-GWAS workflow, we first need to define the covariates and phenotypes from the extracted fields. The following Julia code snippet creates such a covariate file. Here is what id does:

- It makes sure the `FID` and `IID` are defined and correspond to the `eid` field.
- It defines the current participants' age.
- It defines colorectal cancer has having any of the C18, C19 or C20 ICD10 code.
- It extracts the BMI column

```julia
using CSV
using DataFrames
using Dates

function parse_ICD10_col!(parsed_col, col, codes)
    for index in eachindex(parsed_col)
        val = col[index]
        val === missing && continue
        parsed_col[index] =  parsed_col[index] | any(code -> startswith(val, code), codes)
    end
    return parsed_col
end

function parse_cancer(cols...; codes)
    cancer_col = zeros(Int, length(first(cols)))
    for col in cols
        parse_ICD10_col!(cancer_col, col, codes)
    end
    return cancer_col
end

cancer_registry_cols = [
    "participant.p40006_i0",
    "participant.p40006_i1",
    "participant.p40006_i2",
    "participant.p40006_i3",
    "participant.p40006_i4",
    "participant.p40006_i5",
    "participant.p40006_i6",
    "participant.p40006_i7",
    "participant.p40006_i8",
    "participant.p40006_i9",
    "participant.p40006_i10",
    "participant.p40006_i11",
    "participant.p40006_i12",
    "participant.p40006_i13",
    "participant.p40006_i14",
    "participant.p40006_i15",
    "participant.p40006_i16",
    "participant.p40006_i17",
    "participant.p40006_i18",
    "participant.p40006_i19",
    "participant.p40006_i20",
    "participant.p40006_i21"
]

raw_covariates = CSV.read("raw_covariates.csv", DataFrame)

covariates = select(raw_covariates,
    "participant.eid" => "FID",
    "participant.eid" => "IID",
    "participant.p22001" => "SEX",
    "participant.p34" => (x -> Int(Dates.year(now())) .- x) => "AGE",
    "participant.p23104_i0" => "BMI",
    cancer_registry_cols => ((cols...) -> parse_cancer(cols...; codes=["C18", "C19", "C20"])) => "COLORECTAL_CANCER"
)

CSV.write("covariates.csv", covariates)
```

Finally we can upload the dataset to the UK Biobank RAP using `dx upload`:

```bash
dx upload --path /wdl_gwas_covariates.csv covariates.csv
```

## Extracting Imputed and Typed Genotypes

WDL-GWAS requires both imputed and typed genotypes. For that we will use the TOPMed imputed genotypes which should be provided in your `/Bulk/Imputation/` folder.

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