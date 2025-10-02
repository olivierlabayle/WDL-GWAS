# WDL-GWAS

## Purpose

This [repository](https://github.com/olivierlabayle/WDL-GWAS) provides a WDL workflow that can be used to run reproducible genome-wide association studies and fine-mapping on the UKBiobank RAP and elsewhere.

!["Workflow"](assets/wdl-gwas-workflow.png)

## Installation

All dependencies require:

- Java, which can be installed with [sdkman](https://sdkman.io/install/).
- This repository which you can obtain [here](https://github.com/olivierlabayle/WDL-GWAS/releases).

### Installation for Local Usage

Running the workflow locally requires cromwell, [this page](https://cromwell.readthedocs.io/en/stable/tutorials/FiveMinuteIntro/) explains how to install it. Then, we recommend setting the environment variable `CROMWELL_PATH` to point to the newly downloaded jar file.

### Installation for UK Biobank Usage

In order to run the WDL-GWAS workflow on the UK Biobank RAP you need to setup [dxCompiler](https://github.com/dnanexus/dxCompiler#setup). This requires a few things:

- Installing the [dx toolkit](https://documentation.dnanexus.com/downloads).
- Downloading dxCompiler from the [releases page](https://github.com/dnanexus/dxCompiler/releases). Then, we recommend setting the environment variable `DX_COMPILER_PATH` to point to the newly downloaded jar file.

You can then login using your credentials via:

```bash
dx login
```