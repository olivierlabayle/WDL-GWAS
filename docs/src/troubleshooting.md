# Troubleshooting

## Understanding Local Workflow failures

Cromwell's logs are quite dense and difficult to navigate, when the workflow errors, look for a filename ending in `execution/stderr`. This file will contain more information on the error.

## Understanding UK Biobank RAP failures

### Navigating the Logs.

First make sure you have read [this guide](https://dnanexus.gitbook.io/uk-biobank-rap/working-on-the-research-analysis-platform/managing-jobs/troubleshooting-guide) and know how to search for a problem within the produced logs.

### Running Interactively on UKBiobank RAP

To debug errors and provide constructive error reports, it may be useful to run some code interactively. For this, you can use the [Cloud Workstation](https://documentation.dnanexus.com/developer/cloud-workstation). This [tutorial](https://academy.dnanexus.com/interactivecloudcomputing/cloudworkstation) may also be useful. To start an instance:

```bash
dx run \
  --instance-type mem2_ssd1_v2_x8 \
  -imax_session_length="10h" \
  -y \
  --ssh app-cloud_workstation
```

To import one or multiple files (typically input files for a failed step):

```bash
dx download FILE_ID_1 FILE_ID_2...
```

Then you can download the docker image and enter a container:

```bash
docker run -it --rm -v $PWD:/mnt/data olivierlabayle/wdl-gwas:main /bin/bash
```

The current directory is mounted to `/mnt/data`. From there, work as usual, for instance to start a Julia REPL:

```bash
julia --project=/opt/PopGen --sysimage=/opt/PopGen/sysimage.so --startup-file=no
```

You can upload generated data to your project with `dx upload`:

```bash
dx upload --path "$DX_PROJECT_CONTEXT_ID:" <FILE>
```

Finally, when you are finished, terminate the job with the appropriate job-id:

```bash
dx terminate job-J1V4870JpYQP94jgb33y45qP
```