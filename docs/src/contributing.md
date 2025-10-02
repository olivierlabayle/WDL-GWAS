# Contributing and Reporting

## Contact

In order to report a bug, ask for a feature request or say hello, please fill an [issue](https://github.com/olivierlabayle/WDL-GWAS/issues).

## Contributing

Contributions are most welcome. This repository is organised as follows:

- The complete workflow is in `workflow.wdl`.
- All dependencies are built within a docker image defined by `docker/Dockerfile`.
- The "glue" code is written in [Julia](https://julialang.org/). Although contributions using other languages are acceptable, we prefer Julia as it keeps dependencies light and the repository clean.

## Running Tests

Running the tests locally requires a set of dependencies to be installed:

- [Julia](https://julialang.org/install/)
- [R](https://www.r-project.org/)
- [SuSiE](https://stephenslab.github.io/susieR/)
- [METAL](https://github.com/statgen/METAL)
- [plink2](https://www.cog-genomics.org/plink/2.0/)
- [cromwell](https://cromwell.readthedocs.io/en/stable/tutorials/FiveMinuteIntro/).

The tests further assume that the environment variable `CROMWELL_PATH` points to the cromwell jar file.

