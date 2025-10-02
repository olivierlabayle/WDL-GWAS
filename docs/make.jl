using Documenter

makedocs(;
    authors="Olivier Labayle <olabayle@gmail.com> and contributors",
    sitename="WDL-GWAS",
    format=Documenter.HTML(;
        canonical="https://olivierlabayle.github.io/WDL-GWAS/stable/",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Tutorials" => [
            "WDL-GWAS Locally" => "tutorials/local_run.md",
            "WDL-GWAS on UKBiobank RAP" => "tutorials/ukb_run.md",
        ],
        "Workflow Parameters" => "parameters.md",
        "Troobleshooting" => "troubleshooting.md",
        "Contributing and Reporting" => "contributing.md",
    ],
)

deploydocs(;
    repo="github.com/olivierlabayle/WDL-GWAS",
    devbranch="main",
    push_preview=true
)
