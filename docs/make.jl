using PopGen
using Documenter

DocMeta.setdocmeta!(PopGen, :DocTestSetup, :(using PopGen); recursive=true)

makedocs(;
    modules=[PopGen],
    authors="Olivier Labayle <olabayle@gmail.com> and contributors",
    sitename="WDL-GWAS",
    format=Documenter.HTML(;
        canonical="https://olivierlabayle.github.io/WDL-GWAS/stable/",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/olivierlabayle/WDL-GWAS",
    devbranch="main",
)
