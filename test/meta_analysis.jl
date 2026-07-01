module TestMetaAnalysis

using Test
using PopGen
using DataFrames
using CSV

PKGDIR = pkgdir(PopGen)
TESTDIR = joinpath(PKGDIR, "test")

@testset "Test load_meta_analysis_worklist" begin
    tmpdir = mktempdir()
    gwas_results_list_file = joinpath(tmpdir, "regenie_files_list.txt")
    regenie_files_list = readdir(joinpath(TESTDIR, "assets", "meta_analysis"), join=true)
    open(gwas_results_list_file, "w") do io
        for f in regenie_files_list
            println(io, f)
        end
    end
    # No exclusion of groups
    exclude = []
    worklist = PopGen.load_meta_analysis_worklist(gwas_results_list_file; exclude = exclude)
    @test sort(worklist[!, [:PHENOTYPE, :GROUP]]) == DataFrame(
        PHENOTYPE = ["SEVERE_COVID_19", "SEVERE_PNEUMONIA", "SEVERE_PNEUMONIA", "SEVERE_PNEUMONIA"],
        GROUP = ["EAS", "AFR", "EUR", "SAS"]
    )
    # No exclusion of groups
    exclude =["EAS", "AFR"]
    worklist = PopGen.load_meta_analysis_worklist(gwas_results_list_file; exclude = exclude)
    @test sort(worklist[!, [:PHENOTYPE, :GROUP]]) == DataFrame(
        PHENOTYPE = ["SEVERE_PNEUMONIA", "SEVERE_PNEUMONIA"],
        GROUP = ["EUR", "SAS"]
    )
end

@testset "Test meta_analyse" begin
    # Two phenotypes are meta analysed
    ## - SEVERE_COVID_19 has 1 group
    ## - SEVERE_PNEUMONIA has 3 groups, 1 will be excluded
    tmpdir = mktempdir()
    output_prefix = joinpath(tmpdir, "gwas.meta_analysis")
    gwas_results_list_file = joinpath(tmpdir, "regenie_files_list.txt")
    regenie_files_list = readdir(joinpath(TESTDIR, "assets", "meta_analysis"), join=true)
    maf = 0.1
    open(gwas_results_list_file, "w") do io
        for f in regenie_files_list
            println(io, f)
        end
    end
    copy!(ARGS, 
        ["meta-analyse", 
        gwas_results_list_file,
        "--exclude=SAS",
        "--maf=$maf",
        "--output-prefix=$output_prefix"
    ])
    julia_main()
    expected_cols = Set([
        "ID", "BETA", "SE", "LOG10P", "DIRECTION", 
        "HET_ISQ", "HET_CHISQ", "HET_DF", "LOG10P_HET", 
        "CHROM", "POS", "ALLELE_0", "ALLELE_1", 
        "ALLELE_1_FREQ", "ALLELE_1_FREQ_STD", "ALLELE_1_FREQ_MIN",
        "ALLELE_1_FREQ_MAX",
        "N", "NGROUPS"
    ])
    # Check SEVERE_COVID_19
    meta_covid_19 = CSV.read(joinpath(tmpdir, "gwas.meta_analysis.SEVERE_COVID_19.gwas.tsv"), DataFrame)
    @test Set(names(meta_covid_19)) == expected_cols
    @test meta_covid_19[meta_covid_19.DIRECTION .== "0", :NGROUPS] == [0]
    @test all(meta_covid_19[meta_covid_19.DIRECTION .!= "0", :NGROUPS] .== 1)
    @test all(meta_covid_19.ALLELE_1_FREQ_MIN .> maf)
    @test all(meta_covid_19.ALLELE_1_FREQ_MAX .< 1 - maf)
    # Check SEVERE_PNEUMONIA
    meta_pneumonia = CSV.read(joinpath(tmpdir, "gwas.meta_analysis.SEVERE_PNEUMONIA.gwas.tsv"), DataFrame)
    @test Set(names(meta_pneumonia)) == expected_cols
    meta_pneumonia[meta_pneumonia.ID .== "chr1:14012312:T:C", :NGROUPS] == [1]
    @test all(meta_pneumonia.NGROUPS .<= 2)
    @test all(meta_pneumonia.ALLELE_1_FREQ_MIN .> maf)
    @test all(meta_pneumonia.ALLELE_1_FREQ_MAX .< 1 - maf)
    # Check plots have been created
    @test isfile(string(output_prefix, ".SEVERE_COVID_19.manhattan.png"))
    @test isfile(string(output_prefix, ".SEVERE_COVID_19.qq.png"))
    @test isfile(string(output_prefix, ".SEVERE_PNEUMONIA.manhattan.png"))
    @test isfile(string(output_prefix, ".SEVERE_PNEUMONIA.qq.png"))
end

end

true