module TestHarmonizeGWASResults

using Test
using PopGen
using CSV
using DataFrames

PKGDIR = pkgdir(PopGen)
TESTDIR = joinpath(PKGDIR, "test")

@testset "Test harmonize_gwas_results" begin
    shared_colnames = ["CHROM", "POS", "ID", "ALLELE_0", "ALLELE_1", "ALLELE_1_FREQ", "BETA", "SE", "LOG10P", "N"]
    regenie_only_colnames = ["TEST", "CHISQ", "EXTRA"]
    saige_only_colnames = ["ALLELE_1_COUNT", "MISSING_RATE", "T_STAT", "VAR"]
    tmpdir = mktempdir()
    # Test REGENIE
    gwas_results_file = joinpath(TESTDIR, "assets", "gwas_software_outputs", "all.AGE.chr1_AGE.regenie")
    output_file = joinpath(tmpdir, "regenie.harmonized.tsv")
    copy!(ARGS, [
        "harmonize-gwas-results",
        gwas_results_file,
        "--source-software=regenie",
        "--output=$output_file"
    ])
    julia_main()
    regenie_source = CSV.read(gwas_results_file, DataFrame)
    regenie_harmonized = CSV.read(output_file, DataFrame)
    ## Nothing is lost
    @test size(regenie_source) == size(regenie_harmonized)
    ## Expected colnames
    @test names(regenie_harmonized) == vcat(shared_colnames, regenie_only_colnames)

    # Test SAIGE
    gwas_results_file = joinpath(TESTDIR, "assets", "gwas_software_outputs", "all.AGE.chr1.saige.txt")
    output_file = joinpath(tmpdir, "saige.harmonized.tsv")
    copy!(ARGS, [
        "harmonize-gwas-results",
        gwas_results_file,
        "--source-software=saige",
        "--output=$output_file"
    ])
    julia_main()
    saige_source = CSV.read(gwas_results_file, DataFrame)
    saige_harmonized = CSV.read(output_file, DataFrame)
    ## Nothing is lost
    @test size(saige_source) == size(saige_harmonized)
    ## Expected colnames
    @test names(saige_harmonized) == vcat(shared_colnames, saige_only_colnames)
end

end

true