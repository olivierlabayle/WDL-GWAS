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
    saige_binary_only_colnames = ["N_CASES", "N_CONTROLS", "PVAL_NA", "IS_SPA", "AF_CASES", "AF_CONTROLS", "N_CASES_HOM", "N_CASES_HET", "N_CONTROLS_HOM", "N_CONTROLS_HET"]

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

    # Test SAIGE continuous trait
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

    # Test SAIGE binary trait
    gwas_results_file = joinpath(tmpdir, "all.BINARY_TRAIT.chr1.saige.txt")
    saige_output = Dict(
        "CHR" => [1], 
        "POS" => [2], 
        "MarkerID" => ["rs134"], 
        "Allele1" => ["A"], 
        "Allele2" => ["C"], 
        "AC_Allele2" => [1], 
        "AF_Allele2" => [1], 
        "MissingRate" => [0.1], 
        "BETA" => [1], 
        "SE" => [1], 
        "Tstat" => [1], 
        "var" => [1], 
        "p.value" => [0.1], 
        "p.value.NA" => [0.1], 
        "Is.SPA" => [false], 
        "AF_case" => [0.1], 
        "AF_ctrl" => [0.1], 
        "N_case" => [1], 
        "N_ctrl" => [2], 
        "N_case_hom" => [1], 
        "N_case_het" => [1], 
        "N_ctrl_hom" => [1], 
        "N_ctrl_het" => [1]
    )
    CSV.write(gwas_results_file, saige_output, missingstring="NA")
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
    @test size(saige_harmonized) == (1, length(saige_output) + 1) # N column added
    ## Expected colnames
    @test names(saige_harmonized) == vcat(shared_colnames, saige_only_colnames, saige_binary_only_colnames)
end

end

true