module TestGWASPLots

using Test
using PopGen
using DataFrames
using CSV

PKGDIR = pkgdir(PopGen)
TESTDIR = joinpath(PKGDIR, "test")

@testset "Test harmonize_gwas_results" begin
    # If the LOG10P contains NA it will be read as a string column
    ## The NaNs are filtered for plotting
    results = DataFrame(
        LOG10P = ["NA", "1", "1", "2"],
        CHROM = [1, 1, 2, 2],
        GENPOS = [1000, 2000, 3000, 4000],
        ID = ["rs1", "rs2", "rs3", "rs4"],
        A1FREQ = ["0.1", "0.2", "NA", "0.4"]
    )
    harmonized_results = PopGen.harmonize_gwas_results(results)
    @test harmonized_results.CHR == ["1", "1", "2", "2"]
    @test harmonized_results.BP == results.GENPOS
    @test harmonized_results.SNP == results.ID
    @test harmonized_results.P[1] === NaN
    @test harmonized_results.P[2:end] == [0.1, 0.1, 0.01]

    # If the LOG10P has no NA it will be read as a float column
    results = DataFrame(
        LOG10P = [1, 1, 2, 2],
        CHROM = [1, 1, 2, 2],
        GENPOS = [1000, 2000, 3000, 4000],
        ID = ["rs1", "rs2", "rs3", "rs4"],
        A1FREQ = [0.1, 0.2, 0.3, 0.4]
    )
    harmonized_results = PopGen.harmonize_gwas_results(results)
    @test harmonized_results.CHR == ["1", "1", "2", "2"]
    @test harmonized_results.BP == results.GENPOS
    @test harmonized_results.SNP == results.ID
    @test harmonized_results.P == [0.1, 0.1, 0.01, 0.01]
end

@testset "Test region_plot" begin
    finemapping_results = PopGen.harmonize_finemapping_results(
            CSV.read(
        joinpath(TESTDIR, "assets", "results", "results.all_chr.EUR.SEVERE_COVID_19.finemapping.tsv"), 
        DataFrame; 
        delim="\t"
    ))
    finemapping_results = finemapping_results[finemapping_results.LOCUS_ID .== "rs7515509", :]
    gwas_results = PopGen.harmonize_gwas_results(
            CSV.read(
        joinpath(TESTDIR, "assets", "results", "results.all_chr.EUR.SEVERE_COVID_19.gwas.tsv"), 
        DataFrame; 
        delim="\t"
    ))
    region_data = innerjoin(
        gwas_results,
        DataFrames.select(finemapping_results, [:ID, :REF, :ALT, :PIP, :CS, :LOCUS_ID, :PHASED_R2]), 
        on=[:ID]
    )
    fig = PopGen.region_plot(region_data)
    @test fig !== nothing
end

@testset "Test make_plots" begin
    tmpdir = mktempdir()
    output_prefix = joinpath(tmpdir, "plot")
    copy!(ARGS,[
        "make-plots", 
        joinpath(TESTDIR, "assets", "results", "results.all_chr.EUR.SEVERE_COVID_19.gwas.tsv"), 
        joinpath(TESTDIR, "assets", "results", "results.all_chr.EUR.SEVERE_COVID_19.finemapping.tsv"), 
        "--maf=0.01",
        "--output-prefix=$output_prefix"
        ]
    )
    julia_main()
    @test isfile(string(output_prefix, ".manhattan.png"))
    @test isfile(string(output_prefix, ".qq.png"))
    @test isfile(string(output_prefix, ".rs12732514.locuszoom.png"))
    @test isfile(string(output_prefix, ".rs7515509.locuszoom.png"))
end

@testset "Test region_plot" begin
    finemapping_results = PopGen.harmonize_finemapping_results(
            CSV.read(
        joinpath(TESTDIR, "assets", "results", "results.all_chr.EUR.SEVERE_COVID_19.finemapping.tsv"), 
        DataFrame; 
        delim="\t"
    ))
    finemapping_results = finemapping_results[finemapping_results.LOCUS_ID .== "rs7515509", :]
    gwas_results = PopGen.harmonize_gwas_results(
            CSV.read(
        joinpath(TESTDIR, "assets", "results", "results.all_chr.EUR.SEVERE_COVID_19.gwas.tsv"), 
        DataFrame; 
        delim="\t"
    ))
    region_data = innerjoin(
        gwas_results,
        DataFrames.select(finemapping_results, [:ID, :REF, :ALT, :PIP, :CS, :LOCUS_ID, :PHASED_R2]), 
        on=[:ID]
    )
    fig = PopGen.region_plot(region_data)
    @test fig !== nothing
end

end

true