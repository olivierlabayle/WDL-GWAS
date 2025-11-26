module TestGWASPLots

using Test
using PopGen
using DataFrames
using CSV

PKGDIR = pkgdir(PopGen)
TESTDIR = joinpath(PKGDIR, "test")

@testset "Test genetics_makie_gwas_harmonize" begin
    results = DataFrame(
        LOG10P = [missing, 1, 1, 2],
        CHROM = [1, 1, 2, 2],
        POS = [1000, 2000, 3000, 4000],
        ID = ["rs1", "rs2", "rs3", "rs4"],
        ALLELE_1_FREQ = [0.1, 0.2, missing, 0.4]
    )
    harmonized_results = PopGen.genetics_makie_gwas_harmonize(results)
    @test harmonized_results.CHR == ["1", "1", "2", "2"]
    @test harmonized_results.BP == results.POS
    @test harmonized_results.SNP == results.ID
    @test harmonized_results.P[1] === missing
    @test harmonized_results.P[2:end] == [0.1, 0.1, 0.01]
end

@testset "Test region_plot" begin
    finemapping_results = PopGen.genetics_makie_fp_harmonize(
            CSV.read(
        joinpath(TESTDIR, "assets", "results", "results.all_chr.EUR.SEVERE_COVID_19.finemapping.tsv"), 
        DataFrame; 
        delim="\t"
    ))
    finemapping_results = finemapping_results[finemapping_results.LOCUS_ID .== "rs7515509", :]
    gwas_results = PopGen.genetics_makie_gwas_harmonize(
            CSV.read(
        joinpath(TESTDIR, "assets", "results", "results.all_chr.EUR.SEVERE_COVID_19.gwas.tsv"), 
        DataFrame; 
        delim="\t"
    ))
    region_data = innerjoin(
        gwas_results,
        DataFrames.select(finemapping_results, [:ID, :REF, :ALT, :PIP, :CS, :LOCUS_ID, :UNPHASED_R2, :SUSIE_CONVERGED]), 
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

end

true