module TestGWAS

using Test
using PopGen
using DataFrames
using CSV

PKGDIR = pkgdir(PopGen)
TESTDIR = joinpath(PKGDIR, "test")

@testset "Test merge-covariates-pcs" begin
    tmpdir = mktempdir()

    covariates = DataFrame(
        FID = string.(1:10),
        IID = string.(1:10),
        SUPERPOPULATION = ["EUR", "EUR", "AFR", "AFR", "AMR", "AMR", "EAS", "EAS", "SAS", "SAS"],
        COVID_19 = [1, missing, 0, 1, 0, 1, 1, missing, 0, 1],
        AGE = rand(20:80, 10),
    )
    CSV.write(joinpath(tmpdir, "covariates.csv"), covariates, delim="\t")

    for ancestry in ["AFR", "AMR", "EAS", "EUR", "SAS"]
        for chr in 1:3
            matching_covariates = covariates[covariates.SUPERPOPULATION .== ancestry, :]
            pcs = DataFrame(
                FID = matching_covariates.FID,
                IID = matching_covariates.IID,
                PC1 = randn(2),
                PC2 = randn(2),
            )
            rename!(pcs, "FID" => "#FID")
            CSV.write(joinpath(tmpdir, "pca.$ancestry.chr$(chr)_out.eigenvec"), pcs, delim="\t")
        end
    end

    covariates_file = joinpath(tmpdir, "covariates.csv")
    pcs_prefix = joinpath(tmpdir, "pca")
    copy!(ARGS, [
        "merge-covariates-pcs", 
        covariates_file,
        pcs_prefix,
        "--output", joinpath(tmpdir, "merged_covariates_and_pcs.tsv")
    ])
    julia_main()
    merged_covariates_pcs = CSV.read(joinpath(tmpdir, "merged_covariates_and_pcs.tsv"), DataFrame)
    # Merge did not add any row
    @test nrow(merged_covariates_pcs) == 10
    # 6 new columns 2 for PCS, 3 for get_chr_out_string
    @test names(merged_covariates_pcs) == [
        "FID",
        "IID",
        "SUPERPOPULATION",
        "COVID_19",
        "AGE",
        "CHR1_OUT_PC1",
        "CHR1_OUT_PC2",
        "CHR2_OUT_PC1",
        "CHR2_OUT_PC2",
        "CHR3_OUT_PC1",
        "CHR3_OUT_PC2"
    ]
    # missings are NA
    @test sum(merged_covariates_pcs.COVID_19 .== "NA") == 2
end

@testset "Test merge_chr_results" begin
    tmpdir = mktempdir()
    group = "AFR.SEVERE_COVID_19"
    chrs = 1:2
    gwas_merge_list = []
    for chr in chrs
        # Create dummy Regenie Step 2 results files
        gwas_results = DataFrame(
            CHROM = [chr],
            GENPOS = [1000],
            ID = ["chr$(chr):4132:G:A"],
            ALLELE0 = ["A"],
            ALLELE1 = ["T"],
            A1FREQ = [.5],
            N = [100],
            TEST = ["ADD"],
            BETA = [0.01],
            SE = [0.001],
            CHISQ = [10],
            LOG10P = [1.0],
            EXTRA = [""]
        )
        gwas_output_file = joinpath(tmpdir, "$group.chr$(chr).step2_SEVERE_COVID_19.regenie")
        CSV.write(gwas_output_file, gwas_results)
        push!(gwas_merge_list, gwas_output_file)
    end
    gwas_merge_list_file = joinpath(tmpdir, "gwas_merge_list.txt")
    open(gwas_merge_list_file, "w") do io
        for file in gwas_merge_list
            println(io, file)
        end
    end
    output_prefix = joinpath(tmpdir, "results.all_chr")
    copy!(ARGS, [
        "merge-chr-results",
        gwas_merge_list_file,
        "--output-prefix", output_prefix
    ])
    julia_main()
    gwas_results = CSV.read(output_prefix * ".tsv", DataFrame; delim="\t")
    expected_cols = [
        "CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", 
        "A1FREQ", "N", "TEST", "BETA", "SE", 
        "CHISQ", "LOG10P", "EXTRA"
    ]
    @test nrow(gwas_results) == 2
    @test Set(gwas_results.CHROM) == Set([1, 2])
    @test names(gwas_results) == expected_cols
end

end

true