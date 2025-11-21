module TestFineMapping

using PopGen
using Test
using CSV
using DataFrames

PKGDIR = pkgdir(PopGen)
TESTDIR = joinpath(PKGDIR, "test")

@testset "Test make_clean_sample_file" begin
    tmpdir = mktempdir()
    # The file contains actual sample ids
    sample_file = joinpath(tmpdir, "sample_file.txt")
    open(sample_file, "w") do io
        write(io, "odap100\todap100")
    end
    @test PopGen.make_clean_sample_file(sample_file) == sample_file
    # The file contains a list of files containing sample ids written by REGENIE
    sample_file_1 = joinpath(tmpdir, "gwas.individual_ids.EUR.PHENO_1.txt")
    write(sample_file_1, "odap1\todap1")
    sample_file_2 = joinpath(tmpdir, "gwas.individual_ids.AFR.PHENO_1.txt")
    write(sample_file_2, "odap2\todap2")
    sample_file_3 = joinpath(tmpdir, "gwas.individual_ids.AMR.PHENO_1.txt")
    write(sample_file_3, "odap3\todap3")
    sample_file_4 = joinpath(tmpdir, "gwas.individual_ids.EUR.PHENO_2.txt")
    write(sample_file_4, "odap4\todap4")

    sample_file = joinpath(tmpdir, "samples_file.txt")
    open(sample_file, "w") do io
        for filename in [sample_file_1, sample_file_2, sample_file_3, sample_file_4]
            println(io, filename)
        end
    end

    output_sample_file = PopGen.make_clean_sample_file(sample_file; exclude=["AMR"], phenotype="PHENO_1")
    @test readlines(output_sample_file) == ["odap1\todap1", "odap2\todap2"]
end

@testset "Test loci_from_clumps" begin
    gwas_matched_pgen_prefix = joinpath(TESTDIR, "assets", "imputed", "chr1.qced")
    # These two clumps overlap and will be merged
    sig_clumps = DataFrame(
        Dict(
            "#CHROM" => [1, 1],
            "POS" => [14012312, 92682820],
            "NEG_LOG10_P" => [1, 2],
            "ID" => ["chr1:14012312:T:C", "chr1:92682820:C:T"],
            "SP2" => ["chr1:9694126:C:T,chr1:18100537:G:A,chr1:40310265:G:A", "chr1:22542609:T:C"]
        )
    )
    loci = PopGen.loci_from_clumps(sig_clumps, gwas_matched_pgen_prefix; padding=1)
    @test only(loci) == (lead_id = "chr1:92682820:C:T", neglog10pval = 2, locus_start = 9694126, locus_end = 111622622, chr = 1)
    # Only one clump
    sig_clumps = DataFrame(
        Dict(
            "#CHROM" => [1],
            "POS" => [183905563],
            "NEG_LOG10_P" => [1],
            "ID" => ["chr1:183905563:G:A"],
            "SP2" => ["chr1:231799576:C:T,chr1:111622622:C:A"]
        )
    )
    loci = PopGen.loci_from_clumps(sig_clumps, gwas_matched_pgen_prefix; padding=2)
    @test only(loci) == (lead_id = "chr1:183905563:G:A", neglog10pval = 1, locus_start = 40310265, locus_end = 231799576, chr = 1)
end

@testset "Integration Test: finemap_significant_regions" begin
    gwas_results_file = joinpath(TESTDIR, "assets", "results", "regenie.results.group.phenotype.tsv")
    pgen_prefix = joinpath(TESTDIR, "assets", "imputed", "chr1.qced")
    sample_file = joinpath(TESTDIR, "assets", "results", "sample_file.txt")
    covariates_file = joinpath(TESTDIR, "assets", "covariates", "covariates.csv")
    tmpdir = mktempdir()

    gwas_results = CSV.read(gwas_results_file, DataFrame)
    pvar = CSV.read(pgen_prefix * ".pvar", DataFrame; delim='\t', comment="##")
    # Identifying the variants that were not tested in GWAS
    PopGen.tag_variant_id_missing_from_gwas!(pvar, gwas_results)
    gwas_variants = filter(!endswith(".not_in_gwas"), pvar.ID)
    @test gwas_variants == [
        "chr1:14012312:T:C",
        "chr1:18100537:G:A",
        "chr1:22542609:T:C",
        "chr1:40310265:G:A",
        "chr1:92682820:C:T",
        "chr1:111622622:C:A"
    ]
    # Write those variants to a temporary PGEN fileset
    gwas_matched_pgen_prefix = joinpath(tmpdir, "chr1.gwas_matched")
    PopGen.write_new_pgen_from_gwas_results(pgen_prefix, gwas_matched_pgen_prefix, pvar, sample_file)
    @test countlines(gwas_matched_pgen_prefix * ".psam") == 1130
    @test countlines(gwas_matched_pgen_prefix * ".pvar") == 7
    # Write significant clumps: here none will be found but the function should still return an empty DataFrame
    output_prefix = joinpath(tmpdir, "fine_mapping_test")
    sig_clumps = PopGen.get_significant_clumps(gwas_matched_pgen_prefix, gwas_results_file;
        min_sig_clump_size = 1,
        output = string(output_prefix, ".clumps.tsv"),
        lead_pvalue = 1e-8,
        p2_pvalue = 1e-5,
        r2_threshold = 0.3,
        clump_kb = 100,
        clump_id_field = "ID",
        clump_pval_field = "LOG10P",
        allele_1_field = "ALLELE_1"
    )
    @test nrow(sig_clumps) == 0
    @test names(sig_clumps) == [
        "#CHROM", "POS", "ID", "NEG_LOG10_P", "TOTAL", "NONSIG", 
        "S0.05", "S0.01", "S0.001", "S0.0001", "SP2"
    ]
    @test CSV.read(string(output_prefix, ".clumps.tsv"), DataFrame; delim="\t") == sig_clumps
    # Write significant clumps: here we set very permissive parameters to get a clump
    sig_clumps = PopGen.get_significant_clumps(gwas_matched_pgen_prefix, gwas_results_file;
        min_sig_clump_size = 1,
        output = string(output_prefix, ".clumps.tsv"),
        lead_pvalue = 0.3,
        p2_pvalue = 0.3,
        r2_threshold = 0.,
        clump_kb = 300_000,
        clump_id_field = "ID",
        clump_pval_field = "LOG10P",
        allele_1_field = "ALLELE_1"
    )
    @test nrow(sig_clumps) == 1
    sig_clump = first(sig_clumps)
    @test length(split(sig_clump.SP2, ",")) == 2
    loci = PopGen.loci_from_clumps(sig_clumps, gwas_matched_pgen_prefix)
    locus = only(loci)
    @test locus.locus_start <= sig_clump.POS <= locus.locus_end
    @test locus.lead_id == sig_clump.ID
    @test locus.chr == sig_clump["#CHROM"]
    @test locus.neglog10pval == sig_clump.NEG_LOG10_P
    # To get two clumps
    sig_clumps = PopGen.get_significant_clumps(gwas_matched_pgen_prefix, gwas_results_file;
        min_sig_clump_size = 1,
        output = string(output_prefix, ".clumps.tsv"),
        lead_pvalue = 0.3,
        p2_pvalue = 0.3,
        r2_threshold = 0.,
        clump_kb = 70_000,
        clump_id_field = "ID",
        clump_pval_field = "LOG10P",
        allele_1_field = "ALLELE_1"
    )
    @test nrow(sig_clumps) == 2
end

@testset "Integration Test: finemap_locus" begin
    pgen_prefix = joinpath(TESTDIR, "assets", "imputed", "chr1.qced")
    lead_ID = "chr1:40310265:G:A"
    padding = 2
    locus = (locus_start=18100537, locus_end=183905563, chr=1, lead_id=lead_ID)

    # Test function components
    lead_to_locus_r2 = PopGen.compute_lead_to_locus_r2(locus, pgen_prefix)
    @test nrow(lead_to_locus_r2) == 6
    @test lead_to_locus_r2[lead_to_locus_r2.ID_B .== lead_ID, :UNPHASED_R2] == [1]
    ## Extract dosages for the 3 variants in the locus
    X_df, variants_info = PopGen.genotypes_from_pgen(pgen_prefix, locus)
    @test Set(names(X_df)) == Set([lead_to_locus_r2.ID_B..., "FID", "IID"])
    @test variants_info.ID == lead_to_locus_r2.ID_B
    @test names(variants_info) == ["CHROM", "POS", "ID", "REF", "ALT"]
    for variant_id in variants_info.ID
        @test all(0 .<= X_df[!, variant_id] .<= 2)
    end
    ## Check merge is performed correctly
    y_df = X_df[1:1000, [:FID, :IID]]
    y_df[!, :Y] = rand(Bool, nrow(y_df))
    X, y = PopGen.get_susie_inputs(X_df, y_df)
    @test size(X) == (1000, 6)
    @test X isa Matrix{Float64}
    @test size(y) == (1000,)
    @test y isa Vector{Float64}
    ## Finemapping  works
    finemapping_results = PopGen.susie_finemap(X, y; n_causal=2)
    @test length(finemapping_results[:pip]) == 6
    # Post processing
    PopGen.postprocess_finemapping_results!(variants_info, finemapping_results, lead_to_locus_r2)
    @test names(variants_info) == PopGen.FINEMAPPING_RESULT_COLS
    @test nrow(variants_info) == 6
    @test any(ismissing.(variants_info.UNPHASED_R2)) == false
    # Test full function
    finemapping_results = PopGen.finemap_locus(locus, pgen_prefix, y_df;n_causal=2)
    @test names(finemapping_results) == PopGen.FINEMAPPING_RESULT_COLS
    @test nrow(finemapping_results) == 6
end

@testset "Integration Test: finemap_locus_rss" begin
    pgen_prefix = joinpath(TESTDIR, "assets", "imputed", "chr1.qced")
    lead_ID = "chr1:40310265:G:A"
    padding = 2
    locus = (locus_start=18100537, locus_end=183905563, chr=1, lead_id=lead_ID)
    # Get LD matrix
    lead_to_locus_r2 = PopGen.compute_lead_to_locus_r2(locus, pgen_prefix)
    R, variants = PopGen.get_LD_matrix(pgen_prefix, locus)
    @test R isa Matrix{Float64}
    @test size(R) == (6, 6)
    @test variants == lead_to_locus_r2.ID_B
    # Run susie_rss_finemap
    gwas_results = DataFrame(
        CHROM = ["1", "1", "1", "1", "1", "1"],
        ID = ["chr1:183905563:G:A", "chr1:18100537:G:A", "chr1:22542609:T:C", "chr1:40310265:G:A", "chr1:92682820:C:T", "chr1:111622622:C:A"],
        POS = [183905563, 18100537, 22542609, 40310265, 92682820, 111622622],
        ALLELE_0 = ["G", "G", "T", "G", "C", "C"],
        ALLELE_1 = ["A", "A", "C", "A", "T", "A"],
        BETA = randn(6),
        SE = rand(6),
    )
    variants_info = PopGen.initialize_variants_info_rss(pgen_prefix, variants, gwas_results)
    @test variants_info.ID == variants
    @test names(variants_info) == ["CHROM", "POS", "ID", "REF", "ALT", "BETA", "SE"]
    # Fake phenotype vector
    y = rand(100)
    finemapping_results = PopGen.susie_rss_finemap(R, variants_info, y; n_causal=1)
    @test length(finemapping_results[:pip]) == 6
    # Post processing
    PopGen.postprocess_finemapping_results!(variants_info, finemapping_results, lead_to_locus_r2)
    @test names(variants_info) == PopGen.FINEMAPPING_RESULT_COLS
    @test any(ismissing.(variants_info.UNPHASED_R2)) == false
    # Full function
    finemapping_results = PopGen.finemap_locus_rss(locus, gwas_results, pgen_prefix, y;n_causal=1)
    @test names(finemapping_results) == PopGen.FINEMAPPING_RESULT_COLS
    @test nrow(finemapping_results) == 6
end


end

true