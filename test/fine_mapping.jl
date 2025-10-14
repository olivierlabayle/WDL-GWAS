module TestFineMapping

using PopGen
using Test
using CSV
using DataFrames

PKGDIR = pkgdir(PopGen)
TESTDIR = joinpath(PKGDIR, "test")

@testset "Test locus_from_r2s" begin
    nvars = 20
    group = DataFrame(
        ID_B        = ["var$i" for i in 1:nvars],
        CHROM_B     = fill("1", nvars),
        POS_A       = fill(1325, nvars),
        UNPHASED_R2 = [0.09, 0.08, 0.02, 0.05, 0.15, 0.3,  0.4,  1,    0.6,  0.07, 0.5,  0.6,  0.01, 0.02, 0.05, 0.13, 0.09, 0.08, 0.04, 0.03],
        POS_B       = [900,  950,  980,  1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300]
    )
    # The lead is var8 at pos 1300, r2 is 0.1 and padding is 1
    ## window left boundary is 1150 - 1 variant: 1100
    ## window right boundary is 1900 + 1 variant: 2000
    @test PopGen.locus_from_r2s(group, "var8"; r2_threshold=0.1, padding=1) == ("var8", "1", 1100, 2000)
    # The lead is var8 at pos 1300, r2 is 0.55 and padding is 2
    ## no variant on left has sufficient r2, window left boundary is lead - 2 variant: 1200
    ## window right boundary is 1500 + 2 variant: 1700
    @test PopGen.locus_from_r2s(group, "var8"; r2_threshold=0.55, padding=2) == ("var8", "1", 1200, 1700)
    # The lead is var8 at pos 1300, r2 is 0.8 and padding is 2
    ## no variant on left has sufficient r2, window left boundary is lead - 2 variant: 1200
    ## no variant on right has sufficient r2,window right boundary is lead + 2 variant: 1400
    @test PopGen.locus_from_r2s(group, "var8"; r2_threshold=0.8, padding=2) == ("var8", "1", 1200, 1400)
end

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

@testset "Test get_window_idxs" begin
    variant_pos = [1, 2, 3, 4, 5, 6]
    locus = ("var3", "1", 3, 5)
    @test PopGen.get_window_idxs(locus, variant_pos) == (3, 5)
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
    sig_clumps = PopGen.write_significant_clumps(gwas_matched_pgen_prefix, gwas_results_file;
        min_sig_clump_size = 1,
        output = string(output_prefix, ".clumps.tsv"),
        lead_pvalue = 1e-8,
        p2_pvalue = 1e-5,
        r2_threshold = 0.3,
        clump_kb = 100,
        clump_id_field = "ID",
        clump_pval_field = "LOG10P",
        allele_1_field = "ALLELE1"
    )
    @test nrow(sig_clumps) == 0
    @test names(sig_clumps) == [
        "#CHROM", "POS", "ID", "NEG_LOG10_P", "TOTAL", "NONSIG", 
        "S0.05", "S0.01", "S0.001", "S0.0001", "SP2"
    ]
    @test CSV.read(string(output_prefix, ".clumps.tsv"), DataFrame; delim="\t") == sig_clumps
    # Write significant clumps: here we set very permissive parameters to get a clump
    sig_clumps = PopGen.write_significant_clumps(gwas_matched_pgen_prefix, gwas_results_file;
        min_sig_clump_size = 1,
        output = string(output_prefix, ".clumps.tsv"),
        lead_pvalue = 0.3,
        p2_pvalue = 0.3,
        r2_threshold = 0.,
        clump_kb = 300_000,
        clump_id_field = "ID",
        clump_pval_field = "LOG10P",
        allele_1_field = "ALLELE1"
    )
    @test nrow(sig_clumps) == 1
    @test length(split(first(sig_clumps.SP2), ",")) == 2
    @test CSV.read(string(output_prefix, ".clumps.tsv"), DataFrame; delim="\t") == sig_clumps
    # To get two clumps
    sig_clumps = PopGen.write_significant_clumps(gwas_matched_pgen_prefix, gwas_results_file;
        min_sig_clump_size = 1,
        output = string(output_prefix, ".clumps.tsv"),
        lead_pvalue = 0.3,
        p2_pvalue = 0.3,
        r2_threshold = 0.,
        clump_kb = 70_000,
        clump_id_field = "ID",
        clump_pval_field = "LOG10P",
        allele_1_field = "ALLELE1"
    )
end

@testset "Integration Test: finemap_locus" begin
    finemap_window_kb = 30_000
    pgen_prefix = joinpath(TESTDIR, "assets", "imputed", "chr1.qced")
    lead_ID = "chr1:40310265:G:A"
    r2_threshold = 0.1
    padding = 2
    # Test function components
    ## This is a degenerate case where the last variant is the lead variant
    r2s_with_lead = PopGen.get_r2_with_variant(lead_ID, pgen_prefix; ld_window_kb=finemap_window_kb)
    @test nrow(r2s_with_lead) == 4
    @test last(r2s_with_lead.ID_B) == lead_ID
    @test last(r2s_with_lead.UNPHASED_R2) == 1.0
    ## None of the variants have r2 > 0.1 with the lead
    locus = PopGen.locus_from_r2s(r2s_with_lead, lead_ID; r2_threshold=r2_threshold, padding=padding)
    @test locus == (lead_ID, 1, 18100537, 40310265)
    ## Extract dosages for the 3 variants in the locus
    X_df, variants_info = PopGen.genotypes_from_pgen(pgen_prefix, locus)
    @test Set(names(X_df)) == Set(["chr1:18100537:G:A", "chr1:22542609:T:C", "chr1:40310265:G:A", "FID", "IID"])
    @test variants_info.ID == ["chr1:18100537:G:A", "chr1:22542609:T:C", "chr1:40310265:G:A"]
    @test names(variants_info) == ["CHROM", "POS", "ID", "REF", "ALT"]
    for variant_id in variants_info.ID
        @test all(0 .<= X_df[!, variant_id] .<= 2)
    end
    ## Check merge is performed correctly
    y_df = X_df[1:1000, [:FID, :IID]]
    y_df[!, :Y] = rand(Bool, nrow(y_df))
    X, y = PopGen.get_susie_inputs(X_df, y_df)
    @test size(X) == (1000, 3)
    @test X isa Matrix{Float64}
    @test size(y) == (1000,)
    @test y isa Vector{Float64}
    ## Finemapping  works
    finemapping_results = PopGen.susie_finemap(X, y; n_causal=2)
    @test length(finemapping_results[:pip]) == 3
    # Post processing
    PopGen.postprocess_finemapping_results!(variants_info, finemapping_results, r2s_with_lead)
    @test names(variants_info) == PopGen.FINEMAPPING_RESULT_COLS
    @test nrow(variants_info) == 3
    @test any(ismissing.(variants_info.UNPHASED_R2)) == false
    # Test full function
    finemapping_results = PopGen.finemap_locus(lead_ID, pgen_prefix, y_df;
        n_causal=2,
        finemap_window_kb=finemap_window_kb,
        r2_threshold=r2_threshold,
        padding=padding
    )
    @test names(finemapping_results) == PopGen.FINEMAPPING_RESULT_COLS
    @test nrow(finemapping_results) == 3
end

@testset "Integration Test: finemap_locus_rss" begin
    finemap_window_kb = 30_000
    pgen_prefix = joinpath(TESTDIR, "assets", "imputed", "chr1.qced")
    lead_ID = "chr1:40310265:G:A"
    r2_threshold = 0.1
    padding = 2
    # Get LD matrix
    r2s_with_lead = PopGen.get_r2_with_variant(lead_ID, pgen_prefix; ld_window_kb=finemap_window_kb)
    locus = PopGen.locus_from_r2s(r2s_with_lead, lead_ID; r2_threshold=r2_threshold, padding=padding)
    R, variants = PopGen.get_LD_matrix(pgen_prefix, locus)
    @test R isa Matrix{Float64}
    @test size(R) == (3, 3)
    @test variants == ["chr1:18100537:G:A", "chr1:22542609:T:C", "chr1:40310265:G:A"]
    # Run susie_rss_finemap
    gwas_results = DataFrame(
        CHROM = ["1", "1", "1", "1", "1"],
        ID = ["chr1:14012312:T:C", "chr1:18100537:G:A", "chr1:22542609:T:C", "chr1:40310265:G:A", "toto"],
        GENPOS = [14012312, 18100537, 22542609, 40310265, 111622622],
        ALLELE0 = ["T", "G", "T", "G", "C"],
        ALLELE1 = ["C", "A", "C", "A", "A"],
        BETA = randn(5),
        SE = rand(5),
    )
    variants_info = PopGen.initialize_variants_info_rss(pgen_prefix, variants, gwas_results)
    @test variants_info.ID == ["chr1:18100537:G:A", "chr1:22542609:T:C", "chr1:40310265:G:A"]
    @test names(variants_info) == ["CHROM", "POS", "ID", "REF", "ALT", "BETA", "SE"]
    # Fake phenotype vector
    y = rand(100)
    finemapping_results = PopGen.susie_rss_finemap(R, variants_info, y; n_causal=1)
    @test length(finemapping_results[:pip]) == 3
    # Post processing
    PopGen.postprocess_finemapping_results!(variants_info, finemapping_results, r2s_with_lead)
    @test names(variants_info) == PopGen.FINEMAPPING_RESULT_COLS
    @test any(ismissing.(variants_info.UNPHASED_R2)) == false
    # Full function
    finemapping_results = PopGen.finemap_locus_rss(lead_ID, gwas_results, pgen_prefix, y;
        n_causal=1,
        finemap_window_kb=finemap_window_kb,
        r2_threshold=r2_threshold,
        padding=padding
    )
    @test names(finemapping_results) == PopGen.FINEMAPPING_RESULT_COLS
    @test nrow(finemapping_results) == 3
end


end

true