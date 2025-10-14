const FINEMAPPING_RESULT_COLS = [
    "CHROM", "POS", "ID", "REF", "ALT", "LOCUS_ID", "PIP", "CS", "UNPHASED_R2", "SUSIE_CONVERGED"
]

function tag_variant_id_missing_from_gwas!(pvar, gwas_results)
    # Map GWAS variants to their alleles
    gwas_ids_to_alleles = Dict(
        row.ID => Set([row.ALLELE0, row.ALLELE1]) 
        for row in Tables.namedtupleiterator(DataFrames.select(gwas_results, [:ID, :ALLELE0, :ALLELE1]))
    )
    # Update variants IDs in PVAR if they are not in the GWAS variants
    pvar.ID = map(Tables.namedtupleiterator(pvar)) do row
        if haskey(gwas_ids_to_alleles, row.ID)
            gwas_alleles = gwas_ids_to_alleles[row.ID]
            pvar_alleles = Set([row.ALT, row.REF])
            if gwas_alleles != pvar_alleles
                string(row.ID, ".unmatched_alleles")
            else
                row.ID
            end
        else
            string(row.ID, ".not_in_gwas")
        end
    end
end

function write_new_pgen_from_gwas_results(input_pgen_prefix, output_pgen_prefix, pvar, sample_file)
    # Write PVAr with tagged variants
    tmpdir = mktempdir()
    tagged_pvar_file = joinpath(tmpdir, "updated_variants.pvar")
    CSV.write(tagged_pvar_file, pvar; delim='\t')
    # Write variants to exclude to file
    exclude_list = filter(x -> endswith(x, ".not_in_gwas") || endswith(x, ".unmatched_alleles"), pvar.ID)
    exclude_file = joinpath(tmpdir, "exclude_variants.txt")
    open(exclude_file, "w") do io
        for id in exclude_list
            println(io, id)
        end
    end
    # Run PLINK filtering cmd
    run(`plink2 \
        --pgen $input_pgen_prefix.pgen \
        --psam $input_pgen_prefix.psam \
        --pvar $tagged_pvar_file \
        --make-pgen --out $output_pgen_prefix \
        --exclude $exclude_file \
        --keep $sample_file
    `)
end

function write_significant_clumps(pgen_prefix, gwas_results_file;
    min_sig_clump_size = 3,
    output = "clumps.sig.tsv",
    lead_pvalue = 5e-8,
    p2_pvalue = 5e-5,
    r2_threshold = 0.5,
    clump_kb = 250,
    clump_id_field = "ID",
    clump_pval_field = "LOG10P",
    allele_1_field = "ALLELE_1"
    )
    tmpdir = mktempdir()
    output_clump_prefix = joinpath(tmpdir, "clumps")
    run(`plink2 --pfile $pgen_prefix \
        --clump $gwas_results_file \
        --clump-p1 $lead_pvalue \
        --clump-p2 $p2_pvalue \
        --clump-r2 $r2_threshold \
        --clump-kb $clump_kb \
        --clump-id-field $clump_id_field \
        --clump-log10 \
        --clump-p-field $clump_pval_field \
        --clump-a1-field $allele_1_field \
        --out $output_clump_prefix
    `)
    clumps_file = output_clump_prefix * ".clumps"
    clumps = isfile(clumps_file) ? 
        CSV.read(clumps_file, DataFrame; delim="\t") : 
        DataFrame([col => [] for col in ["#CHROM", "POS", "ID", "NEG_LOG10_P", "TOTAL", "NONSIG", "S0.05", "S0.01", "S0.001", "S0.0001", "SP2"]])
    sig_clumps = filter(
        :SP2 => x -> x !== "." && length(split(x, ",")) >= min_sig_clump_size, 
        clumps
    )
    CSV.write(output, sig_clumps, delim="\t")
    return sig_clumps
end

function locus_from_r2s(group, group_key; r2_threshold=0.1, padding=10)
    lead_idx = findfirst(==(group_key), group.ID_B)
    chr = first(group.CHROM_B)
    # Define locus start
    locus_start_idx = findfirst(>=(r2_threshold), group.UNPHASED_R2)
    if locus_start_idx === nothing || locus_start_idx > lead_idx
        locus_start_idx = lead_idx
    end
    locus_start_idx = max(1, locus_start_idx - padding)
    locus_start_pos = group.POS_B[locus_start_idx]

    # Define locus end
    locus_end_idx = findlast(>=(r2_threshold), group.UNPHASED_R2)
    if locus_end_idx === nothing || locus_end_idx < lead_idx
        locus_end_idx = lead_idx
    end
    locus_end_idx = min(nrow(group), locus_end_idx + padding)
    locus_end_pos = group.POS_B[locus_end_idx]

    return (group_key, chr, locus_start_pos, locus_end_pos)
end

function genotypes_from_pgen(pgen_prefix, locus)
    locus_id, chr, locus_start, locus_end = locus
    tmpdir = mktempdir()
    output_prefix = joinpath(tmpdir, "dosage")
    run(`plink2 \
        --pfile $pgen_prefix \
        --from-bp $locus_start \
        --to-bp $locus_end \
        --chr $chr \
        --export A \
        --out $output_prefix`
    )
    X_df = CSV.read(string(output_prefix, ".raw"), DataFrame;
        delim='\t',
        drop = ["PAT", "MAT", "SEX", "PHENOTYPE"],
        missingstring="NA"
    )
    # Imputes missing values and get column names without the trailing _REF
    new_colnames = map(names(X_df)) do colname
        if colname ∈ ("FID", "IID")
            colname
        else
            μ = mean(skipmissing(X_df[!, colname]))
            X_df[!, colname] = coalesce.(X_df[!, colname], μ)
            join(split(colname, "_")[1:end-1], "_")
        end
    end
    pvar = PopGen.read_pvar(pgen_prefix * ".pvar")
    pvar = pvar[locus_start .<= pvar.POS .<= locus_end, :]
    @assert new_colnames[3:end] == pvar.ID
    rename!(X_df, new_colnames)
    rename!(pvar, "#CHROM" => "CHROM")
    return X_df, pvar
end

function get_window_idxs(locus, variants_pos)
    _, _, locus_start, locus_end = locus
    idx_inf = findfirst(==(locus_start), variants_pos)
    idx_sup = findlast(==(locus_end), variants_pos)
    return idx_inf, idx_sup
end

function load_phenotypes_matching_samples(covariates_file, sample_file, phenotype)
    sample_list = read_samples(sample_file)
    covariates = CSV.read(covariates_file, DataFrame; 
        delim='\t', 
        select=["FID", "IID", phenotype]
    )
    return innerjoin(covariates, sample_list, on=[:FID, :IID])
end

function get_susie_inputs(X_df, y_df)
    data = innerjoin(X_df, y_df, on=[:FID, :IID])
    variant_cols = filter(x -> x ∉ ("FID", "IID"), names(X_df))
    phenotype_col = only(filter(x -> x ∉ ("FID", "IID"), names(y_df)))
    X = data[!, variant_cols] |> Matrix{Float64}
    y = data[!, phenotype_col] |> Vector{Float64}
    return X, y
end

function susie_finemap(X, y; n_causal=10)
    @rput X
    @rput y
    @rput n_causal
    R"""
    library(susieR)
    fitted = susie(X, y, L=n_causal)
    """
    return @rget fitted
end

function get_r2_with_variant(variant_id, pgen_prefix; ld_window_kb=1000)
    tmpdir = mktempdir()
    out_prefix = joinpath(tmpdir, "$variant_id.LD")
    run(`plink2 \
        --pfile $pgen_prefix \
        --r2-unphased \
        --ld-snp $variant_id \
        --ld-window-kb $ld_window_kb \
        --ld-window-r2 0 \
        --out $out_prefix
    `)
    variants_r2 = CSV.read(out_prefix * ".vcor", DataFrame, delim='\t')
    # add the lead variant with R2=1.0 to itself
    row = variants_r2[1, :]
    push!(variants_r2, [row["#CHROM_A"], row.POS_A, row.ID_A, row["#CHROM_A"], row.POS_A, row.ID_A, 1.0])
    return sort!(variants_r2, :POS_B)
end

update_credible_sets!(cs_vector, variant_idx::Int, cs::Symbol) = cs_vector[variant_idx] = parse(Int, replace(string(cs), "L" => ""))

function update_credible_sets!(cs_vector, variant_idxs::AbstractVector, cs::Symbol)
    for variant_idx in variant_idxs
        update_credible_sets!(cs_vector, variant_idx, cs)
    end
end

function update_credible_sets!(cs_vector, susie_results)
    for (cs, variant_idxs) in susie_results[:sets][:cs]
        update_credible_sets!(cs_vector, variant_idxs, cs)
    end
end

function get_credible_sets(susie_results, p)
    cs_vector = Vector{Union{Missing, Int}}(undef, p)
    susie_results[:sets][:cs] isa Nothing && return cs_vector
    update_credible_sets!(cs_vector, susie_results)
    return cs_vector
end

function postprocess_finemapping_results!(variants_info, finemapping_results, ld_variants)
    locus_id = first(ld_variants.ID_A)
    variants_info.PIP = finemapping_results[:pip]
    p = length(variants_info.PIP)
    variants_info.LOCUS_ID = fill(locus_id, p)
    variants_info.SUSIE_CONVERGED = fill(finemapping_results[:converged], p)
    variants_info.CS = get_credible_sets(finemapping_results, p)
    leftjoin!(variants_info, ld_variants[!, [:ID_B, :UNPHASED_R2]], on=:ID => :ID_B)
    return select!(variants_info, :CHROM, :POS, :ID, :REF, :ALT, :LOCUS_ID, :PIP, :CS, :UNPHASED_R2, :SUSIE_CONVERGED)
end

function finemap_locus(variant_ID, pgen_prefix, y_df;
    n_causal=10,
    finemap_window_kb=1000,
    r2_threshold=0.1,
    padding=10
    )
    r2s_with_lead = get_r2_with_variant(variant_ID, pgen_prefix; ld_window_kb=finemap_window_kb)
    locus = locus_from_r2s(r2s_with_lead, variant_ID; r2_threshold=r2_threshold, padding=padding)
    X_df, variants_info = genotypes_from_pgen(pgen_prefix, locus)
    X, y = get_susie_inputs(X_df, y_df)
    finemapping_results = susie_finemap(X, y; n_causal=n_causal)
    return postprocess_finemapping_results!(variants_info, finemapping_results, r2s_with_lead)
end

function get_LD_matrix(pgen_prefix, locus)
    _, chr, locus_start, locus_end = locus
    tmpdir = mktempdir()
    output_prefix = joinpath(tmpdir, "ld_matrix")
    run(`plink2 \
        --pfile $pgen_prefix \
        --r-unphased square \
        --chr $chr \
        --from-bp $locus_start \
        --to-bp $locus_end \
        --out $output_prefix`
    )
    R = readdlm(string(output_prefix, ".unphased.vcor1"))
    variants = readlines(string(output_prefix, ".unphased.vcor1.vars"))
    rm(tmpdir, recursive=true)
    return R, variants
end

function susie_rss_finemap(R, variants_info, y; n_causal=10)
    var_y, nsamples = var(y), length(y)
    shat = variants_info.SE
    bhat = variants_info.BETA
    @rput R
    @rput shat
    @rput bhat
    @rput var_y
    @rput n_causal
    @rput nsamples
    R"""
    library(susieR)
    fitted = susie_rss(R=R, bhat=bhat, shat=shat, var_y=var_y, L=n_causal, n=nsamples)
    """
    return @rget fitted
end

function initialize_variants_info_rss(pgen_prefix, variants, gwas_results)
    pvar = read_pvar(string(pgen_prefix, ".pvar"))
    variants_info = leftjoin!(
        leftjoin!(
            DataFrame(ID=variants),
            pvar,
            on=:ID
        ),
        gwas_results,
        on=:ID
    )
    
    return select!(variants_info, :CHROM, :POS, :ID, :REF, :ALT, :BETA, :SE)
end

function finemap_locus_rss(variant_ID, gwas_results, pgen_prefix, y;
                n_causal=n_causal,
                finemap_window_kb=finemap_window_kb,
                r2_threshold=0.1,
                padding=10
            )
    r2s_with_lead = get_r2_with_variant(variant_ID, pgen_prefix; ld_window_kb=finemap_window_kb)
    locus = locus_from_r2s(r2s_with_lead, variant_ID; r2_threshold=r2_threshold, padding=padding)
    R, variants = get_LD_matrix(pgen_prefix, locus)
    variants_info = initialize_variants_info_rss(pgen_prefix, variants, gwas_results)
    finemapping_results = susie_rss_finemap(R, variants_info, y; n_causal=n_causal)
    return postprocess_finemapping_results!(variants_info, finemapping_results, r2s_with_lead)
end

function make_clean_sample_file(input_sample_file; exclude=[], phenotype="Y")
    input_lines = readlines(input_sample_file)
    if isfile(first(input_lines))
        tmpdir = mktempdir()
        output_sample_file = joinpath(tmpdir, "clean_samples.txt")
        open(output_sample_file, "w") do output_io
            for file in input_lines
                _, _, file_group, file_pheno, _ = split(basename(file), '.')
                group_needs_exclusion(file_group, exclude) && continue
                phenotype == file_pheno || continue
                for line in readlines(file)
                    println(output_io, line)
                end
            end
        end
        return output_sample_file
    else
        return input_sample_file
    end
end

"""
The sample_file is assumed to be a tab-delimited file with two columns: FID and IID. 
"""
function read_samples(sample_file)
    return CSV.read(sample_file, DataFrame, header=["FID", "IID"], delim="\t")
end

"""
    finemap_significant_regions(
        gwas_results_file,
        pgen_prefix,
        covariates_file,
        sample_file;
        output_prefix = "clumps.sig.tsv",
        min_sig_clump_size = 3,
        lead_pvalue = 5e-8,
        p2_pvalue = 5e-5,
        r2_threshold = 0.5,
        clump_kb = 250,
        clump_id_field = "ID",
        clump_pval_field = "LOG10P",
        allele_1_field = "ALLELE_1",
        finemap_window_kb=1000,
        n_causal = 10
        )

This function performs fine-mapping of significant regions identified from GWAS results.

1. Identifies variants from the PGEN fileset that have gone through GWAS (some may not due to MAF/MAC/...)
2. Identify clumps of significant variants from the GWAS results to tag independent association regions
3. Finemap these regions using SuSiE

# Arguments
- `gwas_results_file::String`: Path to the GWAS results file (TSV format).
- `pgen_prefix::String`: Prefix for the PGEN fileset (without .pgen extension).
- `covariates_file::String`: Path to the covariates file (TSV format).
- `sample_file::String`: Path to the sample IDs file used to generate the GWAS results.
- `output_prefix::String`: Prefix to output the significant clumps (TSV format).
- `min_sig_clump_size::Int`: Minimum number of variants in a clump to be considered significant.
- `lead_pvalue::Float64`: P-value threshold for lead variants in clump.
- `p2_pvalue::Float64`: Secondary p-value threshold for variants in clump.
- `r2_threshold::Float64`: LD r² threshold for clumping.
- `clump_kb::Int`: Distance in kb for clumping.
- `clump_id_field::String`: Column name in GWAS results for variant IDs.
- `clump_pval_field::String`: Column name in GWAS results for p-values.
- `allele_1_field::String`: Column name in GWAS results for allele 1.
- `finemap_window_kb::Int`: Window size in kb for LD calculation around lead variant to define the fine mapping region.
- `n_causal::Int`: Number of causal variants to assume in SuSiE fine-mapping.
- `phenotype::String`: Name of the phenotype column in the covariates file.
- `rss::Bool`: Whether to use summary statistics fine-mapping.
- `exclude_string::String`: Comma-separated list of group name patterns to exclude from the sample file.
"""
function finemap_significant_regions(
    gwas_results_file,
    pgen_prefix,
    covariates_file,
    sample_file;
    output_prefix = "finemapping_results",
    min_sig_clump_size = 3,
    lead_pvalue = 5e-8,
    p2_pvalue = 5e-5,
    r2_threshold = 0.5,
    clump_kb = 250,
    clump_id_field = "ID",
    clump_pval_field = "LOG10P",
    allele_1_field = "ALLELE_1",
    finemap_window_kb=1000,
    n_causal = 10,
    phenotype="Y",
    rss=false,
    exclude_string="",
    padding=10
    )
    # Read GWAS results, PVAR files and samples ids
    @info "Reading GWAS results and PVAR files"
    gwas_results = CSV.read(gwas_results_file, DataFrame)
    pvar = read_pvar(pgen_prefix * ".pvar")
    sample_file = make_clean_sample_file(sample_file; 
        exclude=split(exclude_string, ","), 
        phenotype=phenotype
    )
    # Tag variant IDs that are not in the GWAS results
    @info "Tagging variants not in GWAS results"
    tag_variant_id_missing_from_gwas!(pvar, gwas_results)

    # Make new PGEN fileset with variants from GWAS only
    @info "Writing new PGEN fileset with variants matched to GWAS results"
    tmpdir = mktempdir()
    gwas_matched_pgen_prefix = joinpath(tmpdir, "gwas_matched")
    write_new_pgen_from_gwas_results(pgen_prefix, gwas_matched_pgen_prefix, pvar, sample_file)

    # Find clumps
    @info "Finding clumps in GWAS-matched PGEN fileset"
    sig_clumps = write_significant_clumps(gwas_matched_pgen_prefix, gwas_results_file;
        min_sig_clump_size = min_sig_clump_size,
        output = string(output_prefix, ".clumps.tsv"),
        lead_pvalue = lead_pvalue,
        p2_pvalue = p2_pvalue,
        r2_threshold = r2_threshold,
        clump_kb = clump_kb,
        clump_id_field = clump_id_field,
        clump_pval_field = clump_pval_field,
        allele_1_field = allele_1_field
    )

    # Finemap each clump
    y_df = load_phenotypes_matching_samples(covariates_file, sample_file, phenotype)
    finemapping_results = []
    for clump in eachrow(sig_clumps)
        @info "Fine Mapping locus led by : $(clump.ID)"
        clump_finemapping_results = if rss
            finemap_locus_rss(clump.ID, gwas_results, gwas_matched_pgen_prefix, y_df[!, phenotype];
                n_causal=n_causal,
                finemap_window_kb=finemap_window_kb,
                r2_threshold=r2_threshold,
                padding=padding
            )
        else
            finemap_locus(clump.ID, gwas_matched_pgen_prefix, y_df;
                n_causal=n_causal,
                finemap_window_kb=finemap_window_kb,
                r2_threshold=r2_threshold,
                padding=padding
            )
        end
        push!(finemapping_results, clump_finemapping_results)
    end
    output_file = string(output_prefix, ".tsv")
    output_df = length(finemapping_results) > 0 ? vcat(finemapping_results...) : DataFrame([col => [] for col in FINEMAPPING_RESULT_COLS])
    CSV.write(output_file, output_df, delim="\t")

    return 0
end

