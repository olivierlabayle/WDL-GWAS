function group_and_phenotype_from_regenie_filename(filename)
    return split(replace(filename, 
        "results.all_chr." => "", 
        ".gwas.tsv" => ""
    ), ".")
end

function run_metal_across_phenotypes!(regenie_files; output_prefix="gwas.meta_analysis", method="STDERR")
    tmp_dir = mktempdir()
    for (phenotype_key, group) in pairs(groupby(regenie_files, :PHENOTYPE))
        phenotype = phenotype_key.PHENOTYPE
        metal_script = """
        # === DESCRIBE THE COLUMNS IN THE INPUT FILES ===
        MARKER ID 
        GENOMICCONTROL ON
        WEIGHT N 
        ALLELE ALLELE_1 ALLELE_0 
        FREQ ALLELE_1_FREQ 
        EFFECT BETA 
        STDERR SE 
        PVAL P_VAL
        SCHEME $method
        LOGPVALUE ON
        # === FOR EACH PHENOTYPE PROCESS / ANALYZE / ANALYZE HETEROGENEITY ===
        """
        for (group_file, group_basename) in zip(group.FILE, group.BASENAME)
            # Add P_VAL column expected by METAL
            group_gwas_results = CSV.read(group_file, DataFrame; 
                delim="\t", 
                missingstring="NA",
                select=[:ID, :ALLELE_0, :ALLELE_1, :ALLELE_1_FREQ, :BETA, :LOG10P, :SE, :N]
            )
            transform!(group_gwas_results, :LOG10P => (x -> neg_exp10.(x))  => :P_VAL)
            CSV.write(joinpath(tmp_dir, group_basename), group_gwas_results; 
                delim="\t", 
                header=true, 
                missingstring="NA"
            )
            # Add PROCESS command
            metal_script *= "PROCESS " * joinpath(tmp_dir, group_basename) * "\n"
        end
        metal_script *= "OUTFILE " * string(output_prefix, ".", phenotype, ". .tbl") * "\n"
        metal_script *= "ANALYZE HETEROGENEITY\n"
        group.METAL_FILE .= string(output_prefix, ".", phenotype, ".1.tbl")
        metal_script *= "QUIT"
        meta_script_file = joinpath(tmp_dir, "metal_script.$phenotype.txt")
        open(meta_script_file, "w") do io
            write(io, metal_script)
        end
        run(`metal $meta_script_file`)
    end
    return regenie_files
end

function post_process_metal_output(regenie_files; output_prefix="gwas.meta_analysis")
    for (phenotype_key, group) in pairs(groupby(regenie_files, :PHENOTYPE))
        metal_results = CSV.read(first(group.METAL_FILE), DataFrame; delim="\t")
        select!(metal_results, 
            "MarkerName" => "ID",
            "Effect" => "BETA",
            "StdErr" => "SE",
            "log(P)" => (x -> .-x) => "LOG10P", # -log10(P) is reported by Regenie, we make it compliant
            "Direction" => "DIRECTION",
            "HetISq" => "HET_ISQ",
            "HetChiSq" => "HET_CHISQ",
            "HetDf" => "HET_DF",
            "logHetP" => "LOG10P_HET"
        )
        append_GWAS_info_to_meta_analysis_results!(metal_results, group.FILE)
        CSV.write(string(output_prefix, ".", phenotype_key.PHENOTYPE, ".gwas.tsv"), metal_results; 
            delim="\t", 
            header=true,
            missingstring="NA"
            )
    end
end

function append_GWAS_info_to_meta_analysis_results!(metal_results, phenotype_gwas_files)
    # Get variant info from GWAS results: CHROM, POS, ALLELE_0, ALLELE_1, ALLELE_1_FREQ, N, NGROUPS
    variants_info_dict = Dict{String, Vector{Any}}()
    for phenotype_gwas_file in phenotype_gwas_files
        phenotype_gwas_results = CSV.read(phenotype_gwas_file, DataFrame; delim="\t", missingstring="NA")
        for row in Tables.namedtupleiterator(phenotype_gwas_results[!, [:CHROM, :POS, :ID, :ALLELE_0, :ALLELE_1, :ALLELE_1_FREQ, :N]])
            if haskey(variants_info_dict, row.ID)
                variant_info = variants_info_dict[row.ID]
                variant_info[end] += 1 # update count of groups the variant was observed in
                variant_info[end-1] += row.N # update total N
                variant_info[end-2] = min(variant_info[end-2], row.ALLELE_1_FREQ) # update min ALLELE_1_FREQ
            else
                variants_info_dict[row.ID] = [row.CHROM, row.POS, row.ALLELE_0, row.ALLELE_1, row.ALLELE_1_FREQ, row.N, 1]
            end
        end
    end
    # Update metal results with variant info
    transform!(metal_results, 
        :ID => ByRow(id -> variants_info_dict[id]) => [:CHROM, :POS, :ALLELE_0, :ALLELE_1, :ALLELE_1_FREQ, :N, :NGROUPS]
    )

    return metal_results
end

function load_meta_analysis_worklist(regenie_files_list; exclude = [])
    regenie_files = DataFrame(FILE = readlines(regenie_files_list))
    regenie_files.BASENAME = basename.(regenie_files.FILE)
    transform!(regenie_files, 
        :BASENAME => ByRow(group_and_phenotype_from_regenie_filename) 
        => [:GROUP, :PHENOTYPE]
    )
    return filter!(:GROUP => (group -> !group_needs_exclusion(group, exclude)), regenie_files)
end

"""
    meta_analyse(regenie_files_list; exclude_string="ADMIXED", method="STDERR", output_prefix="gwas.meta_analysis")

Meta-analyse GWAS results from REGENIE using METAL. Groups and phenotypes are inferred from filenames. Results are meta analysed per phenotype across groups.

- regenie_files_list: path to a text file with a list of GWAS result files to meta-analyse
- exclude_string: comma-separated list of strings, if a group name contains any of these strings it will be excluded from meta-analysis (default: "ADMIXED")
- method: METAL meta-analysis method (default: "STDERR")
- output_prefix: prefix for output files. Per phenotype results are written to "<output_prefix>.<phenotype>.gwas.tsv". (default: "gwas.meta_analysis")
"""
function meta_analyse(regenie_files_list; exclude_string="ADMIXED", method="STDERR", output_prefix="gwas.meta_analysis")
    exclude = split(exclude_string, ",")
    regenie_files = load_meta_analysis_worklist(regenie_files_list; exclude = exclude)
    run_metal_across_phenotypes!(regenie_files; output_prefix=output_prefix, method=method)
    post_process_metal_output(regenie_files, output_prefix=output_prefix)

    return 0
end