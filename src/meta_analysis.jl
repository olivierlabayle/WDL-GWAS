function group_and_phenotype_from_regenie_filename(filename)
    return split(replace(filename, 
        "results.all_chr." => "", 
        ".gwas.tsv" => ""
    ), ".")
end

function split_marker_id(marker_id) 
    chr, pos, all0, all1 = split(marker_id, ":")
    chr = parse(Int, replace(chr, "chr" => ""))
    pos = parse(Int, pos)
    return chr, pos, all0, all1
end

function run_metal_across_phenotypes!(regenie_files; output_prefix="gwas.meta_analysis", method="STDERR", maf=0.01)
    tmp_dir = mktempdir()
    for (phenotype_key, group) in pairs(groupby(regenie_files, :PHENOTYPE))
        phenotype = phenotype_key.PHENOTYPE
        metal_script = """
        # === DESCRIBE THE COLUMNS IN THE INPUT FILES ===
        CUSTOMVARIABLE TotalSampleSize
        AVERAGEFREQ ON
        MINMAXFREQ ON
        MARKER ID 
        GENOMICCONTROL ON
        WEIGHT N 
        ALLELE ALLELE_1 ALLELE_0 
        FREQ ALLELE_1_FREQ
        LABEL TotalSampleSize as N
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
            subset!(group_gwas_results, :ALLELE_1_FREQ => x -> maf .< x .< 1 - maf, skipmissing=true)
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

function post_process_metal_output(regenie_files; maf=0.01, output_prefix="gwas.meta_analysis")
    for (phenotype_key, group) in pairs(groupby(regenie_files, :PHENOTYPE))
        metal_results = CSV.read(first(group.METAL_FILE), DataFrame; delim="\t")
        # harmonize names
        select!(metal_results, 
            "MarkerName" => "ID",
            "MarkerName" => ByRow(x -> split_marker_id(x)) => ["CHROM", "POS", "ALLELE_0", "ALLELE_1"],
            "Effect" => "BETA",
            "StdErr" => "SE",
            "log(P)" => (x -> .-x) => "LOG10P",
            "Freq1" => "ALLELE_1_FREQ",
            "FreqSE" => "ALLELE_1_FREQ_STD",
            "MinFreq" => "ALLELE_1_FREQ_MIN",
            "MaxFreq" => "ALLELE_1_FREQ_MAX",
            "Direction" => "DIRECTION",
            "HetISq" => "HET_ISQ",
            "HetChiSq" => "HET_CHISQ",
            "HetDf" => "HET_DF",
            "logHetP" => "LOG10P_HET",
            "TotalSampleSize" => "N",
            "Direction" => (col -> count.(∉(Set(['?', '0'])), col)) => "NGROUPS"
        )
        # Write output file
        phenotype_prefix = string(output_prefix, ".", phenotype_key.PHENOTYPE)
        CSV.write(string(phenotype_prefix, ".gwas.tsv"), metal_results; 
            delim="\t", 
            header=true,
            missingstring="NA"
        )
        # Make plots
        make_gwas_plots(metal_results; maf=maf, output_prefix=phenotype_prefix)
    end
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
function meta_analyse(regenie_files_list; maf=0.01, exclude_string="ADMIXED", method="STDERR", output_prefix="gwas.meta_analysis")
    exclude = split(exclude_string, ",")
    regenie_files = load_meta_analysis_worklist(regenie_files_list; exclude = exclude)
    run_metal_across_phenotypes!(regenie_files; output_prefix=output_prefix, method=method, maf=maf)
    post_process_metal_output(regenie_files; maf=maf, output_prefix=output_prefix)

    return 0
end