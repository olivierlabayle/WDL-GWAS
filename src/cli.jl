function cli_settings()
    s = ArgParseSettings(
        description="WDL-GWAS CLI",
        add_version = true,
        commands_are_required = false,
        version=string(pkgversion(PopGen))
    )

    @add_arg_table! s begin
        "meta-analyse"
            action = :command
            help = "Runs meta-analysis across GWAS results."

        "finemap"
            action = :command
            help = "Runs fine-mapping analysis."

        "gwas-plots"
            action = :command
            help = "Generates GWAS plots."

        "pca-qc"
            action = :command
            help = "Runs PCAb-ased QC on genotypes to exclude outlier variants."

        "make-gwas-groups"
            action = :command
            help = "Generates groups, phenotypes and covariates files."
        
        "merge-covariates-pcs"
            action = :command
            help = "Merges covariates and PCs files."

        "merge-chr-results"
            action = :command
            help = "Merges REGENIE results from different chromosomes."
    end

    @add_arg_table! s["meta-analyse"] begin
        "gwas-results-list"
            arg_type = String
            required = true
            help = "List of REGENIE results files to be meta-analysed."
    
        "--exclude"
            arg_type = String
            help = "List of groups to exclude from the meta-analysis, comma separated."

        "--method"
            arg_type = String
            help = "Meta-analysis method to use (fixed or random)."
            default = "STDERR"

        "--output-prefix"
            arg_type = String
            help = "Prefix to output files."
            default = "gwas.meta_analysis"
    end

    @add_arg_table! s["finemap"] begin
        "gwas-results-file"
            arg_type = String
            required = true
            help = "Path to the GWAS results file."
        "pgen-prefix"
            arg_type = String
            required = true
            help = "Prefix for the PGEN fileset (without .pgen extension)."
        "covariates-file"
            arg_type = String
            required = true
            help = "Path to the covariates file (TSV format)."
        "sample-file"
            arg_type = String
            required = true
            help = "Path to the sample IDs file used to generate the GWAS results."
        "--Xtype"
            arg_type = String
            help = "Type of genotype data to use for fine-mapping, either `:dosages` or `:genotypes`."
            default = "dosages"
        "--output-prefix"
            arg_type = String
            help = "Prefix to output the significant clumps (TSV format)."
            default = "finemapping_results"
        "--min-sig-clump-size"
            arg_type = Int
            help = "Minimum number of variants in a clump to be considered significant."
            default = 3
        "--lead-pvalue"
            arg_type = Float64
            help = "P-value threshold for lead variants in clump."
            default = 5e-8
        "--p2-pvalue"
            arg_type = Float64
            help = "P-value threshold for secondary variants in clump."
            default = 1e-5
        "--r2-threshold"
            arg_type = Float64
            help = "R2 threshold to consider variants in LD for clumping."
            default = 0.5
        "--clump-kb"
            arg_type = Int
            help = "Window size (in kb) to consider variants in LD for clumping."
            default = 250
        "--n-causal"
            arg_type = Int
            help = "Number of causal variants to assume in fine-mapping."
            default = 10
        "--finemap-window-kb"
            arg_type = Int
            help = "Window size (in kb) to compute LD matrix for fine-mapping."
            default = 1000
        "--phenotype"
            arg_type = String
            help = "Name of the phenotype column in the covariates file."
            default = "Y"
        "--rss"
            help = "Whether to run fine-mapping using summary statistics (RSS)"
            action = :store_true
        "--exclude"
            arg_type = String
            help = "List of groups to exclude from the sample file, comma separated."
            default = ""
    end

    @add_arg_table! s["merge-chr-results"] begin
        "merge-list"
            arg_type = String
            required = true
            help = "File containing list of files to be merged."
        "--output-prefix"
            arg_type = String
            help = "Output path"
            default = "results.all_chr"
    end

    @add_arg_table! s["gwas-plots"] begin
        "gwas-results"
            arg_type = String
            required = true
            help = "Path to GWAS results file."

        "finemapping-results"
            arg_type = String
            required = true
            help = "Path to finemapping results file."

        "--maf"
            arg_type = Float64
            help = "Minor allele frequency threshold to filter results."
            default = 0.01
        
        "--output-prefix"
            arg_type = String
            help = "Prefix to output files."
            default = "gwas.plot"
    end

    @add_arg_table! s["merge-covariates-pcs"] begin
        "covariates-file"
            arg_type = String
            required = true
            help = "Path to covariates file."
        
        "pcs-prefix"
            arg_type = String
            required = true
            help = "Prefix to PCs files."

        "--output"
            arg_type = String
            help = "Output file name."
            default = "covariates_and_pcs.csv"
    end

    @add_arg_table! s["make-gwas-groups"] begin
        "covariates-file"
            arg_type = String
            required = true
            help = "Path to covariates file."

        "--groupby"
            arg_type = String
            help = "Comma separated list of variables to use to stratify the GWAS."
            default = nothing

        "--filters"
            arg_type = String
            help = "Filters to apply to the data."
            default = nothing

        "--covariates"
            arg_type = String
            help = "Comma separated list of covariates to include in the output file."
            default = "AGE"

        "--phenotypes"
            arg_type = String
            help = "Comma separated list of phenotypes to include in the output file."
            default = "SEVERE_COVID_19"

        "--output-prefix"
            arg_type = String
            help = "Prefix to output files."
            default = "group"
        
        "--min-cases-controls"
            arg_type = Int
            help = "Minimum group size."
            default = 100
    end

    return s
end

function julia_main()::Cint
    settings = parse_args(ARGS, cli_settings())
    cmd = settings["%COMMAND%"]
    @info "Running GenOMICC Workflows: $cmd"
    cmd_settings = settings[cmd]
    if cmd == "make-gwas-groups"
        make_gwas_groups(
            cmd_settings["covariates-file"];
            groupby_string=cmd_settings["groupby"],
            covariates_string=cmd_settings["covariates"],
            phenotypes_string=cmd_settings["phenotypes"],
            output_prefix=cmd_settings["output-prefix"],
            min_cases_controls=cmd_settings["min-cases-controls"],
            filters_string=cmd_settings["filters"]
        )
    elseif cmd == "merge-covariates-pcs"
        merge_covariates_and_pcs(
            cmd_settings["covariates-file"],
            cmd_settings["pcs-prefix"];
            output=cmd_settings["output"]
        )
    elseif cmd == "gwas-plots"
        gwas_plots(
            cmd_settings["gwas-results"],
            cmd_settings["finemapping-results"];
            maf=cmd_settings["maf"],
            output_prefix=cmd_settings["output-prefix"]
        )
    elseif cmd == "merge-chr-results"
        merge_chr_results(
            cmd_settings["merge-list"];
            output_prefix=cmd_settings["output-prefix"]
        )
    elseif cmd == "finemap"
        finemap_significant_regions(
            cmd_settings["gwas-results-file"],
            cmd_settings["pgen-prefix"],
            cmd_settings["covariates-file"],
            cmd_settings["sample-file"];
            Xtype=cmd_settings["Xtype"],
            output_prefix=cmd_settings["output-prefix"],
            min_sig_clump_size=cmd_settings["min-sig-clump-size"],
            lead_pvalue=cmd_settings["lead-pvalue"],
            p2_pvalue=cmd_settings["p2-pvalue"],
            r2_threshold=cmd_settings["r2-threshold"],
            clump_kb=cmd_settings["clump-kb"],
            n_causal=cmd_settings["n-causal"],
            finemap_window_kb=cmd_settings["finemap-window-kb"],
            phenotype=cmd_settings["phenotype"],
            rss=cmd_settings["rss"],
            exclude_string=cmd_settings["exclude"]
        )
    elseif cmd == "meta-analyse"
        meta_analyse(
            cmd_settings["gwas-results-list"];
            output_prefix=cmd_settings["output-prefix"],
            exclude_string=cmd_settings["exclude"],
            method=cmd_settings["method"],
        )
    else
        throw(ArgumentError(string("Unknown command: ", cmd)))
    end
    return 0
end