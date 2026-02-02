function safe_pvalue_to_log10p(pvalue)
    log10pval = abs(log10(pvalue))
    if log10pval === Inf
        log10pval = 300
    end
    return log10pval
end

function postprocess_plink2_binary_result(e, or, se, lp)
    missing_result = (missing, missing, missing)
    return if e != "."
        missing_result
    elseif any(ismissing, (or, se, lp))
        missing_result
    else
        result = (log(or), se, lp)
        any(isinf, result) ? missing_result : result
    end
end

function postprocess_plink2_continuous_result(e, beta, se, lp)
    missing_result = (missing, missing, missing)
    result = (beta, se, lp)
    return if e != "."
        missing_result
    elseif any(ismissing, result)
        missing_result
    elseif any(isinf, result)
        missing_result
    else
        result
    end
end

function harmonize_gwas_results(gwas_results_file; source_software="saige", output="harmonized_results.tsv")
    gwas_results = CSV.read(gwas_results_file, DataFrame; missingstring="NA")
    harmonized_results = if source_software == "saige"
        # Binary trait
        if "N_case" in names(gwas_results)
            DataFrames.select(gwas_results,
                :CHR => :CHROM,
                :POS,
                :MarkerID => :ID,
                :Allele1 => :ALLELE_0,
                :Allele2 => :ALLELE_1,
                :AF_Allele2 => :ALLELE_1_FREQ,
                :BETA,
                :SE,
                Symbol("p.value") => (x -> safe_pvalue_to_log10p.(x)) => :LOG10P,
                [:N_case, :N_ctrl] => ((nca, nco) -> nca .+ nco) => :N,
                :AC_Allele2 => :ALLELE_1_COUNT,
                :MissingRate => :MISSING_RATE,
                :Tstat => :T_STAT,
                :var => :VAR,
                :N_case => :N_CASES,
                :N_ctrl => :N_CONTROLS,
                Symbol("p.value.NA") => :PVAL_NA,
                Symbol("Is.SPA") => :IS_SPA,
                :AF_case => :AF_CASES,
                :AF_ctrl => :AF_CONTROLS,
                :N_case_hom => :N_CASES_HOM,
                :N_case_het => :N_CASES_HET,
                :N_ctrl_hom => :N_CONTROLS_HOM,
                :N_ctrl_het => :N_CONTROLS_HET
            )
        # Continuous trait
        else
            DataFrames.select(gwas_results,
                :CHR => :CHROM,
                :POS,
                :MarkerID => :ID,
                :Allele1 => :ALLELE_0,
                :Allele2 => :ALLELE_1,
                :AF_Allele2  => :ALLELE_1_FREQ,
                :BETA,
                :SE,
                Symbol("p.value") => (x -> safe_pvalue_to_log10p.(x)) => :LOG10P,
                :N,
                :AC_Allele2 => :ALLELE_1_COUNT,
                :MissingRate => :MISSING_RATE,
                :Tstat => :T_STAT,
                :var => :VAR,
            )
        end
    elseif source_software == "regenie"
        DataFrames.select(gwas_results,
            :CHROM,
            :GENPOS => :POS,
            :ID,
            :ALLELE0 => :ALLELE_0,
            :ALLELE1 => :ALLELE_1,
            :A1FREQ  => :ALLELE_1_FREQ,
            :BETA,
            :SE,
            :LOG10P,
            :N,
            :TEST,
            :CHISQ,
            :EXTRA
        )
    elseif source_software == "plink2"
        # Binary trait
        if "OR" in names(gwas_results)
            DataFrames.select(gwas_results,
                Symbol("#CHROM") => :CHROM,
                :POS,
                :ID,
                :OMITTED => :ALLELE_0,
                :A1 => :ALLELE_1,
                :A1_FREQ  => :ALLELE_1_FREQ,
                [:ERRCODE, :OR, Symbol("LOG(OR)_SE"), :NEG_LOG10_P] => ByRow(postprocess_plink2_binary_result) => [:BETA, :SE, :LOG10P],
                :OBS_CT => :N,
                :TEST,
                Symbol("FIRTH?") => :FIRTH,
                :ERRCODE,
                :Z_STAT
            )
        # Continuous trait
        else
            DataFrames.select(gwas_results,
                Symbol("#CHROM") => :CHROM,
                :POS,
                :ID,
                :OMITTED => :ALLELE_0,
                :A1 => :ALLELE_1,
                :A1_FREQ  => :ALLELE_1_FREQ,
                [:ERRCODE, :BETA, :SE, :NEG_LOG10_P] => ByRow(postprocess_plink2_continuous_result) => [:BETA, :SE, :LOG10P],
                :OBS_CT => :N,
                :TEST,
                :ERRCODE,
                :T_STAT
            )
        end
    else
        throw(ArgumentError("GWAS software $source_software is not supported."))
    end
    CSV.write(output, harmonized_results, delim="\t", missingstring="NA")

    return 0
end