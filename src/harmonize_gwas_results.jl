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
                Symbol("p.value") => (x -> -log10.(x)) => :LOG10P,
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
                Symbol("p.value") => (x -> -log10.(x)) => :LOG10P,
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
    else
        throw(ArgumentError("GWAS software $source_software is not supported."))
    end
    CSV.write(output, harmonized_results, delim="\t", missingstring="NA")

    return 0
end