function merge_chr_results(merge_list_file; output_prefix = "results.all_chr")
    merge_list = readlines(merge_list_file)
    results = mapreduce(f -> CSV.read(f, DataFrame; missingstring="NA"), vcat, merge_list)
    CSV.write(string(output_prefix, ".tsv"), results; 
        delim="\t", 
        header=true,
        missingstring="NA"
    )
    return 0
end

get_chr_out_string(pc_filename) = splitext(splitext(pc_filename)[1])[2][2:end]

function read_loco_pcs(pc_file)
    chr_out = get_chr_out_string(pc_file)
    pcs = CSV.read(pc_file, DataFrame, drop=["#FID"])
    PC_colnames = filter(!=("IID"), names(pcs))
    for PC_colname in PC_colnames
        rename!(pcs, Symbol(PC_colname) => Symbol(string(uppercase(chr_out), "_", PC_colname)))
    end
    return pcs
end

function merge_covariates_and_pcs(covariates_file, pcs_prefix; output="covariates_and_pcs.csv")
    covariates = CSV.read(covariates_file, DataFrame)
    pcs_dir = dirname(pcs_prefix)
    pcs_dir, pcs_prefix = pcs_dir == "" ? (".", "./$pcs_prefix") : (pcs_dir, pcs_prefix)
    pcs_files = filter(startswith(pcs_prefix), readdir(pcs_dir, join=true))
    chrs = unique(get_chr_out_string.(pcs_files))
    chr_out_pcs = map(chrs) do chr
        chr_out_pcs_files = filter(x -> occursin(chr, x), pcs_files)
        mapreduce(read_loco_pcs, vcat, chr_out_pcs_files)
    end
    covariates_and_pcs = innerjoin(covariates, chr_out_pcs..., on=:IID)
    CSV.write(output, covariates_and_pcs, delim="\t", missingstring="NA")
    return 0
end

group_needs_exclusion(group, exclude) =
    any(occursin(p, group) for p in exclude)

read_plink2_df(plink_file; delim='\t', comment="##") =
    CSV.read(plink_file, DataFrame; delim=delim, comment=comment)

read_pvar(pvar_file; delim='\t', comment="##") =
    read_plink2_df(pvar_file; delim=delim, comment=comment)

read_psam(psam_file; delim='\t', comment="##") =
    read_plink2_df(psam_file; delim=delim, comment=comment)
