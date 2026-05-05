function add_user_defined_covariates!(covariates, covariates_string)
    required_covariate_variables = split(covariates_string, ",")
    all_colnames = names(covariates)
    updated_required_covariate_variables = []
    for variable in required_covariate_variables
        if variable ∈ all_colnames
            if eltype(covariates[!, variable]) <: AbstractString
                covariates[!, variable] = categorical(covariates[!, variable])
                mach = machine(OneHotEncoder(drop_last=true), covariates[!, [variable]]) ## Need to drop last to avoid multicolinearity in GWAS
                fit!(mach, verbosity=0)
                Xt = MLJBase.transform(mach)
                for colname in names(Xt)
                    covariates[!, colname] = Xt[!, colname]
                    push!(updated_required_covariate_variables, colname)
                end
            else
                push!(updated_required_covariate_variables, variable)
            end
        elseif occursin("_x_", variable)
            base_variables = split(variable, "_x_")
            issubset(base_variables, all_colnames) || throw(ArgumentError("Some base covariate in $variable was not found."))
            columns = (covariates[!, v] for v in base_variables)
            covariates[!, variable] = .*(columns...)
            push!(updated_required_covariate_variables, variable)
        else
            throw(ArgumentError("Covariate $variable not suported."))
        end
    end
    return updated_required_covariate_variables
end

apply_filters(data, ::Nothing) = data

function apply_filters(data, filters_string::AbstractString)
    filters_strings = split(filters_string, ",")
    for filter_string in filters_strings
        operator = match(r"(<=|>=|<|>|=)", filter_string).match
        col, value = split(filter_string, operator)
        value = if eltype(data[!, col]) <: Union{Missing,Number}
                parse(Float64, value)
            else
                @assert value in data[!, col] "Filtering value $value not found in column $col."
                value
        end
        # There is likely a better way to do this with metaprogramming but can't findout quickly..
        if operator == "="
            data = data[(data[!, col] .!== missing) .& (data[!, col] .== value), :]
        elseif operator == "<"
            data = data[(data[!, col] .!== missing) .& (data[!, col] .< value), :]
        elseif operator == ">"
            data = data[(data[!, col] .!== missing) .& (data[!, col] .> value), :]
        elseif operator == "<="
            data = data[(data[!, col] .!== missing) .& (data[!, col] .<= value), :]
        elseif operator == ">="
            data = data[(data[!, col] .!== missing) .& (data[!, col] .>= value), :]
        else
            throw(ArgumentError("Operator $operator not supported."))
        end
    end
    return data
end

is_binary_column(data, column) =
    issubset(Set(skipmissing(data[!, column])), Set([0, 1]))

function n_cases_controls(data_no_missing, phenotype)
    group_cases_controls_df = combine(groupby(data_no_missing, phenotype, skipmissing=true), nrow)
    group_cases_controls_dict = Dict(
        string(val) => n for (val, n) in 
            zip(group_cases_controls_df[!, phenotype], group_cases_controls_df.nrow)
    )
    ncases = get(group_cases_controls_dict, "1", 0)
    ncontrols = get(group_cases_controls_dict, "0", 0)
    return ncases, ncontrols
end

function write_covariates_and_phenotypes_group(data, covariates_list, genotypes_prefix; 
    group_id="all", 
    phenotypes=["SEVERE_COVID_19"], 
    output_prefix="gwas", 
    min_cases_controls=100,
    filters_string=nothing,
    king_cutoff=0.0884
    )
    data = apply_filters(data, filters_string)
    n_phenotypes_passed = 0
    for phenotype in phenotypes
        # Missingness filtering for the phenotype and covariates
        data_no_missing = dropmissing(data, [phenotype, covariates_list...])

        # Relatedness filtering
        phenotype_prefix = string(output_prefix, ".individuals.", group_id, ".", phenotype)
        temp_pheno_file = string(phenotype_prefix, ".temp.txt")
        pheno_file = string(phenotype_prefix, ".txt")
        CSV.write(
            temp_pheno_file, 
            DataFrames.select(data_no_missing, ["FID", "IID"]), 
            header=false, 
            delim="\t"
        )
        # Perform relatedness filtering with KING
        plink2_io = IOBuffer()
        try 
            run(pipeline(
                `plink2 --bfile $genotypes_prefix --keep $temp_pheno_file --king-cutoff $king_cutoff`, 
                stderr=plink2_io
            ))
        catch
            msg = String(take!(plink2_io))
            println(msg)
            if endswith(msg, "No samples remaining after main filters.\n") || endswith(msg, "Error: --king-cutoff requires at least 2 samples.\n")
                @info "Skipping group $group_id and phenotype $phenotype because no samples remain after relatedness filtering with KING cutoff of $king_cutoff."
                continue
            else
                rethrow()
            end
        else
            run(pipeline(`tail -n +2 plink2.king.cutoff.in.id`, stdout=pheno_file))
            rm("plink2.king.cutoff.in.id")
            rm("plink2.king.cutoff.out.id")
        finally
            close(plink2_io)
            rm(temp_pheno_file)
        end

        # Check counts if binary phenotype
        if is_binary_column(data_no_missing, phenotype)
            unrelated_individuals = CSV.read(pheno_file, DataFrame; header=[:FID, :IID])
            data_unrelated_no_missing = innerjoin(data_no_missing, unrelated_individuals, on = [:FID, :IID])
            ncases, ncontrols = n_cases_controls(data_unrelated_no_missing, phenotype)
            if ncontrols < min_cases_controls || ncases < min_cases_controls
                @info "Skipping phenotype $phenotype for group $group_id because it has fewer than $min_cases_controls cases/controls: (cases: $(ncases), controls: $(ncontrols))."
                rm(pheno_file)
                continue
            end
        end

        n_phenotypes_passed += 1
    end

    return n_phenotypes_passed
end

function read_and_process_covariates(covariates_file;
    covariates_string=nothing
    )
    # Read the covariates file
    covariates = CSV.read(covariates_file, DataFrame; missingstring=["", "NA", "NULL", "NAN"])
    # Add user defined covariates
    required_covariate_variables = add_user_defined_covariates!(covariates, covariates_string)

    return covariates, required_covariate_variables
end

function make_groups_and_covariates(
    genotypes_prefix,
    covariates_file; 
    groupby_string=nothing,
    covariates_string="AGE",
    phenotypes_string="SEVERE_COVID_19",
    filters_string=nothing,
    output_prefix="gwas", 
    min_cases_controls=100,
    king_cutoff=0.0884
    )
    phenotypes = split(phenotypes_string, ",")
    # Define additional covariates
    covariates, required_covariate_variables = read_and_process_covariates(covariates_file; covariates_string=covariates_string)
    # Write new covariates to file
    CSV.write(
        string(output_prefix, ".covariates.csv"), 
        covariates,
        missingstring="NA",
        delim="\t"
    )
    # Write required covariates list to file for REGENIE
    open(string(output_prefix, ".covariates_list.txt"), "w") do io
        for covariate in required_covariate_variables
            println(io, covariate)
        end
    end
    # Make groups
    n_groups_passed = 0
    if groupby_string !== nothing
        groupby_variables = split(groupby_string, ",")
        for (groupkey, group) in pairs(groupby(covariates, groupby_variables, skipmissing=true, sort=true))
            group_id = join(groupkey, "_")
            n_phenotypes_passed = write_covariates_and_phenotypes_group(group, required_covariate_variables, genotypes_prefix;
                group_id=group_id,
                phenotypes=phenotypes,
                output_prefix=output_prefix,
                min_cases_controls=min_cases_controls,
                filters_string=filters_string,
                king_cutoff=king_cutoff
            )
            n_groups_passed += n_phenotypes_passed
        end
    else
        n_groups_passed = write_covariates_and_phenotypes_group(covariates, required_covariate_variables, genotypes_prefix; 
                group_id="all",
                phenotypes=phenotypes, 
                output_prefix=output_prefix, 
                min_cases_controls=min_cases_controls,
                filters_string=filters_string,
                king_cutoff=king_cutoff
        )
    end

    n_groups_passed > 0 || throw(ArgumentError("No group passed the min cases/controls threshold."))

    return 0
end