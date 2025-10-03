module TestPrepareGroups

using Test
using PopGen
using DataFrames
using CSV

PKGDIR = pkgdir(PopGen)
TESTDIR = joinpath(PKGDIR, "test")

@testset "Test is_binary_column" begin
    df = DataFrame(
        A = [0, 1, 0, 1, missing], 
        B = [0, 1, 2, 1, 0], 
        C = ["a", "b", "a", "b", "a"]
    )
    @test PopGen.is_binary_column(df, :A) == true
    @test PopGen.is_binary_column(df, :B) == false
    @test PopGen.is_binary_column(df, :C) == false
end

@testset "Test n_cases_controls" begin
    df = DataFrame(
        PHENO = [0, 1, 1, 0, 1, missing, 0, 0, 1, 1],
        AGE = [25, 35, 45, 55, 65, 75, 85, missing, 50, 60]
    )
    ncases, ncontrols = PopGen.n_cases_controls(df, :PHENO)
    @test ncases == 5
    @test ncontrols == 4
end

@testset "Test apply_filter" begin
    data = DataFrame(
        COHORT = ["FUTURE_HEALTH", "FUTURE_HEALTH", "UKB", "UKB", "UKB", "UKB", "UKB", "UKB", "UKB", "UKB"],
        PRIMARY_DIAGNOSIS = ["COVID-19", missing, "COVID-19", "COVID-19", "PNEUMONIA", "PNEUMONIA", "PNEUMONIA", "PNEUMONIA", "PNEUMONIA", "PNEUMONIA"],
        AGE = [25, 35, 45, 55, 65, 75, 85, missing, 50, 60]
    )
    @test PopGen.apply_filters(data, nothing) === data
    filter_ukb = PopGen.apply_filters(data, "COHORT=UKB")
    @test nrow(filter_ukb) == 8
    @test all(filter_ukb.COHORT .== "UKB")
    filter_covid_ukb = PopGen.apply_filters(data, "PRIMARY_DIAGNOSIS=COVID-19,COHORT=UKB")
    @test nrow(filter_covid_ukb) == 2
    @test all(filter_covid_ukb.COHORT .== "UKB")
    @test all(filter_covid_ukb.PRIMARY_DIAGNOSIS .== "COVID-19")
    @test filter_covid_ukb.AGE == [45, 55]
    filter_age_ukb = PopGen.apply_filters(data, "AGE>=50,AGE<=75,COHORT=UKB")
    @test nrow(filter_age_ukb) == 5
    @test all(filter_age_ukb.COHORT .== "UKB")
    @test filter_age_ukb.AGE == [55, 65, 75, 50, 60]
end

@testset "Test make-groups-and-covariates" begin
    tmpdir = mktempdir()
    output_prefix = joinpath(tmpdir, "gwas")
    covariates_file = joinpath(TESTDIR, "assets", "covariates", "covariates.csv")
    min_cases_controls = 200
    copy!(ARGS, [
        "make-groups-and-covariates", 
        covariates_file,
        "--groupby=SUPERPOPULATION,SEX",
        "--phenotypes=SEVERE_COVID_19",
        "--covariates=AGE,AGE_x_AGE,AGE_x_SEX,COHORT",
        "--output-prefix", output_prefix, 
        "--min-cases-controls", string(min_cases_controls)
    ])
    julia_main()

    updated_covariates = CSV.read(joinpath(tmpdir, "gwas.covariates.csv"), DataFrame)
    # Check covariate file
    expected_covariate_cols = [
        "FID", 
        "IID", 
        "AGE", 
        "SEX",
        "SUPERPOPULATION",
        "AFR",
        "SAS",
        "EAS",
        "AMR",
        "EUR",
        "COHORT",
        "SEVERE_COVID_19",
        "SEVERE_PNEUMONIA",
        "AGE_x_AGE",
        "AGE_x_SEX",
        "COHORT__FUTURE_HEALTH",
        "COHORT__UKB"
    ]
    @test names(updated_covariates) == expected_covariate_cols
    for row in eachrow(updated_covariates)
        @test row.AGE_x_AGE == row.AGE * row.AGE
        if row.SEX === missing
            @test row.AGE_x_SEX === missing
        else
            @test row.AGE_x_SEX == row.AGE * row.SEX
        end
    end
    # Check covariates list
    @test readlines(joinpath(tmpdir, "gwas.covariates_list.txt"),) == ["AGE", "AGE_x_AGE", "AGE_x_SEX", "COHORT__FUTURE_HEALTH", "COHORT__UKB"]
        
    # Check groups files
    case_control_counts = sort(combine(
        groupby(updated_covariates, [:SUPERPOPULATION, :SEX, :SEVERE_COVID_19], skipmissing=true), 
        nrow), 
        :nrow
    )
    groups_failing_min_case_control_df = filter(x -> x.nrow < min_cases_controls, case_control_counts)[!, [:SUPERPOPULATION, :SEX]]
    groups_failing_min_case_control = Set(collect(zip(groups_failing_min_case_control_df.SUPERPOPULATION, groups_failing_min_case_control_df.SEX)))
    @test groups_failing_min_case_control == Set([
        ("ADMIXED", 1),
        ("ADMIXED", 0),
        ("AMR", 0),
        ("EUR", 0),
        ("AFR", 0),
        ("SAS", 0)
    ])
    for ancestry in ["AFR", "AMR", "EAS", "EUR", "SAS"]
        for sex in [0, 1]
            if (ancestry, sex) ∉ groups_failing_min_case_control
                group_key = string(ancestry, "_", sex)
                individuals = CSV.read(joinpath(tmpdir, "gwas.individuals.$group_key.SEVERE_COVID_19.txt"), DataFrame; header=["FID", "IID"])
                joined = innerjoin(updated_covariates, individuals, on = [:FID, :IID])
                @test all(==(ancestry), joined.SUPERPOPULATION)
                @test all(==(sex), joined.SEX)
                @test nrow(dropmissing(joined[!, ["SEVERE_COVID_19", "AGE", "AGE_x_AGE", "AGE_x_SEX", "COHORT__FUTURE_HEALTH", "COHORT__UKB"]])) == nrow(joined)
            end
        end
    end
end

@testset "Test make-groups-and-covariates: no groups with filter" begin
    tmpdir = mktempdir()
    output_prefix = joinpath(tmpdir, "gwas_all")
    covariates_file = joinpath(TESTDIR, "assets", "covariates", "covariates.csv")
    min_cases_controls = 2500
    copy!(ARGS, [
        "make-groups-and-covariates", 
        covariates_file,
        "--output-prefix", output_prefix,
        "--phenotypes=SEVERE_COVID_19,SEVERE_PNEUMONIA",
        "--filters=AGE>=50,AGE<=75",
        "--covariates=AGE",
        "--min-cases-controls", string(min_cases_controls)
    ])
    julia_main()

    # Check covariate file and group 
    covariates = CSV.read(joinpath(tmpdir, "gwas_all.covariates.csv"), DataFrame)
    # SEVERE_COVID_19 is dropped because it has fewer than 2500 cases/controls
    @test !isfile(joinpath(tmpdir, "gwas_all.individuals.all.SEVERE_COVID_19.txt"))
    # The group consists in all individuals
    individuals = sort(CSV.read(
        joinpath(tmpdir, "gwas_all.individuals.all.SEVERE_PNEUMONIA.txt"), 
        DataFrame; 
        header=["FID", "IID"])
    )
    expected_individuals = sort(
        dropmissing(filter(x -> x.AGE >= 50 && x.AGE <= 75, covariates), ["SEVERE_PNEUMONIA", "AGE"]
        )[!, ["FID", "IID"]])
    @test individuals == expected_individuals

    # Check covariates list
    @test readlines(joinpath(tmpdir, "gwas_all.covariates_list.txt"),) == ["AGE"]
end


end

true