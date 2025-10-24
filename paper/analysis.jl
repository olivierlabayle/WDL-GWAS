using Pkg; Pkg.activate("paper")
using CSV
using DataFrames
using CairoMakie
using GeneticsMakie
using HTTP
using JSON
using Colors

include(joinpath("src", "plotting.jl"))

const ENSEMBL_SERVER = "https://rest.ensembl.org"
const DATA_DIR = expanduser("~/Data/WDL_GWAS")
const WDL_GWAS_DIR = joinpath(DATA_DIR, "wdl_gwas")
const GENE_ATLAS_DIR = joinpath(DATA_DIR, "gene_atlas")
const PAN_UKB_DIR = joinpath(DATA_DIR, "panukb")
const PLOTS_DIR = joinpath(DATA_DIR, "plots")
const CHAIN_FILE = joinpath(DATA_DIR, "hg19ToHg38.over.chain.gz")

function parse_pvalue(log10_pval::AbstractString)
    if log10_pval == "NA"
        return NaN
    else
        return parse_pvalue(parse(Float64, log10_pval))
    end
end

wdl_gwas_filename(trait) = joinpath(WDL_GWAS_DIR, string("all.", trait, ".gwas.tsv"))

parse_pvalue(log10_pval::Real) = exp10(-log10_pval)

parse_a1freq(freq::AbstractString) = freq == "NA" ? NaN : parse(Float64, freq)

parse_a1freq(freq::Real) = freq

function harmonize_gwas_results(gwas_results)
    return DataFrames.transform(gwas_results, 
        :CHROM => (x -> string.(x)) => :CHR,
        :GENPOS => :BP,
        :ID => :SNP,
        :LOG10P => (x -> parse_pvalue.(x))  => :P,
        :A1FREQ => (x -> parse_a1freq.(x)) => :A1FREQ
    )    
end

function harmonize_finemapping_results(finemapping_results)
    return DataFrames.transform(finemapping_results, 
        :CHROM => (x -> string.(x)) => :CHR,
        :POS => :BP,
        :ID => :SNP
    )    
end

function manhattan_plot(gwas_results; title="Manhattan Plot", build=38)
    fig = Figure(size = (792, 408))
    ax = Axis(fig[1, 1], xlabel="Chromosome", ylabel="-log10(P)", yscale=CairoMakie.Makie.pseudolog10)
    GeneticsMakie.plotgwas!(ax, gwas_results, build=build)
    hidespines!(ax, :t, :r)
    # Label(fig[1, 1, Top()], fontsize = 20)
    resize_to_layout!(fig)
    return fig
end

function qqplot(gwas_results; title="QQ Plot")
    fig = Figure(size = (792, 792))
    ax = Axis(fig[1, 1])
    GeneticsMakie.plotqq!(ax, gwas_results; ystep = 5)
    hidespines!(ax, :t, :r)
    # Label(fig[1, 1, Top()], text = title, fontsize = 20)
    resize_to_layout!(fig)
    return fig
end

function liftover_gwas(gwas_file, output; )
    tmpdir = mktempdir()
    input_bed_file = joinpath(tmpdir, "GRCh37_positions.bed")
    output_bed_file = joinpath(tmpdir, "GRCh38_positions.bed")
    unmapped_file = joinpath(tmpdir, "unmapped.bed")
    # Load GWAS Results and convert to BED format
    gwas_results = CSV.read(gwas_file, DataFrame)
    relevant_columns = [
        "ref",
        "alt",
        "beta_EUR",
        "se_EUR",
        "neglog10_pval_EUR",
        "low_confidence_EUR"
    ]
    CSV.write(input_bed_file,
        select(gwas_results, 
            "chr", 
            "pos" => "pos_start", 
            "pos" => "pos_end",
            relevant_columns...
        ),
        header=false,
        delim="\t"
    )
    # Liftover
    run(`CrossMap bed --chromid s --unmap-file $unmapped_file $CHAIN_FILE $input_bed_file $output_bed_file`)
    # Reload, format and write output
    gwas_results = CSV.read(
        output_bed_file, 
        DataFrame; 
        header=["chr", "pos_start", "pos_end", relevant_columns...],
        missingstring="NA"
    )
    @assert all(gwas_results.pos_start .== gwas_results.pos_end)
    relevant_gwas_results = select(gwas_results, 
            :chr => :CHROM, 
            :pos_start => :GENPOS, 
            :ref => :PANUKB_REF, 
            :alt => :PANUKB_ALT,
            :beta_EUR => :PANUKB_BETA, 
            :se_EUR => :PANUKB_SE, 
            :neglog10_pval_EUR => :PANUKB_LOG10P,
            :low_confidence_EUR => :PANUKB_LOW_CONFIDENCE
    )
    CSV.write(output,
        relevant_gwas_results, 
        delim="\t"
    )
    rm(tmpdir, recursive=true)
end

"""

gwas_results is assumed to have columns: CHROM, POS, OTHER_COLUMNS...
"""
function liftover_gwas_new(gwas_results, output)
    tmpdir = mktempdir()
    input_bed_file = joinpath(tmpdir, "GRCh37_positions.bed")
    output_bed_file = joinpath(tmpdir, "GRCh38_positions.bed")
    unmapped_file = joinpath(tmpdir, "unmapped.bed")
    # Write GWAS results to BED format
    column_names = names(gwas_results)
    CSV.write(input_bed_file,
        insertcols(gwas_results, 3, :END_POS => gwas_results.POS),
        header=false,
        missingstring="NA",
        delim="\t"
    )
    # Liftover
    run(`CrossMap bed --chromid s --unmap-file $unmapped_file $CHAIN_FILE $input_bed_file $output_bed_file`)
    # Reload, format and write output
    gwas_results = CSV.read(
        output_bed_file, 
        DataFrame; 
        header=["CHROM", "POS", "POS_END", column_names[3:end]...],
        missingstring="NA",
        delim="\t"
    )
    @assert all(gwas_results.POS .== gwas_results.POS_END)
    CSV.write(output,
        select(gwas_results, column_names), 
        delim="\t"
    )
    rm(tmpdir, recursive=true)
end

liftover_panukb_filename(trait) = joinpath(PAN_UKB_DIR, "gwas.$trait.GRCh38.tsv")

function liftover_panukb()
    @info "Lifting over PANUKB"
    trait_map = [
        "icd10-C18-both_sexes.tsv" => "COLORECTAL_CANCER",
        "continuous-23104-both_sexes-irnt.tsv" => "BMI"
    ]
    for (panukb_file, trait) in trait_map
        isfile(liftover_panukb_filename(trait)) && continue
        @info "Lifting over PANUKB trait: $trait"
        gwas_results = CSV.read(joinpath(PAN_UKB_DIR, panukb_file), DataFrame; missingstring="NA")
        gwas_results = select(gwas_results, 
            :chr => :CHROM, 
            :pos => :POS, 
            :ref => :PANUKB_REF, 
            :alt => :PANUKB_ALT,
            :beta_EUR => :PANUKB_BETA, 
            :se_EUR => :PANUKB_SE, 
            :neglog10_pval_EUR => :PANUKB_LOG10P,
            :low_confidence_EUR => :PANUKB_LOW_CONFIDENCE
        )
        liftover_gwas_new(gwas_results, liftover_panukb_filename(trait))
    end
end

function ref_alt_matching_status(row)
    ref_1, alt_1, ref_2, alt_2 = row
    if ref_1 == ref_2 && alt_1 == alt_2
        return "match"
    elseif ref_1 == alt_2 && alt_1 == ref_2
        return "swapped"
    else
        return "mismatch"
    end
end

function load_csv_files(files)
    return mapreduce(vcat, files) do file
        CSV.read(file, DataFrame)
    end
end

function load_ga_variants_info(files)
    return mapreduce(vcat, files) do file
        df = CSV.read(file, DataFrame)
        df.CHROM .= replace(split(basename(file), ".")[3], "chr" => "")
        df
    end
end

liftover_ga_filename(trait) = joinpath(GENE_ATLAS_DIR, "gwas.$trait.GRCh38.tsv")

function liftover_gene_atlas()
    @info "Lifting over geneATLAS"
    variants_info_files = [joinpath(GENE_ATLAS_DIR, "snps.imputed.chr$chr.csv.gz") for chr in 1:22]
    variants_info = load_ga_variants_info(variants_info_files)
    trait_map = [
        "23104-0.0" => "BMI",
        "cancer_c_C18" => "COLORECTAL_CANCER"
    ]
    for (ga_trait, trait) in trait_map
        @info "Lifting over geneATLAS trait: $trait"
        isfile(liftover_ga_filename(trait)) && continue
        gwas_files = [joinpath(GENE_ATLAS_DIR, "imputed.allWhites.$ga_trait.chr$chr.csv.gz") for chr in 1:22]
        gwas_results = load_csv_files(gwas_files)
        gwas_results = innerjoin(gwas_results, variants_info, on=:SNP)
        gwas_results = select(gwas_results,
            "CHROM",
            "Position" => "POS",
            "SNP" => "GA_ID",
            "A1" => "GA_REF",
            "A2" => "GA_ALT",
            "MAF" => "GA_MAF",
            "NBETA-$ga_trait" => "GA_BETA",
            "NSE-$ga_trait" => "GA_SE",
            "PV-$ga_trait" => (x -> - log10.(x)) => "GA_LOG10P",
            "PV-$ga_trait" => "GA_PV",
        )
        liftover_gwas_new(gwas_results, liftover_ga_filename(trait))
    end
end

find_column(gwas_results; type="REF") = only(filter(x -> endswith(x, type), names(gwas_results)))

function compare_pvalues(gwas_results_1, gwas_results_2)
    ref_1, alt_1, log10p_1 = find_column(gwas_results_1; type="REF"), find_column(gwas_results_1; type="ALT"), find_column(gwas_results_1; type="LOG10P")
    ref_2, alt_2, log10p_2 = find_column(gwas_results_2; type="REF"), find_column(gwas_results_2; type="ALT"), find_column(gwas_results_2; type="LOG10P")

    merged_gwas = innerjoin(gwas_results_1, gwas_results_2, on=[:CHROM, :POS]) # TODO: only 17 M remain here, why?

    transform!(merged_gwas,
        AsTable([ref_1, alt_1, ref_2, alt_2]) => ByRow(ref_alt_matching_status) => :REF_ALT_STATUS
    )
    subset!(merged_gwas, :REF_ALT_STATUS => x -> x .!== "mismatch")
    
    plot_data = dropmissing(merged_gwas[!, [log10p_1, log10p_2]])
    
    fig = Figure()
    ax = Axis(fig[1, 1], 
        xlabel=split(log10p_1, "_")[1], 
        ylabel=split(log10p_2, "_")[1],
    )
    scatter!(ax, 
        plot_data[!, log10p_1], 
        plot_data[!, log10p_2];
        markersize=6
    )
    ablines!(ax, 0, 1, color=:red)
    fig
end

function main()
    liftover_gene_atlas()
    liftover_panukb()

    trait = "BMI"
    # Load WDL GWAS
    wdl_gwas = CSV.read(wdl_gwas_filename(trait), DataFrame, types=Dict("CHROM" => String))
    rename!(wdl_gwas, 
        :GENPOS => :POS, 
        :ALLELE0 => :WDLGWAS_REF, 
        :ALLELE1 => :WDLGWAS_ALT,
        :A1FREQ => :WDLGWAS_ALT_FREQ,
        :BETA => :WDLGWAS_BETA,
        :SE => :WDLGWAS_SE,
        :LOG10P => :WDLGWAS_LOG10P
    )
    wdl_gwas_frequent = wdl_gwas[0.01 .< wdl_gwas.WDLGWAS_ALT_FREQ .< 0.99, :]
    # Load PANUKB GWAS
    panukb_gwas = CSV.read(liftover_panukb_filename(trait), DataFrame)
    panukb_gwas = CSV.read(joinpath(PAN_UKB_DIR, "gwas.$trait.GRCh38.tsv"), DataFrame)

    # Load geneATLAS GWAS
    ga_gwas = CSV.read(liftover_ga_filename(trait), DataFrame)
    subset!(ga_gwas, :GA_PV => x -> x .!= 0 .&& x .!= 1, :GA_MAF => x -> x .> 0.01)

    compare_pvalues(wdl_gwas_frequent, ga_gwas)

    for trait in ("BMI", "COLORECTAL_CANCER")
        @info "Plotting Trait: " trait
        gwas_results = harmonize_gwas_results(CSV.read(joinpath(WDL_GWAS_DIR, string("all.", trait, ".gwas.tsv")), DataFrame))
        fig = manhattan_plot(gwas_results; title=" ")
        save(joinpath(PLOTS_DIR, string(trait, ".manhattan.png")), fig)
        fig = qqplot(gwas_results; title=" ")
        save(joinpath(PLOTS_DIR, string(trait, ".qq.png")), fig)
    end

    # Finemapping
    trait = "BMI"
    gwas_results = harmonize_gwas_results(CSV.read(joinpath(WDL_GWAS_DIR, string("all.", trait, ".gwas.tsv")), DataFrame))
    finemapping_results = harmonize_finemapping_results(CSV.read(joinpath(WDL_GWAS_DIR, string("all.", trait, ".finemapping.tsv")), DataFrame, delim="\t"))
    bmi_fp = finemapping_results[finemapping_results.LOCUS_ID .== "16:53767042:T:C", :]
    region_data = innerjoin(
            gwas_results,
            DataFrames.select(bmi_fp, [:ID, :REF, :ALT, :PIP, :CS, :LOCUS_ID, :UNPHASED_R2, :SUSIE_CONVERGED]), 
            on=[:ID]
    )

    fig = region_plot(region_data)

    rs1421085_locus[rs1421085_locus.ID .== "16:53767042:T:C", :]
end