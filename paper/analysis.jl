using Pkg; Pkg.activate("paper")
using CSV
using DataFrames
using CairoMakie
using GeneticsMakie

const data_dir = expanduser("~/Data/WDL_GWAS")

const GENE_ATLAS_DIR = joinpath(data_dir, "gene_atlas")

function parse_pvalue(log10_pval::AbstractString)
    if log10_pval == "NA"
        return NaN
    else
        return parse_pvalue(parse(Float64, log10_pval))
    end
end

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

function liftover_gwas(gwas_file, chain_file, output)
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
    run(`CrossMap bed --chromid s --unmap-file $unmapped_file $chain_file $input_bed_file $output_bed_file`)
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

function liftover_gwases(data_dir)
    chain_file = joinpath(data_dir, "hg19ToHg38.over.chain.gz")
    @info "Lifting over: icd10-C18-both_sexes.tsv"
    liftover_gwas(
        joinpath(data_dir, "icd10-C18-both_sexes.tsv"), 
        chain_file, 
        joinpath(data_dir, "icd10-C18-both_sexes.GRCh38.tsv")
    )
    @info "Lifting over: continuous-23104-both_sexes-irnt.tsv"
    liftover_gwas(
        joinpath(data_dir, "continuous-23104-both_sexes-irnt.tsv"), 
        chain_file, 
        joinpath(data_dir, "continuous-23104-both_sexes-irnt.GRCh38.tsv")
    )
end

function ref_alt_matching_status(row)
    if row.ALLELE0 == row.PANUKB_REF && row.ALLELE1 == row.PANUKB_ALT
        return "match"
    elseif row.ALLELE0 == row.PANUKB_ALT && row.ALLELE1 == row.PANUKB_REF
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

function liftover_gene_atlas(data_dir)
    variants_info_files = [joinpath(GENE_ATLAS_DIR, "snps.imputed.chr$chr.csv.gz") for chr in 1:22]
    trait = "23104-0.0"
    bmi_files = [joinpath(GENE_ATLAS_DIR, "imputed.allWhites.$trait.chr$chr.csv.gz") for chr in 1:22]
    colorectal_cancer_files = [joinpath(GENE_ATLAS_DIR, "imputed.allWhites.cancer_c_C18.chr$chr.csv.gz") for chr in 1:22]

    variants_info = load_csv_files(variants_info_files)
    gwas_results = load_csv_files(bmi_files)
    gwas_results = innerjoin(gwas_results, variants_info, on=:SNP)
    select(
        gwas_results,
        "CHROM",
        "Position" => "START",
        "Position" => "END",
        "A1" => "GA_REF",
        "A2" => "GA_ALT",
        "NBETA-$trait" => "GA_BETA",
        "NSE-$trait" => "GA_SE",
        "PV-$trait" => "GA_PV"
    )
    all(gwas_results.A2 .== gwas_results.ALLELE)

end

function compare_pvalues(wdl_gwas_file, panukb_gwas_file)
    trait = "BMI"
    wdl_gwas_file = joinpath(data_dir, string("all.", trait, ".gwas.tsv"))
    panukb_gwas_file = joinpath(data_dir, "continuous-23104-both_sexes-irnt.GRCh38.tsv")
    # panukb_gwas_file = joinpath(data_dir, "icd10-C18-both_sexes.GRCh38.tsv")

    wdl_gwas = CSV.read(wdl_gwas_file, DataFrame)
    wdl_gwas.CHROM = string.(wdl_gwas.CHROM)
    panukb_gwas = CSV.read(panukb_gwas_file, DataFrame)
    merged_gwas = innerjoin(wdl_gwas, panukb_gwas, on=[:CHROM, :GENPOS]) # TODO: only 17 M remain here, why?
    transform!(merged_gwas, 
        AsTable([:ALLELE0, :ALLELE1, :PANUKB_REF, :PANUKB_ALT]) => ByRow(ref_alt_matching_status) => :REF_ALT_STATUS
    )
    subset!(merged_gwas, :REF_ALT_STATUS => x -> x .!== "mismatch")
    merged_gwas[merged_gwas.REF_ALT_STATUS .== "swapped", :]

    
    plot_data = subset(
        dropmissing(merged_gwas[!, [:PANUKB_LOG10P, :LOG10P, :CHROM, :ALLELE0, :ALLELE1, :PANUKB_LOW_CONFIDENCE, :A1FREQ]]),
        :PANUKB_LOW_CONFIDENCE => x -> x .== false,
        :A1FREQ => x -> x .> 0.01,
        :CHROM => x -> x .== "1"
    )
    fig = Figure()
    ax = Axis(fig[1, 1], 
        xlabel="PanUKB GWAS", 
        ylabel="WDL GWAS",
        # xscale=CairoMakie.Makie.pseudolog10,
        # yscale=CairoMakie.Makie.pseudolog10
    )
    scatter!(ax, 
        plot_data.PANUKB_LOG10P, 
        plot_data.LOG10P;
        markersize=6
    )
    ablines!(ax, 0, 1, color=:red)
    fig

end

for trait in ("BMI", "COLORECTAL_CANCER")
    @info "Plotting Trait: " trait
    gwas_results = harmonize_gwas_results(CSV.read(joinpath(data_dir, string("all.", trait, ".gwas.tsv")), DataFrame))
    title = replace(trait, "_" => " ")
    fig = manhattan_plot(gwas_results; title="")
    save(joinpath(data_dir, string(trait, ".manhattan.png")), fig)
    fig = qqplot(gwas_results; title="")
    save(joinpath(data_dir, string(trait, ".qq.png")), fig)
end

# phenotypes_manifest = CSV.read(joinpath(data_dir, "phenotype_manifest.tsv"), DataFrame)

function sandbox()
    s = subset(
        select(phenotypes_manifest, :description, :description_more, :aws_path, :filename, :num_pops, :pops), 
        :description => x -> occursin.("BMI", x)
    )
    s = subset(
        select(phenotypes_manifest, :description, :description_more, :aws_path, :filename, :num_pops, :pops), 
        :description => x -> occursin.("C18 Malignant neoplasm of colon", x)
    )

    colorectal_panukb = CSV.read(joinpath(data_dir, "icd10-C18-both_sexes.tsv"), DataFrame)

end