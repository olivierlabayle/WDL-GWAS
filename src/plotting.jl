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
    fig = Figure(size = (1200, 600))
    ax = Axis(fig[1, 1], xlabel="Chromosome", ylabel="-log10(P)")
    GeneticsMakie.plotgwas!(ax, gwas_results, build=build)
    hidespines!(ax, :t, :r)
    Label(fig[1, 1, Top()], text = title, fontsize = 20)
    resize_to_layout!(fig)
    return fig
end

function get_gencode_annotations(;gencode_file = "gencode.v49.annotation.gtf.gz")
    if !isfile(gencode_file)
        Downloads.download("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz", gencode_file)
    end
    gencode = CSV.read(gencode_file, DataFrame; 
        delim = "\t", 
        comment = "#", 
        header = ["seqnames", "source", "feature", "start", "end", "score", "strand", "phase", "info"]
    )
    GeneticsMakie.parsegtf!(gencode)
    return select!(gencode, :seqnames, :feature, :start, :end, :strand, :gene_id, :gene_name, :gene_type, :transcript_id)
end

function qqplot(gwas_results; title="QQ Plot")
    fig = Figure(size = (800, 800))
    ax = Axis(fig[1, 1])
    GeneticsMakie.plotqq!(ax, gwas_results; ystep = 5)
    hidespines!(ax, :t, :r)
    Label(fig[1, 1, Top()], text = title, fontsize = 20)
    resize_to_layout!(fig)
    return fig
end

function get_genomic_features(region; features=["gene", "transcript", "cds", "exon", "regulatory", "motif"])
    ext = string("/overlap/region/human/", region, "?", join(map(f -> "feature=$(f)", features), ";"))
    headers=Dict("Content-Type" => "application/json", "Accept" => "application/json")
    r = HTTP.get(ENSEMBL_SERVER*ext, headers)
    return JSON.parse(String(r.body))
end

function region_plot(region_data)
    chr = only(unique((region_data[!, "CHR"])))
    region = string(chr, ":", minimum(region_data.BP), "-", maximum(region_data.BP))
    genomic_features = get_genomic_features(region; features=["gene"])
    credible_sets = sort(unique(skipmissing(region_data.CS)))
    lead_snps = split(only(unique(region_data.LOCUS_ID)), "_&_")
    markersize = 10
    sig_markersize = 14
    threshold = 5e-8
    fig = Figure(size = (1000, 800))
    # GWAS P-values
    ax1 = Axis(fig[1, 1]; ylabel="-log10(P)", 
        xticksvisible=false, 
        xticklabelsvisible=false, 
        xgridvisible=false,
        ygridvisible=false
    )
    hidespines!(ax1, :t, :r)
    ## Plot region data points
    log10ps = -log10.(region_data.P)
    scatter!(ax1, 
        collect(region_data.BP),
        log10ps,
        markersize=markersize,
        color=region_data.PHASED_R2,
        colormap = :heat
    )
    ## Add significance threshold
    hlines!(ax1, -log10(threshold), color=:green)
    ## Add Lead SNPs
    lead_snps_info = region_data[region_data.ID .== lead_snps, :]
    vlines!(ax1, lead_snps_info.BP, color=:red, linestyle=:dash)
    text!(ax1, lead_snps_info.BP, -log10.(lead_snps_info.P), text=lead_snps, color=:black, align = (:left, :bottom), fontsize=12)
    # Fine Mapping PIPs
    credible_sets_colors = distinguishable_colors(
        length(credible_sets), 
        [Colors.RGB(1,1,1), Colors.RGB(0,0,0)], 
        dropseed=true
    )
    cs_to_color = Dict{Any, Any}(credible_sets .=> credible_sets_colors)
    ax2 = Axis(fig[2, 1]; ylabel="PIP",
        xticksvisible=false, 
        xticklabelsvisible=false, 
        xgridvisible=false,
        ygridvisible=false
    )
    hidespines!(ax2, :t, :r)
    scatter!(ax2,
        collect(region_data.BP), 
        collect(region_data.PIP); 
        color=(:grey, 0.5), 
        markersize=markersize, 
        colormap=:tab10
    )
    for (cs, cs_color) in cs_to_color
        cs_variants = findall(x -> (x !== missing) && (x == cs), region_data.CS)
        pips = region_data.PIP[cs_variants]
        bps = region_data.BP[cs_variants]
        scatter!(ax2, 
            bps, 
            pips; 
            color=cs_color, 
            markersize=sig_markersize,
            marker=:star5,
            colormap=:tab10
        )
    end
    pips_legend_elements = [(string("CS ", cs), PolyElement(color=cs_color, colorrange=1:10, colormap=:tab10)) for (cs, cs_color) in enumerate(credible_sets_colors)]
    Legend(fig[2, 2], getindex.(pips_legend_elements, 2), getindex.(pips_legend_elements, 1), "Fine Mapping CSs", framevisible = false)
    # Genomic annotations
    ax3 = Axis(fig[3, 1]; 
        xlabel=string("Chr", chr),
        ylabel="Genes",
        yticksvisible=false,
        yticklabelsvisible=false,
        xgridvisible=false,
        ygridvisible=false
    )
    hidespines!(ax3, :t, :r)
    for (y_coord, feature) in enumerate(genomic_features)
        lines!(ax3, 
            [feature["start"], feature["end"]], 
            [y_coord, y_coord], 
            color=:blue,
            linewidth=6
        )
        name = haskey(feature, "external_name") ? feature["external_name"] : feature["gene_id"]
        x_coord = (feature["start"] + feature["end"])/2
        text!(ax3, x_coord, y_coord + 0.2, text=name, align = (:center, :bottom), color=:black, fontsize=12)
    end
    return fig
end

function make_plots(gwas_file, finemapping_file; maf=0.01, output_prefix = "gwas.plot")
    group, phenotype, _ = split(basename(gwas_file), ".")
    gwas_results = harmonize_gwas_results(CSV.read(gwas_file, DataFrame, delim="\t"))
    maf_filtered_gwas_results = filter(
        x -> x.A1FREQ > maf, 
        gwas_results
    )
    # Plot Manhattan
    title = string(phenotype, "\n", group)
    fig = manhattan_plot(maf_filtered_gwas_results; title=title, build=38)
    save(string(output_prefix, ".manhattan.png"), fig)
    # Plot QQ
    fig = qqplot(maf_filtered_gwas_results; title=title)
    save(string(output_prefix, ".qq.png"), fig)
    # Plot locuszoom for top hits
    finemapping_results = harmonize_finemapping_results(CSV.read(finemapping_file, DataFrame, delim="\t"))
    for (locus_key, locus_group) in pairs(groupby(finemapping_results, :LOCUS_ID))
        region_data = innerjoin(
            gwas_results,
            DataFrames.select(locus_group, [:ID, :REF, :ALT, :PIP, :CS, :LOCUS_ID, :PHASED_R2]), 
            on=[:ID]
        )
        fig = region_plot(region_data)
        save(string(output_prefix, ".", locus_key.LOCUS_ID, ".locuszoom.png"), fig)
    end
    return 0
end