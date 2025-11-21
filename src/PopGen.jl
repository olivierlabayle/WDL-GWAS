module PopGen

using ArgParse
using CSV
using DataFrames
using Downloads
using HTTP
using CairoMakie
using GeneticsMakie
using Colors
using Tables
using RCall
using Statistics
using JSON
using DelimitedFiles
using CategoricalArrays
using MLJBase
using MLJTransforms

const ENSEMBL_SERVER = "https://rest.ensembl.org"

include("cli.jl")
include("prepare_groups.jl")
include("plotting.jl")
include("meta_analysis.jl")
include("fine_mapping.jl")
include("utils.jl")
include("harmonize_gwas_results.jl")

export julia_main

end