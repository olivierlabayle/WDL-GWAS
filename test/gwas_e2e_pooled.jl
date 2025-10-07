module TestGWASE2E2

# This end to end test runs with the following conditions:
# - a pooled analysis defined by the SUPERPOPULATION column
# - a continuous phenotype: AGE
# - finemapping thresholds are not met

using Test
using PopGen
using DataFrames
using CSV

PKGDIR = pkgdir(PopGen)
TESTDIR = joinpath(PKGDIR, "test")

config = if Sys.isapple()
    "-Dconfig.file=config/cromwell.macOS-dev.conf"
else
    "-Dconfig.file=config/cromwell.local.conf"
end

cmd = Cmd([
    "java", config,
    "-jar", ENV["CROMWELL_PATH"],
    "run", joinpath(PKGDIR, "workflows", "gwas.wdl"),
    "--inputs", joinpath(TESTDIR, "assets", "config", "gwas.pooled.json"),
    "--options", joinpath(TESTDIR, "assets", "config", "gwas.pooled.options.json")
])

# Run the workflow from the package directory
cd(PKGDIR) do
    run(cmd)
end

results_dirs = readdir(joinpath(PKGDIR, "gwas_pooled_outputs/gwas/"), join=true)
results_dir = results_dirs[argmax(mtime(d) for d in results_dirs)]

# GWAS results
gwas_results = CSV.read(
    joinpath(results_dir, "call-merge_gwas_group_chr_results", "shard-0", "execution", "all.AGE.gwas.tsv"),
    DataFrame
)
@test Set(gwas_results.CHROM) == Set([1, 2, 3])
@test nrow(gwas_results) > 5

# No Finemapping results (only headers)
@test countlines(joinpath(results_dir, "call-merge_fp_group_chr_results", "shard-0", "execution", "all.AGE.finemapping.tsv")) == 1

# Plots
plots_dir = joinpath(results_dir, "call-gwas_group_plots", "shard-0", "execution")
plots_subdir = readdir(plots_dir, join=true)[1]
@test readdir(plots_subdir) == [
    "all.AGE.manhattan.png",
    "all.AGE.qq.png"
]

end

true