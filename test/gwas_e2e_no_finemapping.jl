module TestGWASE2E3

# This end to end test runs with the following conditions:
# - GWAS sofware: plink2
# - a discrete covariate given by the SUPERPOPULATION column
# - a binary phenotype: SEVERE_COVID_19
# - finemapping disabled
# - No loco PCA

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
    "--inputs", joinpath(TESTDIR, "assets", "config", "gwas.no_finemapping.json"),
    "--options", joinpath(TESTDIR, "assets", "config", "gwas.no_finemapping.options.json")
])

# Run the workflow from the package directory
cd(PKGDIR) do
    run(cmd)
end

results_dirs = readdir(joinpath(PKGDIR, "gwas_no_finemapping_outputs/gwas/"), join=true)
results_dir = results_dirs[argmax(mtime(d) for d in results_dirs)]

for shard in (0, 1)
    execution_dir = joinpath(results_dir, "call-make_group_gwas_outputs", "shard-$shard", "execution")
    outputs = readdir(execution_dir)
    # Check summary stats
    summary_stats_file = only(filter(endswith("tsv"), outputs))
    summary_stats = CSV.read(joinpath(execution_dir, summary_stats_file), DataFrame)
    @test nrow(summary_stats) > 10 # not empty
    @test "ERRCODE" in names(summary_stats) # plink2 specific column
    # Check plots
    plots_dir = only(filter(startswith("glob"), outputs))
    plots = readdir(joinpath(execution_dir, plots_dir))
    @test length(plots) == 2
end

end