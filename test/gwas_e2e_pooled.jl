module TestGWASE2E2

# This end to end test runs with the following conditions:
# - GWAS sofware: saige
# - a pooled analysis
# - A filter on: SUPERPOPULATION=EUR and AGE>50
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
gwas_execution_dir = joinpath(results_dir, "call-make_group_gwas_outputs", "shard-0", "execution")
gwas_outputs = readdir(gwas_execution_dir, join=true)
gwas_results = CSV.read(
    joinpath(gwas_execution_dir, "all.AGE.gwas.tsv"),
    DataFrame
)
@test Set(gwas_results.CHROM) == Set([1, 2, 3])
@test nrow(gwas_results) > 5
@test maximum(gwas_results.N) < 400 # Filter applied
## SAIGE colnames
@test names(gwas_results) == ["CHROM", "POS", "ID", "ALLELE_0", "ALLELE_1", "ALLELE_1_FREQ", "BETA", "SE", "LOG10P", "N", "ALLELE_1_COUNT", "MISSING_RATE", "T_STAT", "VAR"]
## Plots
plots_dir_index = findfirst(isdir, gwas_outputs)
@test readdir(gwas_outputs[plots_dir_index]) == [
    "all.AGE.gwas.manhattan.png",
    "all.AGE.gwas.qq.png"
]
# No Finemapping results (only headers)
@test countlines(joinpath(results_dir, "call-make_group_finemapping_outputs", "shard-0", "execution", "all.AGE.finemapping.tsv")) == 1

end

true