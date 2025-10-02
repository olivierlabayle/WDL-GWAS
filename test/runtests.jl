using PopGen
using Test

PKGDIR = pkgdir(PopGen)
TESTDIR = joinpath(PKGDIR, "test")

@testset "WDL-GWAS.jl" begin
    # Unit Tests
    @test include(joinpath(TESTDIR, "meta_analysis.jl"))
    @test include(joinpath(TESTDIR, "plotting.jl"))
    @test include(joinpath(TESTDIR, "fine_mapping.jl"))
    @test include(joinpath(TESTDIR, "utils.jl"))
    @test include(joinpath(TESTDIR, "prepare_groups.jl"))

    # End to end Tests of the WDL workflow
    @test include(joinpath(TESTDIR, "gwas_e2e_groupby.jl"))
end
