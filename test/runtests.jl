using Aqua
using Test
using NCDatasets
using EmpiricalOrthogonalFunctions

testdir = @__DIR__
datadir = joinpath(testdir, "data")

ds = NCDataset(joinpath(datadir,"sst.mnmean.1990-present.nc"))

@testset "EmpiricalOrthogonalFunctions" begin

    sst = ds["sst"][:]
    dsin = copy(sst)

    # @testset "EOF Construction" begin
        eof = EmpiricalOrthogonalFunction(dsin)

        dsin = permutedims(sst,(3,1,2))
        eof = EmpiricalOrthogonalFunction(dsin; timedim=1)

        dsin = permutedims(sst,(1,3,2))
        eof = EmpiricalOrthogonalFunction(dsin; timedim=2)
    # end

    # @testset "EOF Analysis" begin

        # @testset "eigenvals" begin
            @test totalanomalyvar(eof) ≈ 2673.0002f0 atol=1e-3
            @test eigenvalues(eof,n=5)[1] ≈ 1332.6342
            @test variancefraction(eof, n=100)[1] ≈ 0.4985537f0 atol=1e-3
        # end

        # @testset "pcs" begin
            @test pcs(eof,n=1)[1,1] ≈ -10.937908 atol=1e-3
            @test pcs(eof,scaling=:divide)[1,1] ≈ -0.299626 atol=1e-3
            @test pcs(eof,scaling=:multiply, n=10)[1,1] ≈ -399.291
        # end

        # @testset "eofs" begin
            @test eofs(eof,n=5)[1,:] ≈ [0.0279191, -0.00658843, 0.02273, -0.00208263, 0.0347027] atol=1e-3
            @test eofs(eof,n=2,scaling=:divide)[1,:] ≈ [0.000764796, -0.000231014] atol=1e-3
            @test eofs(eof,n=3,scaling=:multiply)[1,:] ≈ [ 1.01919, -0.1879, 0.293734] atol=1e-3
        # end

        # @testset "projection" begin
            @test all((reshape(reconstruct(eof), size(sst)) ≈ sst)...)
            @test projectfield(eof, collect(eof.dataset))[1,1] ≈ 230.678 atol=1e-3
        # end

        # @testset "rotation" begin
            reof = orthorotation(eof,n=4)

            @test eofs(reof)[1,:] ≈ [0.02492, -0.00625398, 0.0254053, 0.00619081] atol=1e-3

        # end
    # end
end

close(ds)

Aqua.test_all(EmpiricalOrthogonalFunctions)

# end tests
