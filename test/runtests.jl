using Test
using Random

#=

@testset "step" begin

   n = 2
   ϵ = 0.1
   q = [[0.5,0.2], [0.5,0.8]]
   v = [[0.0,1.0], [0.0,-1.0]]
   collisions = PCollisionMatrix(n, q, v, ϵ)

   @test step!(n, ϵ, q, v, collisions) ≈ 0.2

   @show q
   @test all(vcat(q...) .≈ [0.5, 0.4,  0.5, 0.6])
   @show v
   @test all(vcat(v...) .≈ [1.0, 0.0, -1.0, 0.0])

   @test step!(n, ϵ, q, v, collisions) ≈ 0.5


end
=#

@testset "overlap" begin

end


@testset "init_particles" begin
 
     import Uchiyama: init_particles
     n = 100
     ϵ = 1 / n
     c = trunc(Int, 1000 * 2ϵ)
     rng = MersenneTwister(1234)
     
     q, v = init_particles(n, ϵ, rng)

end

include("test_tempscoll.jl")
