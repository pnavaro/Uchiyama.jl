using Test
using Random
using LinearAlgebra
using Uchiyama

include("test_periodic_hard_spheres.jl")

@testset "step" begin


   n = 2
   ϵ = 0.1
   q = [[0.5,0.2], [0.5,0.8]]
   v = [[0,1], [0,-1]]

   rng = MersenneTwister(42)
   p = SquareParticles(rng, n, ϵ)
   for i in 1:n
       p.q[i] = q[i]
       p.v[i] = v[i]
   end

   collisions = PeriodicCollisions(p)

   @test step!(p, collisions) ≈ 0.2

   @test p.q[1] ≈ [0.5,0.4] 
   @test p.v[1] ≈ [0,-1]
   @test p.q[2] ≈ [0.5,0.6] 
   @test p.v[2] ≈ [0, 1]

   @test step!(p, collisions) ≈ 0.3

   @test p.q[1] ≈ [0.5,0.1] 
   @test p.v[1] ≈ [0,1]
   @test p.q[2] ≈ [0.5,0.9] 
   @test p.v[2] ≈ [0,-1]

end

@testset "overlap" begin

   import Uchiyama: overlap

   q  = [0.5, 0.5]

   p  = [0.6, 0.6]
   ϵ  = 0.1
   @test overlap( p, q, ϵ ) 
   p  = [0.65, 0.65]
   @test overlap( p, q, ϵ ) 
   p  = [0.65, 0.5]
   @test overlap( p, q, ϵ ) 
   p  = [0.5, 0.65]
   @test overlap( p, q, ϵ ) 
   p  = [0.8, 0.8]
   @test !overlap( p, q, ϵ )
   p  = [0.1, 0.1]
   @test !overlap( p, q, ϵ )
   p  = [0.69, 0.69]
   @test overlap( p, q, ϵ ) 
   p  = [0.41, 0.4]
   @test overlap( p, q, ϵ ) 

end


@testset "init_particles" begin
 
     n = 40
     ϵ = 0.01
     rng = MersenneTwister(1234)
     
     sq = SquareParticles(rng, n, ϵ)
     @test true
     hs = HardSpheres(rng, n, ϵ)
     @test true

end

include("test_tempscoll.jl")
include("test_fantome.jl")
include("Aqua.jl")
