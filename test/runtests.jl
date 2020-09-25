using Test
using Random
using LinearAlgebra

@testset "step" begin

   import Uchiyama: PeriodicCollisions, step!

   n = 2
   ϵ = 0.1
   q = [[0.5,0.2], [0.5,0.8]]
   v = [[0,1], [0,-1]]
   collisions = PeriodicCollisions(n, q, v, ϵ)

   @test step!(n, ϵ, q, v, collisions) ≈ 0.2

   @test q[1] ≈ [0.5,0.4] 
   @test v[1] ≈ [0,-1]
   @test q[2] ≈ [0.5,0.6] 
   @test v[2] ≈ [0, 1]

   @test step!(n, ϵ, q, v, collisions) ≈ 0.3

   @test q[1] ≈ [0.5,0.1] 
   @test v[1] ≈ [0,1]
   @test q[2] ≈ [0.5,0.9] 
   @test v[2] ≈ [0,-1]

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
 
     import Uchiyama: init_particles
     n = 40
     ϵ = 0.01
     rng = MersenneTwister(1234)
     
     q, v = init_particles(rng, n, ϵ)

end

include("test_tempscoll.jl")
include("test_fantome.jl")
