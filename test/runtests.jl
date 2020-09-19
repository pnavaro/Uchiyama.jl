using Uchiyama
using Test
using Random

@testset "tempscoll.jl" begin

   epsilon = 0.1 

   for p in 0.2:0.1:0.8
       x = [p, 0.2]; v = [0.0,  1.0]
       y = [p, 0.8]; w = [0.0, -1.0]
       @test Uchiyama.tempscoll(x, y, v, w, epsilon) ≈ 0.2
       x = [0.2, p]; v = [ 1.0, 0.0]
       y = [0.8, p]; w = [-1.0, 0.0]
       @test Uchiyama.tempscoll(x, y, v, w, epsilon) ≈ 0.2
   end

   x = [0.8,0.8]; v = [-1.0,-1.0]
   y = [0.2,0.2]; w = [+1.0,1.0]
   @test Uchiyama.tempscoll(x, y, v, w, epsilon) ≈ 0.5
   x = [0.2,0.8]; v = [+1.0,-1.0]
   y = [0.8,0.2]; w = [-1.0,1.0]
   @test Uchiyama.tempscoll(x, y, v, w, epsilon) ≈ 0.5

end

# @testset "Uchiyama.jl" begin
# 
#     n = 40
#     epsilon = 1 / n
#     
#     sim  = Trajectory(n, epsilon)
#     visu = Visualisation(sim, 1000)
#     anim(visu)
# 
# end
