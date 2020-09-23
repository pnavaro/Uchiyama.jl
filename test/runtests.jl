using Plots
using Random
using Test

include("../src/periodic_uchiyama.jl")

#@testset "run" begin

    n = 10
    ϵ = 0.1
    q = [[0.5       ,0.5       ],
         [0.20006862,0.38139954],
         [0.32267135,0.72687046],
         [0.68044674,0.78095695],
         [0.77473372,0.51989917],
         [0.4949439 ,0.23201753],
         [0.74204115,0.28248482],
         [0.63101615,0.3650337 ],
         [0.5846235 ,0.64552143],
         [0.45656816,0.76644607]]

    v = [[ 0, 1],
         [ 1, 0],
         [-1, 0],
         [ 0,-1],
         [ 1, 0],
         [ 0, 1],
         [ 0, 1],
         [-1, 0],
         [-1, 0],
         [ 0,-1]]

    collisions = PCollisionMatrix(n, q, v, ϵ)
    step!(n, ϵ, q, v, collisions)

#    anim = @animate for _ in 1:100
#
#        dt = step!(n, ϵ, q, v, collisions)
#        p = plot(size=(500,500), xlims=(0,1), ylims=(0,1), 
#        grid=false, axis=nothing, legend=false, framestyle = :box, widen = false)
#        scatter!(getindex.(q,1), getindex.(q,2), markershape = :square, 
#                 markersize=50, aspect_ratio=:equal)
#
#    end
#
#    gif(anim, joinpath(@__DIR__, "periodic_uchiyama.gif"), fps = 15)

#end
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


@testset "init_particles" begin
 
     n = 10
     ϵ = 1 / n
     c = trunc(Int, 1000 * 2ϵ)
     rng = MersenneTwister(1234)
     
     q, v = init_particles(n, ϵ, rng)

     plot(size=(500,500), xlims=(0,1), ylims=(0,1), 
          grid=false, axis=nothing, legend=false, framestyle = :box, widen = false)
     scatter!(getindex.(q,1), getindex.(q,2), markershape = :square, markersize=c, aspect_ratio=:equal)
     savefig("particles.png")
 
end

@testset "tempscoll" begin

   ϵ = 0.1 

   for p in 0.2:0.1:0.8
       x = [p, 0.2]; v = [0.0,  1.0]
       y = [p, 0.8]; w = [0.0, -1.0]
       @test tempscoll(x, y, v, w, ϵ) ≈ 0.2
       x = [0.2, p]; v = [ 1.0, 0.0]
       y = [0.8, p]; w = [-1.0, 0.0]
       @test tempscoll(x, y, v, w, ϵ) ≈ 0.2
   end

   x = [0.8,0.8]; v = [-1.0,-1.0]
   y = [0.2,0.2]; w = [+1.0,1.0]
   @test tempscoll(x, y, v, w, ϵ) ≈ 0.5
   x = [0.2,0.8]; v = [+1.0,-1.0]
   y = [0.8,0.2]; w = [-1.0,1.0]
   @test tempscoll(x, y, v, w, ϵ) ≈ 0.5

end
tempscoll(q[l] + offset[i], q[k], v[l], v[k], ϵ)
=#
