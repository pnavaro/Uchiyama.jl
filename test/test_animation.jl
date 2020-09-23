using Plots
using Random
using Test

include("../src/periodic_uchiyama.jl")

n = 50
ϵ = 1 / n
c = trunc(Int, 500ϵ)
rng = MersenneTwister(1234)
     
q, v = init_particles(n, ϵ, rng)

plot(size=(500,500), xlims=(0,1), ylims=(0,1), 
      grid=false, axis=nothing, legend=false, framestyle = :box, widen = false)
  scatter!(getindex.(q,1), getindex.(q,2), markershape = :square, markersize=c, aspect_ratio=:equal)
savefig("particles.png")
 
collisions = PCollisionMatrix(n, q, v, ϵ)

anim = @animate for _ in 1:100

     dt = step!(n, ϵ, q, v, collisions)
     p = plot(size=(500,500), xlims=(0,1), ylims=(0,1), 
     grid=false, axis=nothing, legend=false, framestyle = :box, widen = false)
     scatter!(getindex.(q,1), getindex.(q,2), markershape = :square, 
              markersize=c, aspect_ratio=:equal)

end

gif(anim, joinpath(@__DIR__, "periodic_uchiyama.gif"), fps = 15)

