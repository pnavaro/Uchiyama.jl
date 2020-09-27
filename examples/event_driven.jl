# -*- coding: utf-8 -*-
using Plots
using Random
using ProgressMeter
using Revise


# +
using Uchiyama

n = 50 # number of particles
ϵ = 0.02
c = trunc(Int, 200ϵ) # marker size

rng = MersenneTwister(1234)

q, v = box_particles(rng, n, ϵ)
particles = ParticleCollisions(n, q, v, ϵ)
walls = WallCollisions(n, q, v, ϵ)

# +

steps = 500
pbar = Progress(steps)

anim = @animate for _ in 1:steps
    
     dt = step!(n, ϵ, q, v, particles, walls)

     p = plot(size  = (200,200), 
              xlims = (-0.5,0.5), 
              ylims = (-0.5,0.5), 
              grid  = false, 
              axis  = nothing, framestyle= :box, legend=false, widen = false)
     
     scatter!( getindex.(q,1) .- getindex.(q,2), 
               getindex.(q,2) .+ getindex.(q,1), 
               markershape  = :square, 
               markersize   = c, 
               aspect_ratio = :equal)
     next!(pbar)

end
# -

gif(anim, joinpath(@__DIR__, "event_driven.gif"), fps = 20)


