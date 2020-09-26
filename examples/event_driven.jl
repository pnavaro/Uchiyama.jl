# -*- coding: utf-8 -*-
using Plots
using Random
using ProgressMeter
using Revise


# +
using Uchiyama

n = 40 # number of particles
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
              axis  = nothing, framestyle= :none, legend=false, widen = false)
     
     plot!( [x->x+0.5, x->-x+0.5, x->x-0.5, x->-x-0.5], lc = [0,0,0,0])

     scatter!( getindex.(q,1), 
               getindex.(q,2), 
               markershape  = :diamond, 
               markersize   = c, 
               aspect_ratio = :equal)
     next!(pbar)

end
# -

gif(anim, joinpath(@__DIR__, "event_driven.gif"), fps = 20)


