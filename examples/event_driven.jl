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

particles = SquareParticles(rng, n, ϵ, option = :box)
pc = ParticleCollisions(particles)
bc = BoxCollisions(particles)

# +

steps = 1000
pbar = Progress(steps)

anim = @animate for _ in 1:steps
    
     dt = step!(particles, pc, bc)

     p = plot(size  = (200,200), 
              xlims = (-0.5,0.5), 
              ylims = (-0.5,0.5), 
              grid  = false, 
              axis  = nothing, framestyle= :box, legend=false, widen = false)
     
     scatter!( getindex.(particles.q,1) .- getindex.(particles.q,2), 
               getindex.(particles.q,2) .+ getindex.(particles.q,1), 
               markershape  = :square, 
               markersize   = c, 
               aspect_ratio = :equal)
     next!(pbar)

end every 10
# -

gif(anim, joinpath(@__DIR__, "event_driven.gif"), fps = 10)


