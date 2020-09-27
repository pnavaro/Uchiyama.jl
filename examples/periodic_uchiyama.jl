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

squares = Squares(rng, n, ϵ)

collisions = PeriodicCollisions(squares)

steps = 500
pbar = Progress(steps)

anim = @animate for _ in 1:steps
    
     dt = step!(squares, collisions)

     p = plot(size  = (200,200), 
              xlims = (0,1), 
              ylims = (0,1), 
              grid  = false, 
              axis  = nothing, legend=false, framestyle = :box, widen = false)

     scatter!( getindex.(squares.q,1), 
               getindex.(squares.q,2), 
               markershape  = :diamond, 
               markersize   = c, 
               aspect_ratio = :equal)
     next!(pbar)

end
# -

gif(anim, joinpath(@__DIR__, "periodic_uchiyama.gif"), fps = 20)
