# # Uchiyama model periodic

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/periodic_uchiyama.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/periodic_uchiyama.ipynb)

import Pkg
Pkg.add(url="https://github.com/pnavaro/Uchiyama.jl")

#-

using Plots
using Random
using Uchiyama


n = 50 # number of particles
ϵ = 0.02
c = trunc(Int, 200ϵ) # marker size

rng = MersenneTwister(1234)

squares = SquareParticles(rng, n, ϵ)

collisions = PeriodicCollisions(squares)

steps = 500

anim = @animate for _ in 1:steps
    
     dt = step!(squares, collisions)

     p = plot(size  = (200,200), 
              xlims = (0,1), 
              ylims = (0,1), 
              grid  = false, 
              axis  = nothing, legend=false, framestyle = :none, widen = false)

     scatter!( getindex.(squares.q,1), 
               getindex.(squares.q,2), 
               markershape  = :diamond, 
               markersize   = c, 
               aspect_ratio = :equal)

end

gif(anim, joinpath(@__DIR__, "periodic_uchiyama.gif"), fps = 20)


