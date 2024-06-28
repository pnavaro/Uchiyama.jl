# # Uchiyama model in box

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/event_driven.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/event_driven.ipynb)

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

particles = SquareParticles(rng, n, ϵ, option = :box)
pc = ParticleCollisions(particles)
bc = BoxCollisions(particles)

steps = 1000

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

end every 10

gif(anim, joinpath(@__DIR__, "event_driven.gif"), fps = 10)
