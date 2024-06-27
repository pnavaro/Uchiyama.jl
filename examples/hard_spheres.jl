# -*- coding: utf-8 -*-
# # Hard spheres in a box
#
# md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/hard_spheres.ipynb)
# md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/hard_spheres.ipynb)

using Plots
using Random
using ProgressMeter
using Uchiyama


# +


n = 50 # number of particles
ϵ = 0.02
c = trunc(Int, 200ϵ) # marker size

rng = MersenneTwister(1234)

hs = HardSpheres(rng, n, ϵ)
pc = ParticleCollisions(hs)
bc = BoxCollisions(hs)

# +

steps = 1000
pbar = Progress(steps)

anim = @animate for _ in 1:steps
    
    dt = step!(hs, pc, bc)

    p = plot(size  = (200,200), 
             xlims = (0,1), 
             ylims = (0,1), 
             grid  = false, 
             axis  = nothing, 
             framestyle= :box, 
             legend=false, 
             widen = false)
    
    scatter!( getindex.(hs.q,1), 
              getindex.(hs.q,2), 
              markershape  = :circle, 
              markersize   = c, 
              aspect_ratio = :equal)

    next!(pbar)

end every 10
# -

gif(anim, joinpath(@__DIR__, "hard_spheres.gif"), fps = 10)


