using Plots
using Random
using ProgressMeter
using Revise


# +
using Uchiyama

n = 100 # number of particles
ϵ = 0.02
c = trunc(Int, 200ϵ) # marker size

rng = MersenneTwister(1234)

hs = HardSpheres(rng, n, ϵ)
pc = PeriodicCollisions(hs)

# +

steps = 100
pbar = Progress(steps)

anim = @animate for _ in 1:steps
    
    dt = step!(hs, pc)

    p = plot(size  = (200,200), 
             xlims = (0,1), 
             ylims = (0,1), 
             grid  = false, 
             axis  = nothing, 
             framestyle= :none, 
             legend=false, 
             widen = false)
    
    scatter!( getindex.(hs.q,1), 
              getindex.(hs.q,2), 
              markershape  = :circle, 
              markersize   = c, 
              aspect_ratio = :equal)

    next!(pbar)

end
# -

gif(anim, joinpath(@__DIR__, "periodic_hard_spheres.gif"), fps = 20)
