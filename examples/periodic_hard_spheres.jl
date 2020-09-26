using Plots
using Random
using ProgressMeter
using Revise


# +
using Uchiyama

n = 200 # number of particles
ϵ = 4 / n
c = trunc(Int, 200ϵ) # marker size

rng = MersenneTwister(1234)

q, v = hs_particles(rng, n, ϵ)
pc = PCollisionMatrix(n, q, v, ϵ)

# +

steps = 100
pbar = Progress(steps)

anim = @animate for _ in 1:steps
    
    dt = step!(n, ϵ, q, v, pc)

    p = plot(size  = (200,200), 
             xlims = (0,1), 
             ylims = (0,1), 
             grid  = false, 
             axis  = nothing, 
             framestyle= :none, 
             legend=false, 
             widen = false)
    
    scatter!( getindex.(q,1), 
              getindex.(q,2), 
              markershape  = :circle, 
              markersize   = c, 
              aspect_ratio = :equal)

    next!(pbar)

end
# -

gif(anim, joinpath(@__DIR__, "periodic_hard_spheres.gif"), fps = 20)
