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

q, v = hs_particles(rng, n, ϵ)
pc = PCollision(n, q, v, ϵ)
wc = WCollision(n, q, v, ϵ)

# +

steps = 500
pbar = Progress(steps)

anim = @animate for _ in 1:steps
    
    dt = step!(n, ϵ, q, v, pc, wc)

    p = plot(size  = (200,200), 
             xlims = (0,1), 
             ylims = (0,1), 
             grid  = false, 
             axis  = nothing, 
             framestyle= :box, 
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

gif(anim, joinpath(@__DIR__, "hard_spheres.gif"), fps = 20)
