using Plots
using Uchiyama
using Random


n = 50 # number of particles
ϵ = 1 / n
c = trunc(Int, 500ϵ) # marker size

rng = MersenneTwister(1234)
     
q, v = init_particles(rng, n, ϵ)

collisions = PeriodicCollisions(n, q, v, ϵ)

anim = @animate for _ in 1:500

     dt = step!(n, ϵ, q, v, collisions)

     p = plot(size  = (500,500), 
              xlims = (0,1), 
              ylims = (0,1), 
              grid  = false, 
              axis  = nothing, legend=false, framestyle = :box, widen = false)

     scatter!( getindex.(q,1), 
               getindex.(q,2), 
               markershape  = :diamond, 
               markersize   = c, 
               aspect_ratio = :equal)

end

gif(anim, joinpath(@__DIR__, "periodic_uchiyama.gif"), fps = 20)
