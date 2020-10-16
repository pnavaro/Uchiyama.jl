var documenterSearchIndex = {"docs":
[{"location":"periodic_hard_spheres/","page":"Hard spheres periodic","title":"Hard spheres periodic","text":"EditURL = \"https://github.com/pnavaro/Uchiyama.jl/blob/master/examples/periodic_hard_spheres.jl\"","category":"page"},{"location":"periodic_hard_spheres/","page":"Hard spheres periodic","title":"Hard spheres periodic","text":"-- coding: utf-8 --","category":"page"},{"location":"periodic_hard_spheres/","page":"Hard spheres periodic","title":"Hard spheres periodic","text":"using Plots\nusing Random\nusing ProgressMeter\nusing Revise","category":"page"},{"location":"periodic_hard_spheres/","page":"Hard spheres periodic","title":"Hard spheres periodic","text":"","category":"page"},{"location":"periodic_hard_spheres/","page":"Hard spheres periodic","title":"Hard spheres periodic","text":"using Uchiyama\n\nn = 100 # number of particles\nϵ = 0.02\nc = trunc(Int, 200ϵ) # marker size\n\nrng = MersenneTwister(1234)\n\nhs = HardSpheres(rng, n, ϵ)\npc = PeriodicCollisions(hs)","category":"page"},{"location":"periodic_hard_spheres/","page":"Hard spheres periodic","title":"Hard spheres periodic","text":"","category":"page"},{"location":"periodic_hard_spheres/","page":"Hard spheres periodic","title":"Hard spheres periodic","text":"steps = 1000\npbar = Progress(steps)\n\nanim = @animate for _ in 1:steps\n\n    dt = step!(hs, pc)\n\n    p = plot(size  = (200,200),\n             xlims = (0,1),\n             ylims = (0,1),\n             grid  = false,\n             axis  = nothing,\n             framestyle= :none,\n             legend=false,\n             widen = false)\n\n    scatter!( getindex.(hs.q,1),\n              getindex.(hs.q,2),\n              markershape  = :circle,\n              markersize   = c,\n              aspect_ratio = :equal)\n\n    next!(pbar)\n\nend every 10","category":"page"},{"location":"periodic_hard_spheres/","page":"Hard spheres periodic","title":"Hard spheres periodic","text":"","category":"page"},{"location":"periodic_hard_spheres/","page":"Hard spheres periodic","title":"Hard spheres periodic","text":"gif(anim, joinpath(@__DIR__, \"periodic_hard_spheres.gif\"), fps = 20)","category":"page"},{"location":"periodic_hard_spheres/","page":"Hard spheres periodic","title":"Hard spheres periodic","text":"","category":"page"},{"location":"periodic_hard_spheres/","page":"Hard spheres periodic","title":"Hard spheres periodic","text":"This page was generated using Literate.jl.","category":"page"},{"location":"hard_spheres/","page":"Hard spheres in box","title":"Hard spheres in box","text":"EditURL = \"https://github.com/pnavaro/Uchiyama.jl/blob/master/examples/hard_spheres.jl\"","category":"page"},{"location":"hard_spheres/#Hard-spheres-in-a-box","page":"Hard spheres in box","title":"Hard spheres in a box","text":"","category":"section"},{"location":"hard_spheres/","page":"Hard spheres in box","title":"Hard spheres in box","text":"(Image: ) (Image: )","category":"page"},{"location":"hard_spheres/","page":"Hard spheres in box","title":"Hard spheres in box","text":"using Plots\nusing Random\nusing ProgressMeter\nusing Revise","category":"page"},{"location":"hard_spheres/","page":"Hard spheres in box","title":"Hard spheres in box","text":"","category":"page"},{"location":"hard_spheres/","page":"Hard spheres in box","title":"Hard spheres in box","text":"using Uchiyama\n\nn = 50 # number of particles\nϵ = 0.02\nc = trunc(Int, 200ϵ) # marker size\n\nrng = MersenneTwister(1234)\n\nhs = HardSpheres(rng, n, ϵ)\npc = ParticleCollisions(hs)\nbc = BoxCollisions(hs)","category":"page"},{"location":"hard_spheres/","page":"Hard spheres in box","title":"Hard spheres in box","text":"","category":"page"},{"location":"hard_spheres/","page":"Hard spheres in box","title":"Hard spheres in box","text":"steps = 1000\npbar = Progress(steps)\n\nanim = @animate for _ in 1:steps\n\n    dt = step!(hs, pc, bc)\n\n    p = plot(size  = (200,200),\n             xlims = (0,1),\n             ylims = (0,1),\n             grid  = false,\n             axis  = nothing,\n             framestyle= :box,\n             legend=false,\n             widen = false)\n\n    scatter!( getindex.(hs.q,1),\n              getindex.(hs.q,2),\n              markershape  = :circle,\n              markersize   = c,\n              aspect_ratio = :equal)\n\n    next!(pbar)\n\nend every 10","category":"page"},{"location":"hard_spheres/","page":"Hard spheres in box","title":"Hard spheres in box","text":"","category":"page"},{"location":"hard_spheres/","page":"Hard spheres in box","title":"Hard spheres in box","text":"gif(anim, joinpath(@__DIR__, \"hard_spheres.gif\"), fps = 10)","category":"page"},{"location":"hard_spheres/","page":"Hard spheres in box","title":"Hard spheres in box","text":"","category":"page"},{"location":"hard_spheres/","page":"Hard spheres in box","title":"Hard spheres in box","text":"This page was generated using Literate.jl.","category":"page"},{"location":"periodic_uchiyama/","page":"Uchiyama model periodic","title":"Uchiyama model periodic","text":"EditURL = \"https://github.com/pnavaro/Uchiyama.jl/blob/master/examples/periodic_uchiyama.jl\"","category":"page"},{"location":"periodic_uchiyama/","page":"Uchiyama model periodic","title":"Uchiyama model periodic","text":"-- coding: utf-8 --","category":"page"},{"location":"periodic_uchiyama/","page":"Uchiyama model periodic","title":"Uchiyama model periodic","text":"using Plots\nusing Random\nusing ProgressMeter\nusing Revise","category":"page"},{"location":"periodic_uchiyama/","page":"Uchiyama model periodic","title":"Uchiyama model periodic","text":"","category":"page"},{"location":"periodic_uchiyama/","page":"Uchiyama model periodic","title":"Uchiyama model periodic","text":"using Uchiyama\n\nn = 50 # number of particles\nϵ = 0.02\nc = trunc(Int, 200ϵ) # marker size\n\nrng = MersenneTwister(1234)\n\nsquares = SquareParticles(rng, n, ϵ)\n\ncollisions = PeriodicCollisions(squares)\n\nsteps = 500\npbar = Progress(steps)\n\nanim = @animate for _ in 1:steps\n\n     dt = step!(squares, collisions)\n\n     p = plot(size  = (200,200),\n              xlims = (0,1),\n              ylims = (0,1),\n              grid  = false,\n              axis  = nothing, legend=false, framestyle = :box, widen = false)\n\n     scatter!( getindex.(squares.q,1),\n               getindex.(squares.q,2),\n               markershape  = :diamond,\n               markersize   = c,\n               aspect_ratio = :equal)\n     next!(pbar)\n\nend","category":"page"},{"location":"periodic_uchiyama/","page":"Uchiyama model periodic","title":"Uchiyama model periodic","text":"","category":"page"},{"location":"periodic_uchiyama/","page":"Uchiyama model periodic","title":"Uchiyama model periodic","text":"gif(anim, joinpath(@__DIR__, \"periodic_uchiyama.gif\"), fps = 20)","category":"page"},{"location":"periodic_uchiyama/","page":"Uchiyama model periodic","title":"Uchiyama model periodic","text":"","category":"page"},{"location":"periodic_uchiyama/","page":"Uchiyama model periodic","title":"Uchiyama model periodic","text":"This page was generated using Literate.jl.","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = Uchiyama","category":"page"},{"location":"#Uchiyama","page":"Home","title":"Uchiyama","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Uchiyama’s particle motion is computed using an Event-Driven algorithm. The Event-Driven method is concerned with the times at which events, in this case collisions, take place. The algorithm is the following: you establish a list of the collisions to come if the particles moved only in straight lines. This list is ordered according to the time at which the collisions will take place, the first element being the closest collision in time. The particles are then displaced until this time and the collision is carried out. The list of future collisions is then updated. Indeed, some of them may be invalidated by the collision that has just taken place since the velocities and therefore the direction of the two particles involved have changed, new ones may also become possible for the same reasons. Then, again, we move to the nearest collision in time and so on ...","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [Uchiyama]","category":"page"},{"location":"event_driven/","page":"Uchiyama model in box","title":"Uchiyama model in box","text":"EditURL = \"https://github.com/pnavaro/Uchiyama.jl/blob/master/examples/event_driven.jl\"","category":"page"},{"location":"event_driven/","page":"Uchiyama model in box","title":"Uchiyama model in box","text":"-- coding: utf-8 --","category":"page"},{"location":"event_driven/","page":"Uchiyama model in box","title":"Uchiyama model in box","text":"using Plots\nusing Random\nusing ProgressMeter\nusing Revise","category":"page"},{"location":"event_driven/","page":"Uchiyama model in box","title":"Uchiyama model in box","text":"","category":"page"},{"location":"event_driven/","page":"Uchiyama model in box","title":"Uchiyama model in box","text":"using Uchiyama\n\nn = 50 # number of particles\nϵ = 0.02\nc = trunc(Int, 200ϵ) # marker size\n\nrng = MersenneTwister(1234)\n\nparticles = SquareParticles(rng, n, ϵ, option = :box)\npc = ParticleCollisions(particles)\nbc = BoxCollisions(particles)","category":"page"},{"location":"event_driven/","page":"Uchiyama model in box","title":"Uchiyama model in box","text":"","category":"page"},{"location":"event_driven/","page":"Uchiyama model in box","title":"Uchiyama model in box","text":"steps = 1000\npbar = Progress(steps)\n\nanim = @animate for _ in 1:steps\n\n     dt = step!(particles, pc, bc)\n\n     p = plot(size  = (200,200),\n              xlims = (-0.5,0.5),\n              ylims = (-0.5,0.5),\n              grid  = false,\n              axis  = nothing, framestyle= :box, legend=false, widen = false)\n\n     scatter!( getindex.(particles.q,1) .- getindex.(particles.q,2),\n               getindex.(particles.q,2) .+ getindex.(particles.q,1),\n               markershape  = :square,\n               markersize   = c,\n               aspect_ratio = :equal)\n     next!(pbar)\n\nend every 10","category":"page"},{"location":"event_driven/","page":"Uchiyama model in box","title":"Uchiyama model in box","text":"","category":"page"},{"location":"event_driven/","page":"Uchiyama model in box","title":"Uchiyama model in box","text":"gif(anim, joinpath(@__DIR__, \"event_driven.gif\"), fps = 10)","category":"page"},{"location":"event_driven/","page":"Uchiyama model in box","title":"Uchiyama model in box","text":"","category":"page"},{"location":"event_driven/","page":"Uchiyama model in box","title":"Uchiyama model in box","text":"This page was generated using Literate.jl.","category":"page"}]
}
