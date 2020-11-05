\begin{section}{title="About this Package", name="About"}

\lead{Uchiyama.jl is based upon  [Nathalie Ayi](https://www.ljll.math.upmc.fr/~ayi)'s work}

You have two types of particle:

* Uchiyama particle model
* Hard sphere

and two types of boundary condition:

* walls
* periodic

\end{section}


<!-- ==============================
     GETTING STARTED
     ============================== -->
\begin{section}{title="Getting started"}

In order to get started, just add the package (with **Julia ≥ 1.3**) and

```julia-repl
julia> pkg'add https://github.com/pnavaro/Uchiyama.jl'
julia> using Uchiyama
```

\end{section}

\begin{section}{title="Example"}

\lead{
    50 particles moving in a square box
}

```
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

```

\end{section}
