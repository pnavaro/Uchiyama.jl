# -*- coding: utf-8 -*-
using Random
using Plots
using Distances
using ProgressMeter

# +
function time_driven(nstep)

    np = 50
    rng = MersenneTwister(1234)
    q = -0.5 .+ rand(2,np)
    v = randn(2,np)

    dt = 0.01

    bounds = [-2, 2, -2, 2]
    ϵ = 0.04

    bar = Progress(nstep)

    anim = @animate for step in 1:nstep
        
        next!(bar)
        
        q .+= dt .* v

        
        d = pairwise(Euclidean(), q)
        p = findall(x -> x < 2ϵ, d)

        ind1 = getindex.(p,1) 
        ind2 = getindex.(p,2)

        unique = ind1 .< ind2

        ind1 = ind1[unique]
        ind2 = ind2[unique]


        for (i1, i2) in zip(ind1, ind2)

            q1 = q[:, i1]
            q2 = q[:, i2]

            v1 = v[:, i1]
            v2 = v[:, i2]

            q_rel = q1 - q2
            v_rel = v1 - v2

            v_cm = 0.5 * (v1+v2) 

            qq_rel = q_rel'q_rel
            vq_rel = v_rel'q_rel
            v_rel = 2 * q_rel * vq_rel / qq_rel - v_rel

            v[:, i1] .= v_cm + v_rel * 0.5
            v[:, i2] .= v_cm - v_rel * 0.5

        end

        crossed_x1 = (q[1,:] .< (bounds[1] + ϵ))
        crossed_x2 = (q[1,:] .> (bounds[2] - ϵ))
        crossed_y1 = (q[2,:] .< (bounds[3] + ϵ))
        crossed_y2 = (q[2,:] .> (bounds[4] - ϵ))

        q[1, crossed_x1] .= bounds[1] + ϵ
        q[1, crossed_x2] .= bounds[2] - ϵ

        q[2, crossed_y1] .= bounds[3] + ϵ
        q[2, crossed_y2] .= bounds[4] - ϵ

        v[1, crossed_x1] .*= -1
        v[1, crossed_x2] .*= -1
        v[2, crossed_y1] .*= -1
        v[2, crossed_y2] .*= -1

        plt = plot(size  = (200,200), 
                   xlims = (-2,2), 
                   ylims = (-2,2), 
                   grid  = false, 
                   axis  = nothing, 
                   framestyle= :box, 
                   legend=false, 
                   widen = false)
    
        c = trunc(Int, 100ϵ)
        scatter!( view(q,1,:), view(q,2,:), markershape  = :circle, 
              markersize   = c , 
              aspect_ratio = :equal)


    end every 10

    gif(anim, joinpath(@__DIR__, "time_driven.gif"), fps = 10)

end
# -

@time time_driven(1000)


