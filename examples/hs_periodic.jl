using Plots
using LinearAlgebra

function main( np, nstep )

    ϵ = 1/np
    
    q = [zeros(2) for i in 1:np]
    q[1] = [0.5,0.5]
    
    for klm = 2:np
        frein = 0
        overlap = 1
        p = zeros(2)
        while ((overlap == 1) && frein < 7000)
            p .= 2ϵ * ones(2) + (1 - 4ϵ) * rand(2)
            overlap2 = 0
            for subh = 1:(klm-1)
                if norm(p - q[subh],2) < 2ϵ
                    overlap2 = 1
                end
            end
            overlap = overlap2
            frein += 1
            if frein == 10000
                frein
                @error "frein trop grand"
                exit
            end
        end
        q[klm] = p
    end
    
    offset = [[0, 0], [1,  0], [-1,  0], 
              [0, 1], [0, -1], [-1,  1], 
              [1, 1], [1, -1], [-1, -1]]
    
    function compute_dt( qi, qj, vi, vj )
    
        δz = qi .- qj
        δv = vi .- vj
    
        c = (δv'δz).^2 - (δv'δv) * (δz'δz - (2ϵ)^2)
    
        if δv'δz >= 0
            δt = Inf
        elseif c < 0
            δt = Inf
        else
            δt = -(δv'δz + sqrt(c)) / (δv'δv)
        end
    
        return δt
    
    end
    
    v = [randn(2) for j = 1:np]
    
    Collisions = Inf .* ones(np,np)
    Fantome = zeros(Int, np,np)
    
    for k=1:np, l=k+1:np

        t = [compute_dt(q[l]+offset[i],q[k],v[l],v[k]) for i in 1:9]

        Collisions[k,l] = minimum(t)
        if isinf(Collisions[k,l])
            Fantome[k,l] = 0
        else
            Fantome[k,l] = argmin(t)
        end
    end

    time = 0
    
    anim = @animate for step in 1:nstep
    
        m, n = Tuple(argmin(Collisions))
        dt = Collisions[m,n]
        time += dt

        if Fantome[m,n] > 0
            qa = [q[n] + offset[k] + dt*v[n] for k in 1:9]
            qb = q[m] + dt .* v[m]
            vr = argmin([norm(qx-qb) for qx in qa])
    
            q[n] = qa[vr]
            q[m] = qb
    
            J = (dot(v[n]-v[m],q[n]-q[m]))/ 2ϵ
            v[m] = v[m] + J * (q[n]-q[m]) / 2ϵ
            v[n] = v[n] - J * (q[n]-q[m]) / 2ϵ
            q[n] = mod.(q[n],1)
            q[m] = mod.(q[m],1)
        end
         
        for i in [1:min(m,n)-1;min(m,n)+1:max(m,n)-1;max(m,n)+1:np]
            q[i] = q[i] .+ dt .* v[i]
            q[i] = mod.(q[i],1)
        end

        plt = plot(size  = (200,200), 
                   xlims = (0,1), 
                   ylims = (0,1), 
                   grid  = false, 
                   axis  = nothing, 
                   framestyle= :none, 
                   legend=false, 
                   widen = false)
    
        c = trunc(Int, 200ϵ)
        scatter!( getindex.(q,1), getindex.(q,2), markershape  = :circle, 
              markersize   = c , 
              aspect_ratio = :equal)

         
        Collisions .-= dt
        Collisions[m,:] .= Inf
        Collisions[n,:] .= Inf
        Collisions[:,m] .= Inf
        Collisions[:,n] .= Inf
        
        for j in [m,n], l in [1:j-1;j+1:np]

            t=[compute_dt(q[l]+offset[k],q[j],v[l],v[j]) for k in 1:9]
            Collisions[j,l]=minimum(t)
    
            if isinf(Collisions[m,l])
                Fantome[j,l]=0
            else 
                Fantome[j,l]=argmin(t)
            end
     
        end
        
    end

    gif(anim, joinpath(@__DIR__, "hs_periodic.gif"), fps = 10)

    return time

end

@time main(100, 1000)
