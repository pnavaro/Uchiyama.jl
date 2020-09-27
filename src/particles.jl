abstract type Particles end

function compute_dt!(dt, i, particles :: Particles)

    q_i = particles.q[i]
    v_i = particles.v[i]
    for j in 1:particles.n
        if i != j
            tcoll = particles(particles.q[j], q_i, particles.v[j], v_i)
            if !isinf(tcoll)
                dt[i, j] = tcoll
            end
        end
    end
end


export HardSpheres

struct HardSpheres <: Particles

    n
    q
    v
    ϵ

    function HardSpheres(rng, n, ϵ)
    
        q = Vector{Float64}[]
    
        J4 = [1 0; 0 1]
    
        push!(q, [0.5000, 0.5000])
        push!(q, [0.1622, 0.4505])
        push!(q, [0.3112, 0.2290])
        push!(q, [0.5285, 0.9133])
        push!(q, [0.1656, 0.1524])
        push!(q, [0.6020, 0.8258])
        push!(q, [0.2630, 0.5383])
        push!(q, [0.3000, 0.8000])
        push!(q, [0.6892, 0.0782])
        push!(q, [0.7482, 0.4427])
    
        p = zeros(2)
    
        for klm in 11:n
            frein = 0
            overlap = 1
            while (overlap == 1) && frein < 10000
    
                p0 = [2ϵ + (1 - 4ϵ) * rand(rng), 2ϵ + (1 - 4ϵ) * rand(rng)]
                p = J4 * p0
                overlap2 = 0
                for subh in 1:(klm-1)
                    if norm(p - q[subh], 2) < 2ϵ
                        overlap2 = 1
                    end
                end
    
                overlap = overlap2
                frein += 1
            end
    
            frein == 10000 && (@error "Echec tirage initial")
    
            push!(q, p)
    
        end
    
        v = [ randn(rng, 2) for j in 1:n ]
    
        new( n, q, v, ϵ )
    
    end

end

function (p :: HardSpheres)(x, y, v, w)

    ϵ = p.ϵ
    δv = v - w
    δz = x - y
    c = (δv'δz).^2 - (δv'δv) * (δz'δz - (2ϵ)^2)

    if (δv'δz) >= 0
        δt = Inf
    elseif c < 0
        δt = Inf
    else
        δt = -(δv'δz + sqrt(c)) / (δv'δv)
    end
    return δt

end


function (hs :: HardSpheres)(x, v)

    ϵ = hs.ϵ

    if v[1] > 0
        s1 = (1 - ϵ - x[1]) / v[1]
    elseif v[1] < 0
        s1 = (ϵ - x[1]) / v[1]
    elseif v[1] == 0
        s1 = Inf
    end

    if v[2] > 0
        s2 = (1 - ϵ - x[2]) / v[2]
    elseif v[2] < 0
        s2 = (ϵ - x[2]) / v[2]
    elseif v[2] == 0
        s2 = Inf
    end

    r = min(s1, s2)

    if r == s1
        s = s1
        i = 1
    else
        s = s2
        i = 2
    end

    return s, i

end
