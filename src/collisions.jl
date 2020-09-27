export ParticleCollisions

struct ParticleCollisions

    dt

    function ParticleCollisions( p :: Particles)

        dt = zeros(p.n, p.n)
        fill!(dt, Inf)

        for k in 1:p.n
            for l in (k+1):p.n
                dt[k, l] = compute_dt(p, l, k)
            end
        end
        new( dt )
    end

end

export BoxCollisions

struct BoxCollisions

    dt

    function BoxCollisions( p :: HardSpheres)

        dt = zeros(p.n, 2)
        fill!(dt, Inf)

        for k in 1:p.n
            t, m = compute_dt(p, k)
            dt[k, m] = t
        end

        new(dt)

    end

    function BoxCollisions( p :: Squares )

        dt = zeros(p.n, 4)
        fill!(dt, Inf)

        for k in 1:p.n
            t, m = compute_dt(p, k)
            !isinf(t) && (dt[k, m] = t)
        end

        new(dt)

    end

end

function reset!(bc :: BoxCollisions, i)
    bc.dt[i, :] .= Inf
end

function compute_dt!(pc :: ParticleCollisions, i, particles :: Particles)

    for j in 1:particles.n
        if i != j
            tcoll = compute_dt(particles, j, i)
            if !isinf(tcoll)
                dt[i, j] = tcoll
            end
        end
    end
end

function dt_min_position( collisions )
    p = argmin(collisions.dt)
    dt_min = collisions.dt[p]
    return dt_min, p[1], p[2]
end


export PeriodicCollisions

struct PeriodicCollisions

    dt :: Array{Float64,2}
    fantome :: Array{Int, 2}

    function PeriodicCollisions( p :: Particles)

        #aij is stored in AP(i+j(j-1)/2) for $i \leq j$;

        dt = zeros(Float64, (p.n, p.n))
        fantome = zeros(Int, (p.n, p.n))
        fill!(dt, Inf)
        fill!(fantome, 0)

        for i in 1:p.n
            for j in (i+1):p.n
                dt_local = Inf
                k = 0
                while (isinf(dt_local) && k < 9)
                    k += 1
                    dt_local = compute_dt(p, j, i, k)
                end

                dt[i, j] = dt_local
                fantome[i, j] = k

            end
        end

        new( dt, fantome )

    end

end

function compute_dt!(pc :: PeriodicCollisions, i, p :: Particles)

    for j in 1:p.n
        if i != j
            dt_local = Inf
            k = 0
            while (isinf(dt_local) && k < 9)
                k += 1
                dt_local = compute_dt(p, j, i, k)
            end

            pc.dt[i,j] = dt_local
            pc.fantome[i,j] = k
        end
    end

end

function dt_min_position(pc :: PeriodicCollisions)
    p = argmin(pc.dt)
    dt_min = pc.dt[p]
    num_fant = pc.fantome[p]
    return dt_min, num_fant, p[1], p[2]
end

function reset!(pc :: PeriodicCollisions, i)
    pc.dt[i, :] .= Inf
    pc.dt[:, i] .= Inf
    pc.fantome[i, :] .= 0
    pc.fantome[:, i] .= 0
end
