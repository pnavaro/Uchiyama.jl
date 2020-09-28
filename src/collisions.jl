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

function reset!(pc :: ParticleCollisions, i)
    pc.dt[i, :] .= Inf
    pc.dt[:, i] .= Inf
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

    function BoxCollisions( p :: SquareParticles )

        dt = zeros(p.n, 4)
        fill!(dt, Inf)

        for k in 1:p.n
            t, m = compute_dt(p, k)
            dt[k, m] = t
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
            pc.dt[i, j] = compute_dt(particles, j, i)
        end
    end
end

function dt_min_position( collisions )
    p = argmin(collisions.dt)
    dt_min = collisions.dt[p]
    return dt_min, p[1], p[2]
end

