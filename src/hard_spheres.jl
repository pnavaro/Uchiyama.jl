function hstempscoll(x, y, v, w, ϵ)

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

function hstempscollmur(x, v, ϵ)

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


function compute_dt!(dt, i, hs :: HardSpheres)

    q_i = hs.q[i]
    v_i = hs.v[i]
    for j in 1:hs.n
        if i != j
            tcoll = hstempscoll(hs.q[j], q_i, hs.v[j], v_i, hs.ϵ)
            if !isinf(tcoll)
                dt[i, j] = tcoll
            end
        end
    end
end

export ParticleCollisions

struct ParticleCollisions

    dt

    function ParticleCollisions( hs :: HardSpheres)

        dt = zeros(hs.n, hs.n)
        fill!(dt, Inf)

        for k in 1:hs.n
            for l in (k+1):hs.n
                dt[k, l] = hstempscoll(hs.q[l], hs.q[k], hs.v[l], hs.v[k], hs.ϵ)
            end
        end
        new( dt )
    end
end

export BoxCollisions

struct BoxCollisions

    dt

    function BoxCollisions( hs :: HardSpheres)

        dt = zeros(hs.n, 2)
        fill!(dt, Inf)

        for k in 1:hs.n
            t, m = hstempscollmur(hs.q[k], hs.v[k], hs.ϵ)
            dt[k, m] = t
        end

        new(dt)

    end

end

const hs_wall_rebound = [
    v -> [-v[1],  v[2]],
    v -> [ v[1], -v[2]]
]



function step!(hs :: HardSpheres, collisions::ParticleCollisions, walls::BoxCollisions)

    ϵ = hs.ϵ
    n = hs.n

    tempsp, i1, i2 = dt_min_position(collisions)

    tempsm, j1, j2 = dt_min_position(walls)

    if tempsp < tempsm

        for i in 1:n
            hs.q[i] = hs.q[i] + tempsp * hs.v[i]
        end

        J = ((hs.v[i2] - hs.v[i1])'*(hs.q[i2] - hs.q[i1])) / 2ϵ
        hs.v[i1] = hs.v[i1] + J * (hs.q[i2] - hs.q[i1]) / 2ϵ
        hs.v[i2] = hs.v[i2] - J * (hs.q[i2] - hs.q[i1]) / 2ϵ

        collisions.dt .-= tempsp
        reset!(collisions.dt, i1)
        reset!(collisions.dt, i2)

        walls.dt .-= tempsp
        walls.dt[i1,:] .= Inf
        walls.dt[i2,:] .= Inf

        compute_dt!(collisions.dt, i1, hs)
        compute_dt!(collisions.dt, i2, hs)

        t, i = hstempscollmur(hs.q[i1], hs.v[i1], ϵ)
        !isinf(t) && ( walls.dt[i1, i] = t)

        t, i = hstempscollmur(hs.q[i2], hs.v[i2], ϵ)
        !isinf(t) && ( walls.dt[i2, i] = t)

    else

        for i in 1:n
            hs.q[i] = hs.q[i] + tempsm * hs.v[i]
        end

        hs.v[j1] = hs_wall_rebound[j2](hs.v[j1])

        collisions.dt .-= tempsm
        reset!(collisions.dt, j1)

        walls.dt .-= tempsm
        walls.dt[j1,:] .= Inf

        compute_dt!(collisions.dt, j1, hs)

        t, i = hstempscollmur(hs.q[j1], hs.v[j1], ϵ)
        !isinf(t) && (walls.dt[j1, i] = t)

    end

    return min(tempsp, tempsm)

end
