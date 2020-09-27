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


function hscompute_dt(i, n, q, v, ϵ, dt)

    q_i = q[i]
    v_i = v[i]
    for j in 1:n
        if i != j
            tcoll = hstempscoll(q[j], q_i, v[j], v_i, ϵ)
            if !isinf(tcoll)
                dt[i, j] = tcoll
            end
        end
    end
end

export PCollision

struct PCollision

    dt

    function PCollision( n, q, v, ϵ)

        dt = zeros(n, n)
        fill!(dt, Inf)

        for k in 1:n
            for l in (k+1):n
                dt[k, l] = hstempscoll(q[l], q[k], v[l], v[k], ϵ)
            end
        end
        new( dt )
    end
end

export WCollision

struct WCollision

    dt

    function WCollision( n, q, v, ϵ)

        dt = zeros(n, 2)
        fill!(dt, Inf)

        for k in 1:n
            t, m = hstempscollmur(q[k], v[k], ϵ)
            dt[k, m] = t
        end

        new(dt)

    end

end

const hs_wall_rebound = [
    v -> [-v[1],  v[2]],
    v -> [ v[1], -v[2]]
]


export hs_particles

function hs_particles(rng, n, ϵ)

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

    return q, v

end

function step!(n, ϵ, q, v, collisions::PCollision, walls::WCollision)


    tempsp, i1, i2 = dt_min_position(collisions)

    tempsm, j1, j2 = dt_min_position(walls)

    if tempsp < tempsm

        for i in 1:n
            q[i] = q[i] + tempsp * v[i]
        end

        J = ((v[i2] - v[i1])'*(q[i2] - q[i1])) / 2ϵ
        v[i1] = v[i1] + J * (q[i2] - q[i1]) / 2ϵ
        v[i2] = v[i2] - J * (q[i2] - q[i1]) / 2ϵ

        collisions.dt .-= tempsp
        reset!(collisions.dt, i1)
        reset!(collisions.dt, i2)

        walls.dt .-= tempsp
        walls.dt[i1,:] .= Inf
        walls.dt[i2,:] .= Inf

        hscompute_dt(i1, n, q, v, ϵ, collisions.dt)
        hscompute_dt(i2, n, q, v, ϵ, collisions.dt)

        t, i = hstempscollmur(q[i1], v[i1], ϵ)
        !isinf(t) && ( walls.dt[i1, i] = t)

        t, i = hstempscollmur(q[i2], v[i2], ϵ)
        !isinf(t) && ( walls.dt[i2, i] = t)

    else

        for i in 1:n
            q[i] = q[i] + tempsm * v[i]
        end

        v[j1] = hs_wall_rebound[j2](v[j1])

        collisions.dt .-= tempsm
        reset!(collisions.dt, j1)

        walls.dt .-= tempsm
        walls.dt[j1,:] .= Inf

        hscompute_dt(j1, n, q, v, ϵ, collisions.dt)

        t, i = hstempscollmur(q[j1], v[j1], ϵ)
        !isinf(t) && (walls.dt[j1, i] = t)

    end

    return min(tempsp, tempsm)

end
