
function tempscollmur(x, v, ϵ)
    t = Inf
    i = 1

    if x[2] > 0 && v[1] == 1
        t = (0.5 - (x[1] + x[2]) - ϵ) / (v[1] + v[2])
        i = 1
    elseif x[1] > 0 && v[2] == 1
        t = (0.5 - (x[1] + x[2]) - ϵ) / (v[1] + v[2])
        i = 1
    elseif x[2] > 0 && v[1] == -1
        t = (0.5 - ϵ + (x[1] - x[2])) / (v[2] - v[1])
        i = 2
    elseif x[1] < 0 && v[2] == 1
        t = (0.5 - ϵ + (x[1] - x[2])) / (v[2] - v[1])
        i = 2
    elseif x[2] < 0 && v[1] == -1
        t = (0.5 - ϵ + x[1] + x[2]) / -(v[1] + v[2])
        i = 3
    elseif x[1] < 0 && v[2] == -1
        t = (0.5 - ϵ + x[1] + x[2]) / -(v[1] + v[2])
        i = 3
    elseif x[2] < 0 && v[1] == 1
        t = (0.5 - ϵ - (x[1] - x[2])) / (v[1] - v[2])
        i = 4
    elseif x[1] > 0 && v[2] == -1
        t = (0.5 - ϵ - (x[1] - x[2])) / (v[1] - v[2])
        i = 4
    end

    return t, i 
end


function compute_dt!(i, n, q, v, ϵ, dt)

    q_i = q[i]
    v_i = v[i]
    for j in 1:n
        if i != j
            tcoll = tempscoll(q[j], q_i, v[j], v_i, ϵ)
            if !isinf(tcoll)
                dt[i, j] = tcoll
            end
        end
    end

end

export ParticleCollisions

struct ParticleCollisions

    dt

    function ParticleCollisions(n, q, v, ϵ)

        dt = zeros(n, n)

        fill!(dt, Inf)

        for k in 1:n
            for l in (k+1):n
                dt[k, l] = tempscoll(q[l], q[k], v[l], v[k], ϵ)
            end
        end

        new( dt )

    end

end

function dt_min_position( collisions )
    p = argmin(collisions.dt)
    dt_min = collisions.dt[p]
    return dt_min, p[1], p[2]
end

export WallCollisions

struct WallCollisions

    dt 

    function WallCollisions( n, q, v, ϵ)
        dt = zeros(n, 4)

        fill!(dt, Inf)

        for k in 1:n
            t, m = tempscollmur(q[k], v[k], ϵ)
            !isinf(t) && (dt[k, m] = t)
        end

        new(dt)

    end

end

const wall_rebound = [
    v -> v[1] ==  1 ? - rot * v :   rot * v,
    v -> v[2] ==  1 ? - rot * v :   rot * v,
    v -> v[2] == -1 ?   rot * v : - rot * v,
    v -> v[2] == -1 ? - rot * v :   rot * v ]



function step!(n, ϵ, q, v, collisions, walls)

    tempsp, i1, i2 = dt_min_position(collisions)

    tempsm, j1, j2 = dt_min_position(walls) 

    if tempsp < tempsm

        for i in 1:n
            q[i] = q[i] + tempsp * v[i]
        end

        if v[i1]'v[i2] == 0

            v[[i1, i2]] = v[[i2, i1]]  # swap velocities

        elseif v[i1]'v[i2] == -1

            rot = [0 -1; 1 0]

            if (rot * v[i1])'*(q[i2] - q[i1]) < 0
                v[i1] = rot * v[i1]
                v[i2] = rot * v[i2]
            else
                v[i1] = -rot * v[i1]
                v[i2] = -rot * v[i2]
            end
        end

        collisions.dt .= collisions.dt .- tempsp
        reset!(collisions.dt, i1)
        reset!(collisions.dt, i2)

        walls.dt .= walls.dt .- tempsp
        walls.dt[i1,:] .= Inf
        walls.dt[i2,:] .= Inf

        compute_dt!(i1, n, q, v, ϵ, collisions.dt)
        compute_dt!(i2, n, q, v, ϵ, collisions.dt)

        t, i = tempscollmur(q[i1], v[i1], ϵ)
        if !isinf(t)
            walls.dt[i1, i] = t
        end

        t, i = tempscollmur(q[i2], v[i2], ϵ)
        if !isinf(t)
            walls.dt[i2, i] = t
        end

    else

        for i in 1:n
            q[i] = q[i] + tempsm * v[i]
        end

        v[j1] = wall_rebound[j2](v[j1])

        collisions.dt .-= tempsm
        reset!(collisions.dt, j1)

        walls.dt .-= tempsm
        walls.dt[j1,:] .= Inf

        compute_dt!(j1, n, q, v, ϵ, collisions.dt)

        t, i = tempscollmur(q[j1], v[j1], ϵ)
        if !isinf(t)
            walls.dt[j1, i] = t
        end

    end

    return min(tempsp, tempsm)

end

export box_particles

function box_particles(rng, n, ϵ)

    l1 = 0.5
    l2 = 0.5

    cost = 1 / sqrt(2)
    J4 = [cost cost; -cost cost]

    q = Vector{Float64}[]
    push!(q, [ 0,           0])
    push!(q, [ 0.0685359 ,  0.15134794])
    push!(q, [-0.37969185, -0.02715423])
    push!(q, [-0.23362135,  0.07949096])
    push!(q, [-0.03214348,  0.07088258])
    push!(q, [ 0.04120227,  0.33615893])
    push!(q, [-0.1507703 ,  0.32089688])
    push!(q, [-0.2       , -0.2])
    push!(q, [ 0.25      , -0.05])
    push!(q, [ 0.1       , -0.2])

    p = zeros(2)

    for klm in 11:n
        frein = 0
        overlap = 1
        while (overlap == 1) && frein < 10000

            p0 = [2 * l1 * rand() - l1, 2 * l2 * rand() - l2]
            p = J4 * p0
            p = p * (1 - 4 * ϵ) / sqrt(2)
            overlap2 = 0
            for subh in 1:(klm-1)
                if norm(p - q[subh], 1) < 2ϵ
                    overlap2 = 1
                end
            end

            overlap = overlap2
            frein += 1
        end

        if frein == 10000
            print("Echec tirage initial")
        end
        push!(q, p)
    end

    vitesses = [[1, 0], [0, 1], [-1, 0], [0, -1]]
    v = [vitesses[rand(1:end)] for i in 1:n]

    return q, v

end
