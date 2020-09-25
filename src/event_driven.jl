function tempscoll(x, y, v, w, ϵ)

    deltat1 = Inf
    deltat2 = Inf

    if (y - x)'v > 0 > (y - x)'w

        dx = y[1] - x[1]
        dy = y[2] - x[2]
        adx = abs(dx)
        ady = abs(dy)
        adxy = abs(adx - ady)

        if w[1] == 1 && v[1] == -1 && ady < 2ϵ
            deltat1 = (2ϵ - dx - ady) / 2
            deltat2 = (- 2ϵ - dx + ady) / 2
        elseif v[1] == 1 && w[1] == -1 && ady < 2ϵ
            deltat1 = (2ϵ + dx - ady) / 2
            deltat2 = (- 2ϵ + dx + ady) / 2
        elseif w[2] == 1 && v[2] == -1 && adx < 2ϵ
            deltat1 = (2ϵ - adx - dy) / 2
            deltat2 = (- 2ϵ + adx - dy) / 2
        elseif v[2] == 1 && w[2] == -1 && adx < 2ϵ
            deltat1 = (2ϵ - adx + dy) / 2
            deltat2 = (- 2ϵ + adx + dy) / 2
        elseif w[1] == 1 && v[2] == 1 && adxy < 2ϵ
            deltat1 = (dy - dx - 2ϵ) / 2
            deltat2 = (dy - dx + 2ϵ) / 2
        elseif v[1] == 1 && w[2] == 1 && adxy < 2ϵ
            deltat1 = (-dy + dx - 2ϵ) / 2
            deltat2 = (-dy + dx + 2ϵ) / 2
        elseif w[1] == 1 && v[2] == -1 && adxy < 2ϵ
            deltat1 = (- dx - dy - 2ϵ) / 2
            deltat2 = (- dx - dy + 2ϵ) / 2
        elseif v[1] == 1 && w[2] == -1 && adxy < 2ϵ
            deltat1 = (dx + dy - 2ϵ) / 2
            deltat2 = (dx + dy + 2ϵ) / 2
        elseif w[1] == -1 && v[2] == 1 && adxy < 2ϵ
            deltat1 = (dx + dy + 2ϵ) / 2
            deltat2 = (dx + dy - 2ϵ) / 2
        elseif v[1] == -1 && w[2] == 1 && adxy < 2ϵ
            deltat1 = (- dx - dy + 2ϵ) / 2
            deltat2 = (- dx - dy - 2ϵ) / 2
        elseif w[1] == -1 && v[2] == -1 && adxy < 2ϵ
            deltat1 = (dx - dy + 2ϵ) / 2
            deltat2 = (dx - dy - 2ϵ) / 2
        elseif v[1] == -1 && w[2] == -1 && adxy < 2ϵ
            deltat1 = (- dx + dy + 2ϵ) / 2
            deltat2 = (- dx + dy - 2ϵ) / 2
        end

    end

    return min(deltat1, deltat2)

end


function tempscollmur(x, v, ϵ)
    t = Inf
    i = 0

    if x[2] > 0 && v[1] == 1
        t = (0.5 - (x[1] + x[2]) - ϵ) / (v[1] + v[2])
        i = 1
    elseif x[0] > 0 && v[2] == 1
        t = (0.5 - (x[1] + x[2]) - ϵ) / (v[1] + v[2])
        i = 1
    elseif x[2] > 0 && v[1] == -1
        t = (0.5 - ϵ + (x[1] - x[2])) / (v[2] - v[1])
        i = 2
    elseif x[1] < 0 && v[2] == 1
        t = (0.5 - ϵ + (x[1] - x[2])) / (v[2] - v[1])
        i = 2
    elseif x[2] < 0 && v[0] == -1
        t = (0.5 - ϵ + x[0] + x[2]) / -(v[1] + v[2])
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

    return t, i - 1
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


struct ParticleCollision

    dt

    function ParticleCollision(n, q, v, ϵ)

        dt = zeros(n, n)

        fill!(dt, Inf)
        self.dt_min = Inf

        for k in 1:n
            for l in (k+1):n
                dt[k, l] = tempscoll(q[l], q[k], v[l], v[k], ϵ)
            end
        end

        new( dt )

    end

end

function dt_min( collisions )
    p = argmin(collisions.dt)
    dt_min = self.dt[p]
    return dt_min, p[1], p[2]
end


struct WallCollision

    dt 

    function WallCollision( n, q, v, ϵ)
        dt = zeros(n, 4)

        fill!(dt, Inf)

        for k in 1:n
            t, m = tempscollmur(q[k], v[k], ϵ)
            dt[k, m] = t
        end

        new(dt)

    end

end


const wall_rebound = [
    v -> v[1] ==  1 ? - rot * v :   rot * v,
    v -> v[2] ==  1 ? - rot * v :   rot * v,
    v -> v[2] == -1 ?   rot * v : - rot * v,
    v -> v[2] == -1 ? - rot * v :   rot * v ]


function init_particles(rng, n, ϵ,)

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

    for klm in 10:n
        frein = 0
        overlap = 1
        while (overlap == 1) && frein < 10000

            p0 = array([[2 * l1 * rand() - l1, 2 * l2 * rand() - l2]])
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

function step!(n, ϵ, q, v, collisions, walls)

    tempsp, i1, i2 = dt_min_position(collisions)

    tempsm, j1, j2 = dt_min_position(walls) 

    if tempsp < tempsm

        q += tempsp * v

        if v[i1]'v[i2] == 0

            v[[i1, i2]] = v[[i2, i1]]  # swap velocities

        elseif v[i1]'v[i2] == -1

            rot = [0 -1; 1 0]

            if (rot * v[i1])'*(q[i2] - q[i1]) < 0
                v[i1] = rot * v[i1]
                v[i2] = rot * v[i2]
            elseif (rot * v[i1])'*(q[i2] - q[i1]) > 0
                self.v[i1] = -rot * v[i1]
                self.v[i2] = -rot * v[i2]
            end
        end

        collisions.dt = collisions.dt - tempsp
        reset(collisions, i1)
        reset(collisions, i2)

        murs.dt = murs.dt - tempsp
        reset!(murs, i1)
        reset!(murs, i2)

        compute_dt(i1, n, q, v, ϵ, collisions.dt)
        compute_dt(i2, n, q, v, ϵ, collisions.dt)

        t, i = tempscollmur(q[i1], v[i1], ϵ)
        if !isinf(t)
            self.murs.dt[i1, i] = t
        end

        t, i = tempscollmur(q[i2], v[i2], ϵ)
        if !isinf(t)
            self.murs.dt[i2, i] = t
        end

    else

        q += tempsm * self.v

        v[j1] = wall_rebound[j2](v[j1])

        collisions.dt .-= tempsm
        reset!(collisions, j1)

        murs.dt .-= tempsm
        reset!(murs, j1)

        compute_dt(j1, n, q, v, ϵ, collisions.dt)

        t, i = tempscollmur(q[j1], v[j1], ϵ)
        if !isinf(t)
            murs.dt[j1, i] = t
        end

    end

    return  min(tempsp, tempsm)

end
