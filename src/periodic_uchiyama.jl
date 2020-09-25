using LinearAlgebra

function overlap( p, q, k, ϵ )
   for j in 1:(k-1)
       (norm(p .- q[j], 1) < 2ϵ) && (return true)
   end
   return false
end

function init_particles(n, ϵ, rng)

    q = [zeros(2) for i in 1:n]
    q[1] = [0.5,0.5]

    
    for k in 2:n

        p = copy(q[1])

        while overlap(p, q, k, ϵ)
            p .= rand(2)
        end

        q[k] .= p

    end

    vitesses = [[1, 0], [0, 1], [-1, 0], [0, -1]]
    v = [vitesses[rand(1:end)] for i in 1:n]

    return q, v

end


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


const offset = [[0, 0], [1,  0], [-1,  0], 
                [0, 1], [0, -1], [-1,  1], 
                [1, 1], [1, -1], [-1, -1]]

function compute_dt!(i, n, q, v, ϵ, dt, fantome)

    for j in 1:n
        if i != j
            dt_local = Inf
            k = 0
            while (isinf(dt_local) && k < 9)
                k += 1
                dt_local = tempscoll(q[j] .+ offset[k], q[i], v[j], v[i], ϵ)
            end

            dt[i,j] = dt_local
            fantome[i,j] = k
        end
    end

end


struct PCollisionMatrix

    dt :: Array{Float64,2}
    fantome :: Array{Int, 2}

    function PCollisionMatrix(n, q, v, ϵ)

        dt = zeros(Float64, (n, n))
        fantome = zeros(Int, (n, n))
        fill!(dt, Inf)
        fill!(fantome, 0)

        for i in 1:n
            compute_dt!(i, i, q, v, ϵ, dt, fantome)
        end


        new( dt, fantome )

    end

end

function dt_min_position(self)
    p = argmin(self.dt)
    dt_min = self.dt[p]
    num_fant = self.fantome[p]
    return dt_min, num_fant, p[1], p[2]
end

function reset!(dt, i)
    dt[i, :] .= Inf
    dt[:, i] .= Inf
end

const rot = [0 -1; 1 0]

function step!(n, ϵ, q, v, collisions)

    """ Compute smaller time step between two collisions """

    dt, num_fant, i1, i2 = dt_min_position(collisions)

    for i in 1:n
        qnew = q[i] .+ dt .* v[i]
        qnew = qnew .% [1.,1.]
        q[i] .= qnew
    end

    """ Collide the two particles i1, i2 """

    if v[i1]'v[i2] == 0

        tmp = v[i1] 
        v[i1] .= copy(v[i2])
        v[i2] .= tmp  # swap velocities

    elseif v[i1]'v[i2] == -1

        if (rot * v[i1])'*(q[i2] .+ offset[num_fant] .- q[i1]) < 0
            v[i1] .= rot * v[i1]
            v[i2] .= rot * v[i2]
        elseif (rot * v[i1])'*(q[i2] .+ offset[num_fant] .- q[i1]) > 0
            v[i1] .= - rot * v[i1]
            v[i2] .= - rot * v[i2]
        end

    end

    collisions.dt .-= dt
    reset!(collisions.dt, i1)
    reset!(collisions.dt, i2)

    compute_dt!(i1, n, q, v, ϵ, collisions.dt, collisions.fantome)
    compute_dt!(i2, n, q, v, ϵ, collisions.dt, collisions.fantome)

    return dt

end

