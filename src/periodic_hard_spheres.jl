function hscompute_dt(i, n, q, v, ϵ, dt, fantome)

    for j in 1:n
        if i != j
            dt_local = Inf
            k = 0

            while (isinf(dt_local) && k < 9)
                k += 1
                dt_local = hstempscoll(q[j] + offset[k], q[i], v[j], v[i], ϵ)
            end

            if !isinf(dt_local)
                dt[i,j] = dt_local
                fantome[i,j] = k
            end
        end
    end
end

export PCollisionMatrix
            
struct PCollisionMatrix

    dt
    fantome

    function PCollisionMatrix(n, q, v, ϵ)
        dt = zeros(n, n)
        fantome = zeros(Int8,(n, n))
        fill!(dt, Inf)
        fill!(fantome, 0)

        for k in 1:n
            for l in (k+1):n
                dt_local = Inf
                i = 0
                while (isinf(dt_local) && i < 9)
                    i += 1
                    dt_local = hstempscoll(q[l] + offset[i], q[k], v[l], v[k], ϵ)
                end

                if !isinf(dt_local)
                    dt[k, l] = dt_local
                    fantome[k, l] = i
                end
            end
        end

        new(dt, fantome)
    end

end

function dt_min_position(self :: PCollisionMatrix)
    p = argmin(self.dt)
    dt_min = self.dt[p]
    num_fant = self.fantome[p]
    return dt_min, num_fant, p[1], p[2]
end

function reset!(self :: PCollisionMatrix, i)
    self.dt[i, :] .= Inf
    self.dt[:, i] .= Inf
    self.fantome[:, i] .= 0
    self.fantome[i, :] .= 0
end

function step!(n, ϵ, q, v, collisions :: PCollisionMatrix)

    dt, num_fant, i1, i2 = dt_min_position(collisions)

    for i in 1:n
        qnew = q[i] .+ dt .* v[i]
        q[i] = mod.(qnew, 1.)
    end

    J = ((v[i2] - v[i1])'*(q[i2] .+ offset[num_fant] .- q[i1])) / 2ϵ

    v[i1] = v[i1] + J * (q[i2] .+ offset[num_fant] .- q[i1]) / 2ϵ
    v[i2] = v[i2] - J * (q[i2] .+ offset[num_fant] .- q[i1]) / 2ϵ

    collisions.dt .-= dt
    reset!( collisions, i1)
    reset!( collisions, i2)

    hscompute_dt(i1, n, q, v, ϵ, collisions.dt, collisions.fantome)
    hscompute_dt(i2, n, q, v, ϵ, collisions.dt, collisions.fantome)

    return dt

end
