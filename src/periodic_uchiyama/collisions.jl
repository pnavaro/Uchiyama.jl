function tempscoll(x, y, v, w, epsilon)

    deltat1 = Inf
    deltat2 = Inf

    if (y - x) *  v > 0 > (y - x) * w

        dx = y[1] - x[1]
        dy = y[2] - x[2]
        adx = abs(dx)
        ady = abs(dy)
        adxy = abs(adx - ady)
        depsilon = 2.0 * epsilon

        if w[1] == 1 && v[1] == -1 && ady < depsilon
            deltat1 = (depsilon - dx - ady) / 2
            deltat2 = (- depsilon - dx + ady) / 2
        elseif v[1] == 1 && w[1] == -1 && ady < depsilon
            deltat1 = (depsilon + dx - ady) / 2
            deltat2 = (- depsilon + dx + ady) / 2
        elseif w[2] == 1 && v[1] == -1 && adx < depsilon
            deltat1 = (depsilon - adx - dy) / 2
            deltat2 = (- depsilon + adx - dy) / 2
        elseif v[2] == 1 && w[2] == -1 && adx < depsilon
            deltat1 = (depsilon - adx + dy) / 2
            deltat2 = (- depsilon + adx + dy) / 2
        elseif w[1] == 1 && v[2] == 1 && adxy < depsilon
            deltat1 = (dy - dx - depsilon) / 2
            deltat2 = (dy - dx + depsilon) / 2
        elseif v[1] == 1 && w[2] == 1 && adxy < depsilon
            deltat1 = (-dy + dx - depsilon) / 2
            deltat2 = (-dy + dx + depsilon) / 2
        elseif w[1] == 1 && v[2] == -1 && adxy < depsilon
            deltat1 = (- dx - dy - depsilon) / 2
            deltat2 = (- dx - dy + depsilon) / 2
        elseif v[1] == 1 && w[2] == -1 && adxy < depsilon
            deltat1 = (dx + dy - depsilon) / 2
            deltat2 = (dx + dy + depsilon) / 2
        elseif w[1] == -1 && v[2] == 1 && adxy < depsilon
            deltat1 = (dx + dy + depsilon) / 2
            deltat2 = (dx + dy - depsilon) / 2
        elseif v[1] == -1 && w[2] == 1 && adxy < depsilon
            deltat1 = (- dx - dy + depsilon) / 2
            deltat2 = (- dx - dy - depsilon) / 2
        elseif w[1] == -1 && v[2] == -1 && adxy < depsilon
            deltat1 = (dx - dy + depsilon) / 2
            deltat2 = (dx - dy - depsilon) / 2
        elseif v[1] == -1 && w[2] == -1 && adxy < depsilon
            deltat1 = (- dx + dy + depsilon) / 2
            deltat2 = (- dx + dy - depsilon) / 2
        end

    end

    return min(deltat1, deltat2)

end

offset = [[0, 0], [1, 0], [-1, 0], [0, 1], [0, -1], [-1, 1], [1, 1], [1, -1], [-1, -1]]

function compute_dt(i, q, v, epsilon, dt, fantome):

    n = length(q)

    for j in 1:n
        if i != j
            dt_local = Inf
            k = 1

            while (isinf(dt_local) && k < 10)
                dt_local = tempscoll(q[j] + offset[k], q[i], v[j], v[i], epsilon)
                k += 1
            end

            dt[i,j] = dt_local
            fantome[i,j] = k
        end
    end

end


struct PCollisionMatrix

    dt :: Array{Float64,2}
    dt_min :: Float64
    fantome :: Array{Int, 2}

    function PCollision(npart, q, v, epsilon)

        dt = zeros(Float64, (npart, npart))
        fantome = zeros(Int, (npart, npart))
        fill!(dt, Inf)
        fill!(fantome, 0)
        dt_min = Inf

        for k in 1:npart
            for l in (k+1):npart
                dt_local = Inf
                i = 1

                while (isinf(dt_local) and i < 10)
                    dt_local = tempscoll(q[l] + offset[i], q[k], v[l], v[k], epsilon)
                    i += 1
                end

                dt[k, l] = dt_local
                fantome[k, l] = i
            end
        end

        new( dt, dt_min, fantome )

    end

end

function dt_min_position(self)
    i, j = unravel_index(argmin(self.dt), (self.npart, self.npart))
    self.dt_min = self.dt[i, j]
    self.fant = self.fantome[i, j]
    return i, j
end

function reset(self, i)
    self.dt[i, :] .= Inf
    self.dt[:, i] .= Inf
end


