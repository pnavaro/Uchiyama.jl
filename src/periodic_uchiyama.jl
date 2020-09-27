using Random
using LinearAlgebra

function tempscoll(x, y, v, w, ϵ)

    deltat1 = Inf
    deltat2 = Inf

    if (y - x)'v > 0.0 > (y - x)'w

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

function compute_dt!(i, squares, dt, fantome)

    for j in 1:squares.n
        if i != j
            dt_local = Inf
            k = 0
            while (isinf(dt_local) && k < 9)
                k += 1
                dt_local = tempscoll(squares.q[j] .+ offset[k], 
                                     squares.q[i], 
                                     squares.v[j], 
                                     squares.v[i], 
                                     squares.ϵ)
            end

            dt[i,j] = dt_local
            fantome[i,j] = k
        end
    end

end


export PeriodicCollisions

struct PeriodicCollisions

    dt :: Array{Float64,2}
    fantome :: Array{Int, 2}

    function PeriodicCollisions( squares :: Squares)

        #aij is stored in AP(i+j(j-1)/2) for $i \leq j$;

        dt = zeros(Float64, (squares.n, squares.n))
        fantome = zeros(Int, (squares.n, squares.n))
        fill!(dt, Inf)
        fill!(fantome, 0)

        for i in 1:squares.n
            for j in (i+1):squares.n
                dt_local = Inf
                k = 0
                while (isinf(dt_local) && k < 9)
                    k += 1
                    dt_local = tempscoll(squares.q[j] + offset[k], 
                                         squares.q[i], 
                                         squares.v[j], 
                                         squares.v[i], 
                                         squares.ϵ)
                end

                dt[i, j] = dt_local
                fantome[i, j] = k

            end
        end

        new( dt, fantome )

    end

end

function dt_min_position(self :: PeriodicCollisions)
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

export step!

function step!(squares, collisions)

    """ Compute smaller time step between two collisions """

    dt, num_fant, i1, i2 = dt_min_position(collisions)

    for i in 1:squares.n
        qnew = squares.q[i] .+ dt .* squares.v[i]
        squares.q[i] = mod.(qnew, 1.)
    end

    """ Collide the two particles i1, i2 """

    if squares.v[i1]'squares.v[i2] == 0

        squares.v[[i1,i2]] = squares.v[[i2,i1]] # swap velocities

    elseif squares.v[i1]'squares.v[i2] == -1

        if (rot * squares.v[i1])'*(squares.q[i2] .+ offset[num_fant] .- squares.q[i1]) < 0
            squares.v[i1] = rot * squares.v[i1]
            squares.v[i2] = rot * squares.v[i2]
        elseif (rot * squares.v[i1])'*(squares.q[i2] .+ offset[num_fant] .- squares.q[i1]) > 0
            squares.v[i1] = - rot * squares.v[i1]
            squares.v[i2] = - rot * squares.v[i2]
        else # very rare case
            squares.v[[i1,i2]] = squares.v[[i2,i1]] # swap velocities
        end

    end

    collisions.dt .-= dt
    collisions.dt[:, i1] .= Inf
    collisions.dt[:, i2] .= Inf

    compute_dt!(i1, squares, collisions.dt, collisions.fantome)
    compute_dt!(i2, squares, collisions.dt, collisions.fantome)

    return dt

end

