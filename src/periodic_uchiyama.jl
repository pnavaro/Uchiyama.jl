using Random
using LinearAlgebra

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
                dt_local = squares(squares.q[j] .+ offset[k], 
                                   squares.q[i], 
                                   squares.v[j], 
                                   squares.v[i])
            end

            dt[i,j] = dt_local
            fantome[i,j] = k
        end
    end

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

