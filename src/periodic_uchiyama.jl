const rot = [0 -1; 1 0]

export step!

function step!(sq, pc)

    dtmin, fant, i1, i2 = dt_min_position(pc)

    for i in 1:sq.n
        sq.q[i] = mod.(sq.q[i] .+ dtmin .* sq.v[i], 1)
    end

    if sq.v[i1]'sq.v[i2] == 0

        sq.v[[i1,i2]] = sq.v[[i2,i1]] # swap velocities

    elseif sq.v[i1]'sq.v[i2] == -1

        if (rot * sq.v[i1])'*(sq.q[i2] .+ offset[fant] .- sq.q[i1]) < 0
            sq.v[i1] = rot * sq.v[i1]
            sq.v[i2] = rot * sq.v[i2]
        elseif (rot * sq.v[i1])'*(sq.q[i2] .+ offset[fant] .- sq.q[i1]) > 0
            sq.v[i1] = - rot * sq.v[i1]
            sq.v[i2] = - rot * sq.v[i2]
        else # very rare case
            sq.v[[i1,i2]] = sq.v[[i2,i1]] # swap velocities
        end

    end

    pc.dt .-= dtmin
    reset!(pc, i1)
    reset!(pc, i2)

    compute_dt!(pc, i1, sq)
    compute_dt!(pc, i2, sq)

    return dtmin

end

