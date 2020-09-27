const wall_rebound = [
    v -> v[1] ==  1 ? - rot * v :   rot * v,
    v -> v[2] ==  1 ? - rot * v :   rot * v,
    v -> v[2] == -1 ?   rot * v : - rot * v,
    v -> v[2] == -1 ? - rot * v :   rot * v ]



function step!(sq :: Squares, collisions, walls)

    n = sq.n
    ϵ = sq.ϵ
    tempsp, i1, i2 = dt_min_position(collisions)

    tempsm, j1, j2 = dt_min_position(walls) 

    if tempsp < tempsm

        for i in 1:n
            sq.q[i] = sq.q[i] + tempsp * sq.v[i]
        end

        if sq.v[i1]'sq.v[i2] == 0

            sq.v[[i1, i2]] = sq.v[[i2, i1]]  # swap velocities

        elseif sq.v[i1]'sq.v[i2] == -1

            rot = [0 -1; 1 0]

            if (rot * sq.v[i1])'*(sq.q[i2] - sq.q[i1]) < 0
                sq.v[i1] = rot * sq.v[i1]
                sq.v[i2] = rot * sq.v[i2]
            else
                sq.v[i1] = -rot * sq.v[i1]
                sq.v[i2] = -rot * sq.v[i2]
            end
        end

        collisions.dt .= collisions.dt .- tempsp
        collisions.dt[:, i1] .= Inf
        collisions.dt[:, i2] .= Inf

        walls.dt .= walls.dt .- tempsp
        walls.dt[i1,:] .= Inf
        walls.dt[i2,:] .= Inf

        compute_dt!(collisions.dt, i1, sq)
        compute_dt!(collisions.dt, i2, sq)

        t, i = sq(sq.q[i1], sq.v[i1])
        !isinf(t) && (walls.dt[i1, i] = t)

        t, i = sq(sq.q[i2], sq.v[i2])
        !isinf(t) && (walls.dt[i2, i] = t)

    else

        for i in 1:n
            sq.q[i] = sq.q[i] + tempsm * sq.v[i]
        end

        sq.v[j1] = wall_rebound[j2](sq.v[j1])

        collisions.dt .-= tempsm
        reset!(collisions.dt, j1)

        walls.dt .-= tempsm
        walls.dt[j1,:] .= Inf

        compute_dt!(collisions.dt, j1, sq)

        t, i = sq(sq.q[j1], sq.v[j1])
        !isinf(t) && (walls.dt[j1, i] = t)

    end

    return min(tempsp, tempsm)

end

