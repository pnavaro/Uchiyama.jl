const wall_rebound = [
    v -> v[1] ==  1 ? - rot * v :   rot * v,
    v -> v[2] ==  1 ? - rot * v :   rot * v,
    v -> v[2] == -1 ?   rot * v : - rot * v,
    v -> v[2] == -1 ? - rot * v :   rot * v ]



function step!(sq :: SquareParticles, pc, bc)

    n = sq.n
    ϵ = sq.ϵ
    tempsp, i1, i2 = dt_min_position(pc)
    tempsm, j1, j2 = dt_min_position(bc) 

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

        pc.dt .= pc.dt .- tempsp
        pc.dt[:, i1] .= Inf
        pc.dt[:, i2] .= Inf

        bc.dt .= bc.dt .- tempsp
        bc.dt[i1,:] .= Inf
        bc.dt[i2,:] .= Inf

        compute_dt!(pc, i1, sq)
        compute_dt!(pc, i2, sq)

        t, i = compute_dt(sq, i1)
        bc.dt[i1, i] = t

        t, i = compute_dt(sq, i2)
        bc.dt[i2, i] = t

    else

        for i in 1:n
            sq.q[i] = sq.q[i] + tempsm * sq.v[i]
        end

        sq.v[j1] = wall_rebound[j2](sq.v[j1])

        pc.dt .-= tempsm
        reset!(pc, j1)

        bc.dt .-= tempsm
        bc.dt[j1,:] .= Inf

        compute_dt!(pc, j1, sq)

        t, i = compute_dt(sq, j1)
        bc.dt[j1, i] = t

    end

    return min(tempsp, tempsm)

end

