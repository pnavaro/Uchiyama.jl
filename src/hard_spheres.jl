const hs_wall_rebound = [
    v -> [-v[1],  v[2]],
    v -> [ v[1], -v[2]]
]



function step!(hs :: HardSpheres, pc::ParticleCollisions, bc::BoxCollisions)

    ϵ = hs.ϵ
    n = hs.n

    tempsp, i1, i2 = dt_min_position(pc)
    tempsm, j1, j2 = dt_min_position(bc)

    if tempsp < tempsm

        for i in 1:n
            hs.q[i] = hs.q[i] + tempsp * hs.v[i]
        end

        J = ((hs.v[i2] - hs.v[i1])'*(hs.q[i2] - hs.q[i1])) / 2ϵ
        hs.v[i1] = hs.v[i1] + J * (hs.q[i2] - hs.q[i1]) / 2ϵ
        hs.v[i2] = hs.v[i2] - J * (hs.q[i2] - hs.q[i1]) / 2ϵ

        pc.dt .-= tempsp
        reset!(pc.dt, i1)
        reset!(pc.dt, i2)

        bc.dt .-= tempsp
        bc.dt[i1,:] .= Inf
        bc.dt[i2,:] .= Inf

        compute_dt!(pc.dt, i1, hs)
        compute_dt!(pc.dt, i2, hs)

        t, i = hs(hs.q[i1], hs.v[i1])
        !isinf(t) && ( bc.dt[i1, i] = t)

        t, i = hs(hs.q[i2], hs.v[i2])
        !isinf(t) && ( bc.dt[i2, i] = t)

    else

        for i in 1:n
            hs.q[i] = hs.q[i] + tempsm * hs.v[i]
        end

        hs.v[j1] = hs_wall_rebound[j2](hs.v[j1])

        pc.dt .-= tempsm
        reset!(pc.dt, j1)

        bc.dt .-= tempsm
        bc.dt[j1,:] .= Inf

        compute_dt!(pc.dt, j1, hs)

        t, i = hs(hs.q[j1], hs.v[j1])
        !isinf(t) && (bc.dt[j1, i] = t)

    end

    return min(tempsp, tempsm)

end
