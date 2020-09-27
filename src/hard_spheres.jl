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
            hs.q[i] = hs.q[i] .+ tempsp .* hs.v[i]
        end

        δv = hs.v[i2] - hs.v[i1]
        δx = hs.q[i2] - hs.q[i1]

        J = (δv'δx) ./ 2ϵ

        hs.v[i1] = hs.v[i1] .+ J .* δx ./ 2ϵ
        hs.v[i2] = hs.v[i2] .- J .* δx ./ 2ϵ

        pc.dt .-= tempsp
        reset!(pc, i1)
        reset!(pc, i2)

        bc.dt .-= tempsp
        reset!(bc, i1)
        reset!(bc, i2)

        compute_dt!(pc, i1, hs)
        compute_dt!(pc, i2, hs)

        t, i = compute_dt(hs, i1)
        bc.dt[i1, i] = t

        t, i = compute_dt(hs, i2)
        bc.dt[i2, i] = t

    else

        for i in 1:n
            hs.q[i] = hs.q[i] + tempsm * hs.v[i]
        end

        hs.v[j1] = hs_wall_rebound[j2](hs.v[j1])

        pc.dt .-= tempsm
        reset!(pc, j1)
        compute_dt!(pc, j1, hs)

        bc.dt .-= tempsm
        reset!(bc, j1)
        t, i = compute_dt(hs, j1)
        bc.dt[j1, i] = t

    end

    return min(tempsp, tempsm)

end
