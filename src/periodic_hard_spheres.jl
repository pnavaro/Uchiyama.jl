function step!(hs :: HardSpheres, pc :: PeriodicCollisions)

    n = hs.n
    ϵ = hs.ϵ

    dtmin, fant, i1, i2 = dt_min_position(pc)

    for i in 1:n
        qnew = hs.q[i] .+ dtmin .* hs.v[i]
        hs.q[i] = mod.(qnew, 1.)
    end

    δv = hs.v[i2] .- hs.v[i1]
    δx = hs.q[i2] .+ offset[fant] .- hs.q[i1]

    J = (δv'δx) ./ 2ϵ

    hs.v[i1] = hs.v[i1] .+ J .* δx ./ 2ϵ
    hs.v[i2] = hs.v[i2] .- J .* δx ./ 2ϵ

    pc.dt .-= dtmin
    reset!( pc, i1)
    reset!( pc, i2)

    compute_dt!(pc, i1, hs)
    compute_dt!(pc, i2, hs)

    return dtmin

end
