function step!(hs :: HardSpheres, pc :: PeriodicCollisions)

    n = hs.n
    ϵ = hs.ϵ

    dtmin, fant, i1, i2 = dt_min_position(pc)

    if fant > 0
       qa = [hs.q[i2] + offset[k] + dtmin*hs.v[i2] for k in 1:9]
       qb = hs.q[i1] + dtmin .* hs.v[i1]
       vr = argmin([norm(qx-qb) for qx in qa])
    
       hs.q[i2] = qa[vr]
       hs.q[i1] = qb
    
       J = (dot(hs.v[i2]-hs.v[i1],hs.q[i2]-hs.q[i1]))/ 2ϵ
       hs.v[i1] = hs.v[i1] + J * (hs.q[i2]-hs.q[i1]) / 2ϵ
       hs.v[i2] = hs.v[i2] - J * (hs.q[i2]-hs.q[i1]) / 2ϵ
       hs.q[i1] = mod.(hs.q[i1],1)
       hs.q[i2] = mod.(hs.q[i2],1)

    end


    for i in [1:min(i1,i2)-1;min(i1,i2)+1:max(i1,i2)-1;max(i1,i2)+1:n]
        qnew = hs.q[i] .+ dtmin .* hs.v[i]
        hs.q[i] = mod.(qnew, 1.)
    end

    pc.dt .-= dtmin
    reset!( pc, i1)
    reset!( pc, i2)

    compute_dt!(pc, i1, hs)
    compute_dt!(pc, i2, hs)

    return dtmin

end
