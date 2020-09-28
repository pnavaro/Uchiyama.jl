export PeriodicCollisions

struct PeriodicCollisions

    dt :: Array{Float64,2}
    fantome :: Array{Int, 2}

    function PeriodicCollisions( p :: Particles)

        #aij is stored in AP(i+j(j-1)/2) for $i \leq j$;

        dt = zeros(Float64, (p.n, p.n))
        fantome = zeros(Int, (p.n, p.n))
        fill!(dt, Inf)

        for i in 1:p.n
            for j in 1:p.n
                if i != j 
                   dt_local = Inf
                   k = 0
                   while (isinf(dt_local) && k < 9)
                       k += 1
                       dt_local = compute_dt(p, j, i, k)
                   end

                   dt[i, j] = dt_local
                   fantome[i, j] = k
                end
            end
        end

        new( dt, fantome )

    end

end

function compute_dt!(pc :: PeriodicCollisions, i, p :: Particles)

    for j in 1:p.n
        if i != j
            dt_local = Inf
            k = 0
            while (isinf(dt_local) && k < 9)
                k += 1
                dt_local = compute_dt(p, j, i, k)
            end

            pc.dt[i,j] = dt_local
            pc.fantome[i,j] = k
        end
    end

end

function dt_min_position(pc :: PeriodicCollisions)
    p = argmin(pc.dt)
    dt_min = pc.dt[p]
    num_fant = pc.fantome[p]
    return dt_min, num_fant, p[1], p[2]
end

function reset!(pc :: PeriodicCollisions, i)
    pc.dt[i, :] .= Inf
    pc.dt[:, i] .= Inf
end
