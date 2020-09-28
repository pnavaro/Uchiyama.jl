function overlap( p, q, ϵ )

   for r in q
       all(abs.(p .- r) .< 2ϵ) && (return true)
   end

   return false

end

function box_particles(rng, n, ϵ)

    l1 = 0.5
    l2 = 0.5

    cost = 1 / sqrt(2)
    J4 = [cost cost; -cost cost]

    q = Vector{Float64}[]
    push!(q, [ 0,           0])
    push!(q, [ 0.0685359 ,  0.15134794])
    push!(q, [-0.37969185, -0.02715423])
    push!(q, [-0.23362135,  0.07949096])
    push!(q, [-0.03214348,  0.07088258])
    push!(q, [ 0.04120227,  0.33615893])
    push!(q, [-0.1507703 ,  0.32089688])
    push!(q, [-0.2       , -0.2])
    push!(q, [ 0.25      , -0.05])
    push!(q, [ 0.1       , -0.2])

    p = zeros(2)

    for klm in 11:n
        frein = 0
        overlap = 1
        while (overlap == 1) && frein < 10000

            p0 = [2 * l1 * rand(rng) - l1, 2 * l2 * rand(rng) - l2]
            p = J4 * p0
            p = p * (1 - 4 * ϵ) / sqrt(2)
            overlap2 = 0
            for subh in 1:(klm-1)
                if norm(p - q[subh], 1) < 2ϵ
                    overlap2 = 1
                end
            end

            overlap = overlap2
            frein += 1
        end

        if frein == 10000
            print("Echec tirage initial")
        end
        push!(q, p)
    end

    vitesses = [[1, 0], [0, 1], [-1, 0], [0, -1]]
    v = [vitesses[rand(rng,1:end)] for i in 1:n]

    return q, v

end

export SquareParticles

struct SquareParticles <: Particles

    n
    q
    v
    ϵ

    function SquareParticles(q :: Vector{Vector{Float64}}, v, ϵ :: AbstractFloat)

        n = length(q)

        @assert n == length(v)

        new(n, q, v, ϵ)

    end


    function SquareParticles(rng, n, ϵ; option = :none )

        if option == :box

           q, v = box_particles(rng, n, ϵ)
        
        else

           q = Vector{Float64}[]

           push!(q, [0.5,0.5])

           for i in 1:n-1

               p = [0.5, 0.5]
               k = 0
               while overlap(p, q, ϵ) && k < 1000
                   rand!(rng, p)
                   k += 1
               end

               if k < 1000
                  push!(q, p)
               else
                  @error "echec du tirage $i"
               end

           end

           vitesses = [[1, 0], [0, 1], [-1, 0], [0, -1]]
           v = [vitesses[rand(rng, 1:end)] for i in 1:n]
        end 

        new(n, q, v, ϵ)

    end

end


function compute_dt(p :: SquareParticles, i, j, k = 1)

    ϵ = p.ϵ
    deltat1 = Inf
    deltat2 = Inf
    x = p.q[i] .+ offset[k]
    y = p.q[j]
    v = p.v[i]
    w = p.v[j]

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

function compute_dt(p :: SquareParticles, j)

    ϵ = p.ϵ
    t = Inf
    x = p.q[j]
    v = p.v[j]
    i = 1

    if x[2] > 0 && v[1] == 1
        t = (0.5 - (x[1] + x[2]) - ϵ) / (v[1] + v[2])
        i = 1
    elseif x[1] > 0 && v[2] == 1
        t = (0.5 - (x[1] + x[2]) - ϵ) / (v[1] + v[2])
        i = 1
    elseif x[2] > 0 && v[1] == -1
        t = (0.5 - ϵ + (x[1] - x[2])) / (v[2] - v[1])
        i = 2
    elseif x[1] < 0 && v[2] == 1
        t = (0.5 - ϵ + (x[1] - x[2])) / (v[2] - v[1])
        i = 2
    elseif x[2] < 0 && v[1] == -1
        t = (0.5 - ϵ + x[1] + x[2]) / -(v[1] + v[2])
        i = 3
    elseif x[1] < 0 && v[2] == -1
        t = (0.5 - ϵ + x[1] + x[2]) / -(v[1] + v[2])
        i = 3
    elseif x[2] < 0 && v[1] == 1
        t = (0.5 - ϵ - (x[1] - x[2])) / (v[1] - v[2])
        i = 4
    elseif x[1] > 0 && v[2] == -1
        t = (0.5 - ϵ - (x[1] - x[2])) / (v[1] - v[2])
        i = 4
    end

    return t, i 

end

