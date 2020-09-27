abstract type Particles end

function overlap( p, q, ϵ )

   for r in q
       all(abs.(p .- r) .< 2ϵ) && (return true)
   end

   return false

end

export Squares

struct Squares <: Particles

    n
    q
    v
    ϵ

    function Squares(rng, n, ϵ )

        q = Vector{Float64}[]

        push!(q, [0.5,0.5])

        for i in 1:n-1

            p = copy(q[1])
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
        n = length(q)
        v = [vitesses[rand(rng, 1:end)] for i in 1:n]

        new(n, q, v, ϵ)

    end

end

export HardSpheres

struct HardSpheres <: Particles

    n
    q
    v
    ϵ

    function HardSpheres(rng, n, ϵ)
    
        q = Vector{Float64}[]
    
        J4 = [1 0; 0 1]
    
        push!(q, [0.5000, 0.5000])
        push!(q, [0.1622, 0.4505])
        push!(q, [0.3112, 0.2290])
        push!(q, [0.5285, 0.9133])
        push!(q, [0.1656, 0.1524])
        push!(q, [0.6020, 0.8258])
        push!(q, [0.2630, 0.5383])
        push!(q, [0.3000, 0.8000])
        push!(q, [0.6892, 0.0782])
        push!(q, [0.7482, 0.4427])
    
        p = zeros(2)
    
        for klm in 11:n
            frein = 0
            overlap = 1
            while (overlap == 1) && frein < 10000
    
                p0 = [2ϵ + (1 - 4ϵ) * rand(rng), 2ϵ + (1 - 4ϵ) * rand(rng)]
                p = J4 * p0
                overlap2 = 0
                for subh in 1:(klm-1)
                    if norm(p - q[subh], 2) < 2ϵ
                        overlap2 = 1
                    end
                end
    
                overlap = overlap2
                frein += 1
            end
    
            frein == 10000 && (@error "Echec tirage initial")
    
            push!(q, p)
    
        end
    
        v = [ randn(rng, 2) for j in 1:n ]
    
        new( n, q, v, ϵ )
    
    end

end
