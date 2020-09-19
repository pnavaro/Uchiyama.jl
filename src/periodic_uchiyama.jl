function init_particles(n, ϵ, rng)

    q = [zeros(2) for i in 1:n]
    l1 = 0.5
    l2 = 0.5
    q[1] = [0.5,0.5]

    J4 = Matrix(I, 2, 2)

    for klm in 2:n

        frein = 0
        overlap = 1

        while (overlap == 1) && frein < 10000

            p0 = rand(2)
            p  = dot(J4, p0)
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
            @error "Echec tirage initial"
        end

        q[klm] = p

    end

    vitesses = [[1, 0], [0, 1], [-1, 0], [0, -1]]
    v = [vitesses[rand(1:end)] for i in 1:n]

    return q, v
end


function tempscoll(x, y, v, w, ϵ)

    deltat1 = Inf
    deltat2 = Inf

    if (y .- x)'v > 0 > (y .- x)'w

        dx = y[1] - x[1]
        dy = y[2] - x[2]
        adx = abs(dx)
        ady = abs(dy)
        adxy = abs(adx - ady)

        if (w[1] == 1 && v[1] == -1 && ady < 2ϵ)
            deltat1 = (2ϵ - dx - ady) / 2
            deltat2 = (- 2ϵ - dx + ady) / 2
        elseif (v[1] == 1 && w[1] == -1 && ady < 2ϵ)
            deltat1 = (2ϵ + dx - ady) / 2
            deltat2 = (- 2ϵ + dx + ady) / 2
        elseif (w[2] == 1 && v[1] == -1 && adx < 2ϵ)
            deltat1 = (2ϵ - adx - dy) / 2
            deltat2 = (- 2ϵ + adx - dy) / 2
        elseif (v[2] == 1 && w[2] == -1 && adx < 2ϵ)
            deltat1 = (2ϵ - adx + dy) / 2
            deltat2 = (- 2ϵ + adx + dy) / 2
        elseif (w[1] == 1 && v[2] == 1 && adxy < 2ϵ)
            deltat1 = (dy - dx - 2ϵ) / 2
            deltat2 = (dy - dx + 2ϵ) / 2
        elseif (v[1] == 1 && w[2] == 1 && adxy < 2ϵ)
            deltat1 = (-dy + dx - 2ϵ) / 2
            deltat2 = (-dy + dx + 2ϵ) / 2
        elseif (w[1] == 1 && v[2] == -1 && adxy < 2ϵ)
            deltat1 = (- dx - dy - 2ϵ) / 2
            deltat2 = (- dx - dy + 2ϵ) / 2
        elseif (v[1] == 1 && w[2] == -1 && adxy < 2ϵ)
            deltat1 = (dx + dy - 2ϵ) / 2
            deltat2 = (dx + dy + 2ϵ) / 2
        elseif (w[1] == -1 && v[2] == 1 && adxy < 2ϵ)
            deltat1 = (dx + dy + 2ϵ) / 2
            deltat2 = (dx + dy - 2ϵ) / 2
        elseif (v[1] == -1 && w[2] == 1 && adxy < 2ϵ)
            deltat1 = (- dx - dy + 2ϵ) / 2
            deltat2 = (- dx - dy - 2ϵ) / 2
        elseif (w[1] == -1 && v[2] == -1 && adxy < 2ϵ)
            deltat1 = (dx - dy + 2ϵ) / 2
            deltat2 = (dx - dy - 2ϵ) / 2
        elseif (v[1] == -1 && w[2] == -1 && adxy < 2ϵ)
            deltat1 = (- dx + dy + 2ϵ) / 2
            deltat2 = (- dx + dy - 2ϵ) / 2
        end

    end

    return min(deltat1, deltat2)

end


const offset = [[0, 0], [1,  0], [-1,  0], 
                [0, 1], [0, -1], [-1,  1], 
                [1, 1], [1, -1], [-1, -1]]

function compute_dt(i, q, v, ϵ, dt, fantome)

    n = length(q)

    for j in 1:n
        if i != j
            dt_local = Inf
            k = 1
            while (isinf(dt_local) && k < 10)
                dt_local = tempscoll(q[j] + offset[k], q[i], v[j], v[i], ϵ)
                k += 1
            end

            dt[i,j] = dt_local
            fantome[i,j] = k
        end
    end

end


struct PCollisionMatrix

    dt :: Array{Float64,2}
    dt_min :: Float64
    fantome :: Array{Int, 2}

    function PCollision(npart, q, v, ϵ)

        dt = zeros(Float64, (npart, npart))
        fantome = zeros(Int, (npart, npart))
        fill!(dt, Inf)
        fill!(fantome, 0)
        dt_min = Inf

        for k in 1:npart
            for l in (k+1):npart
                dt_local = Inf
                i = 1

                while (isinf(dt_local) && i < 10)
                    dt_local = tempscoll(q[l] + offset[i], q[k], v[l], v[k], ϵ)
                    i += 1
                end

                dt[k, l] = dt_local
                fantome[k, l] = i
            end
        end

        new( dt, dt_min, fantome )

    end

end

function dt_min_position(self)
    p = argmin(self.dt)
    self.dt_min = self.dt[p]
    self.fant = self.fantome[p]
    return p[1], p[2]
end

function reset!(dt, i)
    dt[i, :] .= Inf
    dt[:, i] .= Inf
end

struct Particles

    function Particles( npart, ϵ, rng)

        q, v = init_particles(npart, ϵ, iseed)
        collisions = PCollisionMatrix(npart, q, v, ϵ)
        new( npart, ϵ, q, v, collisions)

    end

end



function step(self)

    """ Compute small time step between two collisions """

    i1, i2 = self.collisions.dt_min_position

    tempsp = self.collisions.dt_min
    num_fant = self.collisions.fant-1

    self.q = (self.q + tempsp * self.v) % 1

    """ Algorithme de collision inter particules"""

    if self.v[i1]'self.v[i2] == 0

        self.v[[i1, i2]] = self.v[[i2, i1]]  # swap velocities

    elseif self.v[i1]'self.v[i2] == -1

        rot = [0 -1; 1 0]

        if (rot * self.v[i1])'(self.q[i2] + offset[num_fant] - self.q[i1]) < 0
            self.v[i1] = rot * self.v[i1]
            self.v[i2] = rot * self.v[i2]
        elseif (rot * self.v[i1])'(self.q[i2] + offset[num_fant] - self.q[i1]) > 0
            self.v[i1] = - rot * self.v[i1]
            self.v[i2] = - rot * self.v[i2]
        end

    end

    self.collisions.dt = self.collisions.dt - tempsp
    reset!(self.collisions, i1)
    reset!(self.collisions. i2)

    compute_dt(i1, self.q, self.v, self.ϵ, self.collisions.dt, self.collisions.fantome)

    compute_dt(i2, self.q, self.v, self.ϵ, self.collisions.dt, self.collisions.fantome)

    self.tcoll += tempsp

end

function run(self, nstep=15000, ctime=2)

    istep = 0
    while istep < nstep && self.tcoll < ctime

        istep += 1
        step(self)

        return istep, self.tcoll, self.etime

    end

end
