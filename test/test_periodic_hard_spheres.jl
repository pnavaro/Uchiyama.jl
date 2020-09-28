@testset "periodic hard-spheres" begin

ϵ   = 0.1

q = [[0.5,  0.8], [0.8,  0.5]]
v = [[0.0, -1.0], [-1.0, 0.0]]

p  = HardSpheres( q, v, ϵ )

pc = PeriodicCollisions( p )

@show step!( p, pc )
@show p.v
@show step!( p, pc )
@show p.q
@show p.v

end
