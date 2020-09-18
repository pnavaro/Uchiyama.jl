using Uchiyama
using Test

@testset "Uchiyama.jl" begin

n = 40
epsilon = 1 / n

sim  = Trajectory(n, epsilon)
visu = Visualisation(sim, 1000)
anim(visu)

end
