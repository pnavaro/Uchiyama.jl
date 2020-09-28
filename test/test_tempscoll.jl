@testset "tempscoll" begin

    n = 10
    ϵ = 0.1
    rng = MersenneTwister(42)

    p = SquareParticles(rng, n, ϵ)
    
    q = [[0.5       , 0.5       ],
         [0.20006862, 0.38139954],
         [0.32267135, 0.72687046],
         [0.68044674, 0.78095695],
         [0.77473372, 0.51989917],
         [0.4949439 , 0.23201753],
         [0.74204115, 0.28248482],
         [0.63101615, 0.3650337 ],
         [0.5846235 , 0.64552143],
         [0.45656816, 0.76644607]]
    
    v = [[ 0.,   1.],
         [ 1.,   0.],
         [-1.,   0.],
         [ 0.,  -1.],
         [ 1.,   0.],
         [ 0.,   1.],
         [ 0.,   1.],
         [-1.,   0.],
         [-1.,   0.],
         [ 0.,  -1.]]

    for i in 1:n
        p.q[i] = q[i]
        p.v[i] = v[i]
    end
    
    @test Uchiyama.compute_dt(p, 4, 1, 1) ≈ 0.130701845
    @test Uchiyama.compute_dt(p, 9, 1, 1) ≈ 0.015072465
    @test Uchiyama.compute_dt(p, 4, 2, 1) ≈ 0.339967765
    @test Uchiyama.compute_dt(p, 6, 2, 1) ≈ 0.122128645
    @test Uchiyama.compute_dt(p, 8, 2, 1) ≈ 0.123656685
    @test Uchiyama.compute_dt(p, 7, 3, 3) ≈ 0.412507920
    @test Uchiyama.compute_dt(p, 6, 4, 1) ≈ 0.267221131
    @test Uchiyama.compute_dt(p, 7, 4, 1) ≈ 0.180033270
    @test Uchiyama.compute_dt(p, 8, 5, 2) ≈ 0.405573952
    @test Uchiyama.compute_dt(p, 9, 5, 2) ≈ 0.367756020
    @test Uchiyama.compute_dt(p, 8, 6, 1) ≈ 0.034544210
    @test Uchiyama.compute_dt(p, 8, 7, 7) ≈ 0.885761939

end
