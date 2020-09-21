using Test

include("../src/periodic_uchiyama.jl")

n = 10

ϵ = 0.1


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


v = [[ 0.,   1.] ,
     [ 1.,   0.] ,
     [-1.,   0.] ,
     [ 0.,  -1.] ,
     [ 1.,   0.] ,
     [ 0.,   1.] ,
     [ 0.,   1.] ,
     [-1.,   0.] ,
     [-1.,   0.] ,
     [ 0.,  -1.]]

@test tempscoll(q[4] + offset[1], q[1], v[4], v[1], ϵ) ≈ 0.130701845
@test tempscoll(q[9] + offset[1], q[1], v[9], v[1], ϵ) ≈ 0.015072465
@test tempscoll(q[4] + offset[1], q[2], v[4], v[2], ϵ) ≈ 0.339967765
@test tempscoll(q[6] + offset[1], q[2], v[6], v[2], ϵ) ≈ 0.122128645
@test tempscoll(q[8] + offset[1], q[2], v[8], v[2], ϵ) ≈ 0.123656685
@test tempscoll(q[7] + offset[3], q[3], v[7], v[3], ϵ) ≈ 0.412507920
@test tempscoll(q[6] + offset[1], q[4], v[6], v[4], ϵ) ≈ 0.267221131
@test tempscoll(q[7] + offset[1], q[4], v[7], v[4], ϵ) ≈ 0.180033270
@test tempscoll(q[8] + offset[2], q[5], v[8], v[5], ϵ) ≈ 0.405573952
@test tempscoll(q[9] + offset[2], q[5], v[9], v[5], ϵ) ≈ 0.367756020
@test tempscoll(q[8] + offset[1], q[6], v[8], v[6], ϵ) ≈ 0.034544210
@test tempscoll(q[8] + offset[7], q[7], v[8], v[7], ϵ) ≈ 0.885761939


c = PCollisionMatrix( n, q, v, epsilon)

ref = [0 9 9 1 9 9 9 9 1 1;
       0 0 9 1 9 1 9 1 9 1;
       0 0 0 9 9 9 3 9 9 6;
       0 0 0 0 9 1 1 9 9 9;
       0 0 0 0 0 9 9 2 2 9;
       0 0 0 0 0 0 9 1 9 1;
       0 0 0 0 0 0 0 7 9 9;
       0 0 0 0 0 0 0 0 9 9;
       0 0 0 0 0 0 0 0 0 1;
       0 0 0 0 0 0 0 0 0 0]

@show ref' + ref 
@show c.fantome

@test all( (ref' + ref) .==  c.fantome )

step!(n, ϵ, q, v, c)

