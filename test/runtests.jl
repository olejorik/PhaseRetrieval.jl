using PhaseRetrieval
# using SampledDomains
import Pkg
Pkg.add("`https://github.com/olejorik/SampledDomains.jl");
using Test

@testset "utils.jl" begin
    # Write your tests here.
    q1d = PhaseRetrieval.quadratic(5)
    q2d = PhaseRetrieval.quadratic((5,5))
    @test q1d == q2d ==  [
    8.0  5.0  4.0  5.0  8.0
    5.0  2.0  1.0  2.0  5.0
    4.0  1.0  0.0  1.0  4.0
    5.0  2.0  1.0  2.0  5.0
    8.0  5.0  4.0  5.0  8.0]

    dm = PhaseRetrieval.diskmatrix(10, .8)
    @test dm[1] == [
        0  0  0  0  0  0  0  0  0  0
        0  0  0  0  1  1  0  0  0  0
        0  0  1  1  1  1  1  1  0  0
        0  0  1  1  1  1  1  1  0  0
        0  1  1  1  1  1  1  1  1  0
        0  1  1  1  1  1  1  1  1  0
        0  0  1  1  1  1  1  1  0  0
        0  0  1  1  1  1  1  1  0  0
        0  0  0  0  1  1  0  0  0  0
        0  0  0  0  0  0  0  0  0  0 
            ]
    # @test dm[2] == [
    #     NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
    #     NaN  NaN  NaN  NaN    1    1  NaN  NaN  NaN  NaN
    #     NaN  NaN    1    1    1    1    1    1  NaN  NaN
    #     NaN  NaN    1    1    1    1    1    1  NaN  NaN
    #     NaN    1    1    1    1    1    1    1    1  NaN
    #     NaN    1    1    1    1    1    1    1    1  NaN
    #     NaN  NaN    1    1    1    1    1    1  NaN  NaN
    #     NaN  NaN    1    1    1    1    1    1  NaN  NaN
    #     NaN  NaN  NaN  NaN    1    1  NaN  NaN  NaN  NaN
    #     NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
    # ]    

    @test PhaseRetrieval.lineararray(-1:.5:1, -1:.25:1, 3, 2,1) == [
        -4.0  -2.5  -1.0  0.5  2.0
        -3.5  -2.0  -0.5  1.0  2.5
        -3.0  -1.5   0.0  1.5  3.0
        -2.5  -1.0   0.5  2.0  3.5
        -2.0  -0.5   1.0  2.5  4.0
        -1.5   0.0   1.5  3.0  4.5
        -1.0   0.5   2.0  3.5  5.0
        -0.5   1.0   2.5  4.0  5.5
         0.0   1.5   3.0  4.5  6.0
    ]

    @test PhaseRetrieval.lineararray(5, .5, 1, -1.5) == [
        0.0  0.5  1.0  1.5  2.0
        1.0  1.5  2.0  2.5  3.0
        2.0  2.5  3.0  3.5  4.0
        3.0  3.5  4.0  4.5  5.0
        4.0  4.5  5.0  5.5  6.0
    ]

    @test PhaseRetrieval.subdivide_sum(ones(10,10),5) == [
        25.0  25.0
        25.0  25.0
    ]

    a = SampledDomains.CartesianDomain2D(1:5, -2.5:.5:2.5)
    ap = PhaseRetrieval.aperture(a, 2.5, (3,-0.1))
    b = phwrap.(ap[2])
    # b  == [
    #     NaN  NaN    NaN    NaN    NaN
    #     NaN  NaN    NaN    NaN    NaN
    #     NaN  NaN    NaN    NaN    NaN
    #     NaN  NaN      1.0  NaN    NaN
    #     NaN    1.0    1.0    1.0  NaN
    #     NaN    1.0    1.0    1.0  NaN
    #     NaN    1.0    1.0    1.0  NaN
    #     NaN  NaN      1.0  NaN    NaN
    #     NaN  NaN    NaN    NaN    NaN
    #     NaN  NaN    NaN    NaN    NaN
    #     NaN  NaN    NaN    NaN    NaN
    #  ]
    
    @test isnan.(b) == [
        1  1  1  1  1
        1  1  1  1  1
        1  1  1  1  1
        1  1  0  1  1
        1  0  0  0  1
        1  0  0  0  1
        1  0  0  0  1
        1  1  0  1  1
        1  1  1  1  1
        1  1  1  1  1
        1  1  1  1  1
    ]

   a =  PhaseRetrieval.tile(reshape(1:60,(6,10)), 2)
   @test a[:,:,1] == [ 1 7;2 8] && a[:,:,15] == [53 59; 54 60]
   

end
