using PhaseRetrieval
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
end
