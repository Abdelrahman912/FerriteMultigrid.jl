import FerriteMultigrid: init, AMGCoarseSolver

@testset "MultiLevel Parameters" begin
    K, f, fe_space = poisson(3, 2, 3)
    ## SA-AMG as coarse solver
    pmgsolver = init(K, f, fe_space, SmoothedAggregationCoarseSolver, presmoother=GaussSeidel(; iter=4), postsmoother=GaussSeidel(; iter=2))
    ml = pmgsolver.ml
    @test ml.coarse_solver isa AMGCoarseSolver
    @test ml.coarse_solver.alg isa SmoothedAggregationAMG
    @test ml.levels |> length == 1
    @test ml.presmoother isa GaussSeidel
    @test ml.presmoother.iter == 4
    @test ml.postsmoother isa GaussSeidel
    @test ml.postsmoother.iter == 2

    ## RS-AMG as coarse solver
    pmgsolver = init(K, f, fe_space, RugeStubenCoarseSolver, presmoother=GaussSeidel(; iter=2), postsmoother=GaussSeidel(; iter=4))
    ml = pmgsolver.ml
    @test ml.coarse_solver isa AMGCoarseSolver
    @test ml.coarse_solver.alg isa RugeStubenAMG
    @test ml.levels |> length == 1
    @test ml.presmoother isa GaussSeidel
    @test ml.presmoother.iter == 2
    @test ml.postsmoother isa GaussSeidel
    @test ml.postsmoother.iter == 4

    ## Direct Solver as coarse solver (e.g. Pinv)
    pmgsolver = init(K, f, fe_space, Pinv)
    ml = pmgsolver.ml
    @test ml.coarse_solver isa Pinv

end
