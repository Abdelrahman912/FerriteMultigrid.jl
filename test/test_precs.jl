@testset "LinearSolvePrecs" begin
    # TODO: more tests after merge
    A, b, fe_space = poisson(3, 2, 3)

    prob = LinearProblem(A, b)
    prec_builder = PolynomialMultigridPreconBuilder(fe_space, Pinv)
    strategy = KrylovJL_CG(precs=prec_builder)
    @test A*solve(prob, strategy, atol=1.0e-14) ≈ b rtol = 1.0e-8

end
