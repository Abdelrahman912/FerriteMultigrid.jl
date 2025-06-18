@testset "Poisson Equation Example" begin
    K, f, fe_space = poisson(100, 2, 3)
    #diff_problem = DiffusionMultigrid(1.0)
    u = K\f
    x, res = solve(K, f, fe_space)
    @test x ≈ u
end

