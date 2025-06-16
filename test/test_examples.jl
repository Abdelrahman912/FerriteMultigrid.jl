@testset "Poisson Equation Example" begin
    K, f, fe_space = poisson(3, 2, 3)
    u = K\f
    x, res = solve(K, f, fe_space)
    @test x ≈ u
end
