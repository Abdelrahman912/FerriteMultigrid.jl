@testset "Poisson Equation Example" begin
    configs = [
        PMultigridConfiguration(Galerkin(),DirectProjection()),
        PMultigridConfiguration(Rediscretization(DiffusionMultigrid(1.0)),DirectProjection()),
        PMultigridConfiguration(Galerkin(),StepProjection(1)),
        PMultigridConfiguration(Rediscretization(DiffusionMultigrid(1.0)),StepProjection(1)),
    ]
    ## 1D Poisson equation example ##
    for config in configs
        K, f, fe_space = poisson(1000, 2, 3)
        # 1. default configuration
        x, res = solve(K, f, fe_space,config)
        @test K * x ≈ f
    end

    ## 2D Poisson equation example ##
    for config in configs
        K, f, fe_space = poisson((100,100), 2, 3)
        # 1. default configuration
        x, res = solve(K, f, fe_space,config)
        @test K * x ≈ f
    end
end

