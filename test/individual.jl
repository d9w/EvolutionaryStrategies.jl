using EvolutionaryStrategies
using Cambrian
using Test

cfg = get_config("test.yaml")

@testset "Individual" begin
    ind = ESIndividual(cfg)
    @test ind.fitness[1] == -Inf

    @test rosenbrock(ind)[1] < 0
    @test sphere(ind)[1] < 0
end
