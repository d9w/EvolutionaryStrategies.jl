using EvolutionaryStrategies
using Cambrian
using Test

cfg = get_config("test.yaml")

function test_evo(T::DataType, fitness::Function)
    e = T(cfg, fitness)
    @test length(e.population) == e.cfg.n_population

    step!(e)
    best = e.elites[end]
    @test best.fitness[1] <= 0.0
    @test e.gen == 1

    println(T)
    @timev step!(e)
    new_best = e.elites[end]
    @test new_best.fitness[1] >= best.fitness[1]

    run!(e)
    final_best = e.elites[end]
    @test new_best.fitness[1] >= best.fitness[1]
end

@testset "Exponential Natural Evolution Strategy" begin
    test_evo(xNES, sphere)
    test_evo(xNES, rosenbrock)
end

# @testset "Separable Natural Evolution Strategy" begin
#     test_evo(SNES, sphere)
#     test_evo(SNES, rosenbrock)
# end

# @testset "CMA Evolution Strategy" begin
#     test_evo(CMAES, sphere)
#     test_evo(CMAES, rosenbrock)
# end
