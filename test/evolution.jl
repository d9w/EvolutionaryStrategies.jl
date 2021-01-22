using EvolutionaryStrategies
using Cambrian
using Test

cfg = get_config("test.yaml")

function test_evo(T::DataType, fitness::Function)
    println("Evolving ", T, " on ", fitness)
    e = T(cfg, fitness)
    @test length(e.population) == e.config.n_population

    step!(e)
    best = e.elites[end]
    @test best.fitness[1] <= 0.0
    @test e.gen == 1

    step!(e)
    new_best = e.elites[end]
    @test new_best.fitness[1] >= best.fitness[1]

    @timev run!(e)
    final_best = e.elites[end]
    @test final_best.fitness[1] >= new_best.fitness[1]

    println("Final fitness on ", fitness, ": ", final_best.fitness[1])
    return final_best.fitness[1]
end

@testset "Evolution xNES" begin
    @test abs(test_evo(xNES, sphere)) < 1e-5
    @test abs(test_evo(xNES, rosenbrock)) < 1.0
end

@testset "Evolution sNES" begin
    @test abs(test_evo(sNES, sphere)) < 1e-5
    @test abs(test_evo(sNES, rosenbrock)) < 1.0
end

@testset "CMA Evolution Strategy" begin
    @test abs(test_evo(CMAES, sphere)) < 1e-5
    @test abs(test_evo(CMAES, rosenbrock)) < 1.0
end
