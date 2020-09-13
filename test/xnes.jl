using EvolutionaryStrategies
using Cambrian
using Test

cfg = get_config("test.yaml")

@testset "Exponential Natural Evolution Strategy" begin
    e = xNES(cfg, sphere)
    @test ~any(isnan.(e.state.μ))
    @test ~any(isnan.(e.state.B))

    population = EvolutionaryStrategies.xnes_init(e.cfg, e.state)
    @test length(population) == e.cfg.n_population
    for p in population
        @test all(e.fitness .== -Inf)
    end

    step!(e)
    for p in e.population
        ~@test any(e.fitness .== -Inf)
    end
    @test ~any(isnan.(e.state.μ))
    @test ~any(isnan.(e.state.B))
    best = e.elites[end]
    @test best.fitness[1] <= 0.0
    @test e.gen == 1

    println(T)
    @timev step!(e)
    new_best = e.elites[end]
    @test new_best.fitness[1] >= best.fitness[1]
end
