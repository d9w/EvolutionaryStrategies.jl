using EvolutionaryStrategies
using Cambrian
using Test

cfg = get_config("test.yaml")

@testset "Exponential Natural Evolution Strategy" begin
    e = xNES(cfg, sphere)
    @test ~any(isnan.(e.state.μ))
    @test ~any(isnan.(e.state.B))

    population = EvolutionaryStrategies.xnes_init(e.config, e.state)
    @test length(population) == e.config.n_population
    @test all([p.fitness[1] .== -Inf for p in population])

    step!(e)
    @test all([p.fitness[1] .!= -Inf for p in e.population])

    @test ~any(isnan.(e.state.μ))
    @test ~any(isnan.(e.state.B))
    best = e.elites[end]
    @test best.fitness[1] <= 0.0
    @test e.gen == 1

    @timev step!(e)
    new_best = e.elites[end]
    @test new_best.fitness[1] >= best.fitness[1]
end
