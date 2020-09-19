using EvolutionaryStrategies
using Cambrian
using Test

cfg = get_config("test.yaml")

@testset "Exponential Natural Evolution Strategy" begin
    center = ones(cfg.n_genes)
    e = xNES(cfg, i->sphere(i; center=ones(cfg.n_genes)))
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

    for i in 1:10
        step!(e)
        new_best = e.elites[end]
        @test new_best.fitness[1] >= best.fitness[1]
        @test ~any(isnan.(e.state.μ))
        @test ~any(isnan.(e.state.B))
        @test all([p.fitness[1] .!= -Inf for p in e.population])
    end
    @timev step!(e)
    @test e.elites[end].fitness[1] < 0.1
end
