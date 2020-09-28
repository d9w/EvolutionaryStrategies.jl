export sNES, snes_init, snes_populate, snes_generation, populate, evaluate, generation

# adapted from
# https://github.com/francescoalemanno/NaturalES.jl/blob/master/src/separable_nes.jl

import Cambrian: populate, evaluate, generation

mutable struct sNESState <: ESState
    μ::Array{Float64}
    σ::Array{Float64}
    u::Array{Float64}
    s::Array{Float64}
end

function sNESState(d::Int, n::Int;
                   μ::Array{Float64}=randn(d),
                   σ::Array{Float64}=ones(d))
    u = zeros(n)
    u = max.(0, log(n/2+1) .- log.(1:n))
    u .= u ./ sum(u) .- 1/n
    s = randn(d, n)
    sNESState(μ, σ, u, s)
end

function snes_config(d::Int; n_population = 4+ceil(Int, log(3*d)),
                     ημ = 1.0, ησ = (3 + log(d)) / (5 * sqrt(d)))
    (n_population = n_population, ημ = ημ, ησ = ησ)
end

mutable struct sNES <: Cambrian.AbstractEvolution
    config::NamedTuple
    logger::CambrianLogger
    population::Array{AbstractESIndividual}
    elites::Array{AbstractESIndividual}
    state::sNESState
    fitness::Function
    gen::Int
end

function snes_init(cfg::NamedTuple, state::ESState; T::Type=ESIndividual)
    population = Array{AbstractESIndividual}(undef, cfg.n_population)
    for i in 1:cfg.n_population
        genes = state.μ .+ state.σ .* view(state.s, :, i)
        population[i] = T(genes, -Inf*ones(cfg.d_fitness))
    end
    population
end

function sNES(cfg::NamedTuple, fitness::Function, state::sNESState; T::Type=ESIndividual, logfile=string("logs/", cfg.id, ".csv"))
    logger = CambrianLogger(logfile)
    population = snes_init(cfg, state, T=T)
    elites = deepcopy([population[i] for i in 1:cfg.n_elite])
    sNES(cfg, logger, population, elites, state, fitness, 0)
end

"""
Separable Natural Evolutionary Strategies default constructor
Will use a random starting point using N(0, 1)
If cfg contains keys from snes_config, these will be used
To provide initial state, see sNES(cfg, fitness, state)
"""
function sNES(cfg::NamedTuple, fitness::Function; logfile=string("logs/", cfg.id, ".csv"))
    logger = CambrianLogger(logfile)
    cfg = merge(snes_config(cfg.n_genes), cfg)
    state = sNESState(cfg.n_genes, cfg.n_population)
    sNES(cfg, fitness, state; logfile=logfile)
end

"generate next population"
function snes_populate(e::sNES)
    for i in eachindex(e.population)
        e.population[i].genes .= e.state.μ .+ e.state.σ .* view(e.state.s, :, i)
        e.population[i].fitness .= -Inf
    end
end

"update NES state, called after populate and evaluate"
function snes_generation(e::sNES)
    d = e.config.n_genes
    n = e.config.n_population

    # copy population information
    F = zeros(n)
    for i in eachindex(e.population)
        F[i] = -e.population[i].fitness[1]
    end
    idx = sortperm(F)

    # compute gradients
    ∇μ = zeros(d)
    ∇σ = zeros(d)
    for i in 1:n
        j = idx[i]
        ∇μ .+= e.state.u[i] .* e.state.s[:, j]
        ∇σ .+= e.state.u[i] .* (e.state.s[:, j].^2 .- 1.0)
    end

    # update state variables
    e.state.μ .+= e.config.ημ .* e.state.σ .* ∇μ
    e.state.σ .*= exp.(e.config.ησ/2 .* ∇σ)
    randn!(e.state.s)
    Cambrian.elites_generation(e)
end

populate(e::sNES) = snes_populate(e)
evaluate(e::sNES) = fitness_evaluate(e, e.fitness)
generation(e::sNES) = snes_generation(e)
