export xNES, xnes_init, xnes_populate, xnes_generation, populate, evaluate, generation

# adapted from
# https://github.com/francescoalemanno/NaturalES.jl/blob/master/src/exponential_nes.jl

import Cambrian: populate, evaluate, generation

mutable struct xNESState <: ESState
    μ::Array{Float64}
    σ::Float64
    B::Array{Float64}
    u::Array{Float64}
    s::Array{Float64}
end

function xNESState(d::Int, n::Int;
                   μ::Array{Float64}=randn(d),
                   σ::Float64=1.0,
                   A=diagm(fill(σ, d)))
    σ = abs(det(A))^(1/d)
    B = A ./ σ
    u = zeros(n)
    u = max.(0, log(n/2+1) .- log.(1:n))
    u .= u ./ sum(u) .- 1/n
    s = randn(d, n)
    xNESState(μ, σ, B, u, s)
end

function xnes_config(d::Int; n_population = 4+ceil(Int, log(3*d)),
                     ημ = 1.0, ησ = (9 + 3 * log(d)) / (5 * d * sqrt(d)),
                     ηB = ησ)
    (n_population = n_population, ημ = ημ, ησ = ησ, ηB = ηB)
end

mutable struct xNES <: Cambrian.AbstractEvolution
    config::NamedTuple
    logger::CambrianLogger
    population::Array{AbstractESIndividual}
    elites::Array{AbstractESIndividual}
    state::xNESState
    fitness::Function
    gen::Int
end

function xnes_init(cfg::NamedTuple, state::ESState; T::Type=ESIndividual)
    population = Array{AbstractESIndividual}(undef, cfg.n_population)
    for i in 1:cfg.n_population
        genes = state.μ .+ state.σ .* (state.B * view(state.s, :, i))
        population[i] = T(genes, -Inf*ones(cfg.d_fitness))
    end
    population
end

function xNES(cfg::NamedTuple, fitness::Function, state::xNESState; T::Type=ESIndividual ,logfile=string("logs/", cfg.id, ".csv"))
    logger = CambrianLogger(logfile)
    population = xnes_init(cfg, state, T=T)
    elites = deepcopy([population[i] for i in 1:cfg.n_elite])
    xNES(cfg, logger, population, elites, state, fitness, 0)
end

"""
Exponential Natural Evolutionary Strategies default constructor
Will use a random starting point using N(0, 1)
If cfg contains keys from snes_config, these will be used
To provide initial state, see xNES(cfg, fitness, state)
"""
function xNES(cfg::NamedTuple, fitness::Function;  T::Type=ESIndividual, logfile=string("logs/", cfg.id, ".csv"))
    logger = CambrianLogger(logfile)
    cfg = merge(xnes_config(cfg.n_genes), cfg)
    state = xNESState(cfg.n_genes, cfg.n_population)
    xNES(cfg, fitness, state; T=T, logfile=logfile)
end

"generate next population"
function xnes_populate(e::xNES)
    for i in eachindex(e.population)
        e.population[i].genes .= e.state.μ .+ e.state.σ .* (e.state.B * view(e.state.s, :, i))
        e.population[i].fitness .= -Inf
    end
end

"update NES state, called after populate and evaluate"
function xnes_generation(e::xNES)
    d = e.config.n_genes
    n = e.config.n_population

    # copy population information
    F = zeros(n)
    for i in eachindex(e.population)
        F[i] = -e.population[i].fitness[1]
    end
    idx = sortperm(F)

    # compute gradients
    ∇δ = zeros(d)
    ∇M = zeros(d,d)
    for i in 1:n
        j = idx[i]
        for k2 in 1:d
            ∇δ[k2] += e.state.u[i] * e.state.s[k2,j]
            for k1 in 1:d
                ∇M[k1,k2] += e.state.u[i] * (e.state.s[k1,j] * e.state.s[k2,j] - (k1==k2))
            end
        end
    end
    ∇σ = tr(∇M) / d
    for i in 1:d
        ∇M[i,i] -=  ∇σ
    end

    # update state variables
    e.state.μ .+= e.config.ημ .* e.state.σ .* (e.state.B * ∇δ)
    e.state.σ *= exp(e.config.ησ/2 * ∇σ)
    e.state.B = e.state.B * exp(e.config.ηB/2 .* ∇M)
    randn!(e.state.s)
    Cambrian.elites_generation(e)
end

populate(e::xNES) = xnes_populate(e)
evaluate(e::xNES) = fitness_evaluate(e, e.fitness)
generation(e::xNES) = xnes_generation(e)
