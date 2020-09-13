export xNES, xnes_init, xnes_populate, populate

# adapted from
# https://github.com/francescoalemanno/NaturalES.jl/blob/master/src/exponential_nes.jl

import Cambrian: populate, evaluate, generation

function utilities(n::Integer)
    u = zeros(n)
    u = max.(0, log(n/2+1) .- log.(1:n))
    u ./ sum(u) .- 1/n
end

mutable struct xNESState <: ESState
    μ::Array{Float64}
    σ::Float64
    B::Array{Float64}
end

function xNESState(d::Int, n::Int; μ::Array{Float64}=randn(d),
                   σ::Float64=1.0,
                   A=diagm(fill(σ, d)))
    σ = abs(det(A))^(1/d)
    B = A ./ σ
    xNESState(μ, σ, B)
end

function xnes_config(d::Int)
    (n_population = 4+ceil(Int, log(3*d)),
     ημ = 1.0,
     ησ = (9 + 3 * log(d)) / (5 * d * sqrt(d)),
     ηB = 1.0,
     u = utilities(n_population))
end

mutable struct xNES <: Cambrian.AbstractEvolution
    config::NamedTuple
    logger::CambrianLogger
    population::Array{FloatIndividual}
    elites::Array{FloatIndividual}
    state::xNESState
    fitness::Function
    gen::Int
end

"generate initial population"
function initialize(itype::Type, cfg::NamedTuple)
end


function xnes_init(cfg::NamedTuple, state::ESState)
    population = Array{FloatIndividual}(undef, cfg.n_population)
    for i in 1:cfg.n_population
        genes = e.state.μ .+ e.state.σ .* e.state.B * randn(cfg.n_genes)
        population[i] = FloatIndividual(genes, -Inf*ones(cfg.d_fitness))
    end
    population
end

function xNES(cfg::NamedTuple, fitness::Function, state::xNESState; logfile=string("logs/", cfg.id, ".csv"))
    logger = CambrianLogger(logfile)
    population = xnes_init(cfg, state)
    elites = Array{FloatIndividual}(undef, cfg.n_elite)
    xNES(cfg, logger, population, elites, state, fitness, 0)
end

function xNES(cfg::NamedTuple, fitness::Function; logfile=string("logs/", cfg.id, ".csv"))
    logger = CambrianLogger(logfile)
    cfg = merge(xnes_config(cfg.n_genes), cfg)
    state = xNESState(cfg.n_genes, cfg.n_population)
    xNES(fitness, state, cfg; logfile=logfile)
end

"generate next population"
function xnes_populate(e::xNES)
    for i in eachindex(e.population)
        e.population[i].genes .= e.state.μ .+ e.state.σ .* e.state.B * randn(cfg.n_genes)
        e.population[i].fitness .= -Inf
    end
end

"update NES state"
function xnes_generation(e::xNES)
    d = e.cfg.n_genes
    n = e.cfg.n_population

    # copy population information
    F = zeros(n)
    Z = zeros(d, n)
    for i in eachindex(population)
        F[i] = e.population[i].fitness[1]
        Z[:, i] = e.population[i].genes
    end
    idx = sortperm(F)

    # compute gradients
    ∇δ = zeros(d)
    ∇M = zeros(d,d)
    for i in 1:n
        j = idx[i]
        for k2 in 1:d
            ∇δ[k2] += e.cfg.u[i] * Z[k2,j]
            for k1 in 1:d
                ∇M[k1,k2] += u[i] * (Z[k1,j] * Z[k2,j] - (k1==k2))
            end
        end
    end
    ∇σ = tr(∇M) / d
    for i in 1:d
        ∇M[i,i] -=  ∇σ
    end

    # update state variables
    e.state.μ = e.state.μ .+ e.cfg.ημ .* e.state.σ .* e.state.B * ∇δ
    e.state.σ = e.state.σ * exp(e.cfg.ησ/2 * ∇σ)
    e.state.B = e.state.B * exp(e.cfg.ηB/2 .* ∇M)
    Cambrian.elites_generation(e)
end

populate(e::xNES) = xnes_populate(e)
evaluate(e::xNES) = fitness_evaluate(e, e.fitness)
generation(e::xNES) = xnes_generation(e)
