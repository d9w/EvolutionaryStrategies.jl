export RandomSearch, RandomSearch_populate, RandomSearch_generation, populate, evaluate, generation

# CMA-ES implementation based on pycma: https://github.com/CMA-ES/pycma
# To install it: pip install cma

import Cambrian: populate, evaluate, generation
using PyCall
using Random

# With x random vector: new_genes = ax + b

mutable struct RandomSearch <: Cambrian.AbstractEvolution
    config::NamedTuple
    logger::CambrianLogger
    population::Array{AbstractESIndividual}
    elites::Array{AbstractESIndividual}
    fitness::Function
    gen::Int
    a::Array{Float64}
    b::Array{Float64}
end

function RandomSearch(cfg::NamedTuple, fitness::Function, a::Array{Float64}, b::Array{Float64}; T::Type=ESIndividual ,logfile=string("logs/", cfg.id, ".csv"))
    logger = CambrianLogger(logfile)
    elite = T(cfg)
    e = RandomSearch(cfg, logger, [], [], fitness, 0, a, b)
    e.config = merge((n_population=n_population = 4+3*ceil(Int, log(cfg.n_genes)), ), e.config)
    RandomSearch_populate(e, T)
    e.elites = deepcopy([e.population[i] for i in 1:cfg.n_elite])
    e
end

function RandomSearch(cfg::NamedTuple, fitness::Function, a::Float64, b::Float64; T::Type=ESIndividual ,logfile=string("logs/", cfg.id, ".csv"))
    arr_a = ones(cfg.n_genes) * a
    arr_b = ones(cfg.n_genes) * b
    RandomSearch(cfg, fitness, arr_a, arr_b; T=T, logfile=logfile)
end

function RandomSearch(cfg::NamedTuple, fitness::Function, a::Int64, b::Int64; T::Type=ESIndividual ,logfile=string("logs/", cfg.id, ".csv"))
    RandomSearch(cfg, fitness, a*1.0, b*1.0; T=T, logfile=logfile)
end

function RandomSearch(cfg::NamedTuple, fitness::Function; T::Type=ESIndividual ,logfile=string("logs/", cfg.id, ".csv"))
    RandomSearch(cfg, fitness, cfg.a, cfg.b; T=T, logfile=logfile)
end

"get next population from pycma and create Cambrian individuals"
function RandomSearch_populate(e::RandomSearch, T::Type)
    e.population=[]
    for i in 1:e.config.n_population
        genes = rand(Float64, e.config.n_genes) .* e.config.a .+ e.config.b
        indiv = T(genes, -Inf*ones(e.config.d_fitness))
        push!(e.population, indiv)
    end
    e
end

"get next population from pycma and create Cambrian individuals"
function RandomSearch_populate(e::RandomSearch)
    T = typeof(e.population[1])
    RandomSearch_populate(e, T)
end

"update pycma, called after populate and evaluate"
function RandomSearch_generation(e::RandomSearch)
    Cambrian.elites_generation(e)
end

populate(e::RandomSearch) = RandomSearch_populate(e)
evaluate(e::RandomSearch) = fitness_evaluate(e, e.fitness)
generation(e::RandomSearch) = RandomSearch_generation(e)
