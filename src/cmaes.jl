export CMAES, CMAES_populate, CMAES_generation, populate, evaluate, generation

# CMA-ES implementation based on pycma: https://github.com/CMA-ES/pycma
# To install it: pip install cma

import Cambrian: populate, evaluate, generation
using PyCall

mutable struct CMAES <: Cambrian.AbstractEvolution
    config::NamedTuple
    logger::CambrianLogger
    population::Array{AbstractESIndividual}
    elites::Array{AbstractESIndividual}
    state::Array{Float64}
    fitness::Function
    gen::Int
    pycma::PyObject
end

function CMAES(cfg::NamedTuple, fitness::Function, state::Array{Float64}; T::Type=ESIndividual ,logfile=string("logs/", cfg.id, ".csv"))
    cma = pyimport("cma")
    pycma = cma.CMAEvolutionStrategy(state, cfg.m_rate)
    logger = CambrianLogger(logfile)
    e = CMAES(cfg, logger, [], [], state, fitness, 0, pycma)
    CMAES_populate(e, T)
    e.elites = deepcopy([e.population[i] for i in 1:cfg.n_elite])
    e
end

function CMAES(cfg::NamedTuple, fitness::Function;  T::Type=ESIndividual, logfile=string("logs/", cfg.id, ".csv"))
    logger = CambrianLogger(logfile)
    state = zeros(cfg.n_genes)
    CMAES(cfg, fitness, state; T=T, logfile=logfile)
end

"get next population from pycma and create Cambrian individuals"
function CMAES_populate(e::CMAES, T::Type)
    solutions = e.pycma.ask()
    e.config = merge((n_population=length(solutions), ), e.config)
    e.population=[T(s, -Inf*ones(e.config.d_fitness)) for s in solutions]
end

"get next population from pycma and create Cambrian individuals"
function CMAES_populate(e::CMAES)
    T = typeof(e.population[1])
    CMAES_populate(e, T)
end

"update pycma, called after populate and evaluate"
function CMAES_generation(e::CMAES)
    genes = (i->i.genes).(e.population)
    # Cambrian maximises, pycma minimises, so -1 factor on fitness
    fitnesses = (i->-1*i.fitness[1]).(e.population)

    e.pycma.tell(genes, fitnesses)
    Cambrian.elites_generation(e)
end

populate(e::CMAES) = CMAES_populate(e)
evaluate(e::CMAES) = fitness_evaluate(e, e.fitness)
generation(e::CMAES) = CMAES_generation(e)
