import Cambrian.get_config
using YAML
using Dates

export ESIndividual, AbstractESIndividual, get_config

abstract type AbstractESIndividual <: Cambrian.Individual end

"FloatIndividual : example Individual class using a floating point genotype in [0, 1]"
struct ESIndividual <: AbstractESIndividual
    genes::Array{Float64}
    fitness::Array{Float64}
end

function ESIndividual(cfg::NamedTuple)
    ESIndividual(rand(cfg.n_genes), -Inf*ones(cfg.d_fitness))
end

function ESIndividual(st::String)
    dict = ind_parse(st)
    ESIndividual(Float64.(dict["genes"]), Float64.(dict["fitness"]))
end

## utils

function Cambrian.get_config(cfg_file::String; kwargs...)
    cfg = YAML.load_file(cfg_file)
    for (k, v) in kwargs
        cfg[String(k)] = v
    end
    # generate id, use date if no existing id
    if ~(:id in keys(cfg))
        cfg["id"] = replace(string(Dates.now()), r"[]\.:]" => "_")
    end
    get_config(cfg)
end
