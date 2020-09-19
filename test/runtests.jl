using EvolutionaryStrategies
using Cambrian
using Test

function sphere(i::Individual; center=zeros(length(i.genes)))
    [-sum((i.genes .- center).^2)]
end

function rosenbrock(i::Individual)
    x = i.genes
    y = -(sum([(1.0 - x[i])^2 + 100.0 * (x[i+1] - x[i]^2)^2
               for i in 1:(length(x)-1)]))
    [y]
end

# include("individual.jl")
include("xnes.jl")
# include("evolution.jl")
