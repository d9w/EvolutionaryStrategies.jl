# EvolutionaryStrategies.jl

[![Build Status](https://travis-ci.org/d9w/EvolutionaryStrategies.jl.svg?branch=master)](https://travis-ci.org/d9w/EvolutionaryStrategies.jl) [![Coverage Status](https://coveralls.io/repos/d9w/EvolutionaryStrategies.jl/badge.svg?branch=master)](https://coveralls.io/r/d9w/EvolutionaryStrategies.jl?branch=master) [![codecov](https://codecov.io/gh/d9w/EvolutionaryStrategies.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/d9w/EvolutionaryStrategies.jl)

Julia implementations of Natural Evolutionary Strategies and CMA-ES for the
[Cambrian.jl](https://github.com/d9w/Cambrian.jl) framework.

<img src="es.gif" width="400px" height="auto">

## Installation

`EvolutionaryStrategies.jl` can be installed through the Julia package manager:

```julia
pkg> add EvolutionaryStrategies
```

Tests are also provided:

```julia
pkg> test EvolutionaryStrategies
```

## Algorithms

`EvolutionaryStrategies` implements the Exponential and Separable Natural
Evolutionary Strategies, as described in:

Wierstra, D., Schaul, T., Glasmachers, T., Sun, Y., Peters, J., & Schmidhuber,
J. (2014). Natural evolution strategies. The Journal of Machine Learning
Research, 15(1), 949-980.
[pdf](https://www.ini.rub.de/PEOPLE/glasmtbl/paper/wierstra2014.pdf)

and the Covariance Matrix Adaptation Evolutionary Strategy (CMA-ES), as
described in:

Hansen, N., & Ostermeier, A. (2001). Completely derandomized self-adaptation in
evolution strategies. Evolutionary computation, 9(2), 159-195.
[pdf](http://staff.elka.pw.edu.pl/~jarabas/ALHE/CMAES.pdf)

with additional implementation details based on
[pycma](https://github.com/CMA-ES/pycma).

## Usage

The function to optimize must first be defined:

```julia
fitness(i::Individual) = -sum(i.genes .^ 2)
```

Note that Cambrian by default **maximizes** objective fitness, which is common
in neuroevolution and genetic programming. Evolutionary Strategies often
*minimize* objective functions, but for coherence with Cambrian,
`EvolutionaryStrategies.jl` maximizes. For objective function definitions, you
must negate fitness if aiming to minimize, as demonstrated above.

Then, create and run the desired ES:

```julia
cfg = get_config("cfg/cma-es.yaml")
es = CMAES(cfg, fitness)
run!(es)
```

Examples can be found in the `scripts/` directory.

## Similar packages

Other Evolutionary Strategies resources, notably other Julia packages:

+ [pycma](https://github.com/CMA-ES/pycma)
+ [CMAEvolutionStrategy.jl](https://github.com/jbrea/CMAEvolutionStrategy.jl)
+ [NaturalES.jl](https://github.com/francescoalemanno/NaturalES.jl)
+ [BlackBoxOptim.jl](https://github.com/robertfeldt/BlackBoxOptim.jl)

## Development

Next steps (pull requests are greatly appreciated):

+ ~~Separable NES~~
+ CMA-ES
+ Multi-objective
+ Constraints/boundaries
+ Generalization to other types besides `Float64`
