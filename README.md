# Klimakoffer.jl

[![Docs-dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://klimakoffer.github.io/Klimakoffer.jl/dev)
[![Build Status](https://github.com/klimakoffer/Klimakoffer.jl/workflows/CI/badge.svg)](https://github.com/klimakoffer/Klimakoffer.jl/actions?query=workflow%3ACI)
[![Coveralls](https://coveralls.io/repos/github/klimakoffer/Klimakoffer.jl/badge.svg?branch=main)](https://coveralls.io/github/klimakoffer/Klimakoffer.jl?branch=main)
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5221552.svg)](https://doi.org/10.5281/zenodo.5221552) -->

Klimakoffer is about making the process of climate simulation comprehensible and fun.

**Note: Klimakoffer is currently in pre-alpha stage and anything might change at any time.**

## Installation
If you have not yet installed Julia, please
[follow the instructions for your operating system](https://julialang.org/downloads/platform/).
Klimakoffer works with Julia v1.6.

To obtain Klimakoffer, clone the repository and use the Julia package manager
`Pkg` to install all dependencies:
```shell
git clone git@github.com:klimakoffer/Klimakoffer.jl.git
cd Klimakoffer.jl
julia --project=@. -e 'import Pkg; Pkg.instantiate()'
```
Then, start Julia with the `--project` flag set to your local clone, e.g.,
```shell
julia --project=.
```

## Usage
In the Julia REPL, first load the package Klimakoffer
```julia
julia> using Klimakoffer
```
Then, set the number of time steps per year for the solver
```julia
NT = 48 # this is a good default
```
Now you can create the mesh and model with
```julia
mesh = Mesh()
model = Model(mesh, NT)
```
and combine everything into the discretization:
```julia
discretization = Discretization(mesh, model, NT)
```
Finally, you can solve for the equilibrium temperature with
```julia
GlobTemp = compute_equilibrium!(discretization)
```

For the impatient, this example can also be easily reproduced by just including
the file `equilibrium_temperature_1950.jl` from the `examples` folder:
```julia
include(joinpath("examples", "equilibrium_temperature_1950.jl"))
```


## Documentation
There is not much there yet either, but
[here you go](https://klimakoffer.github.io/Klimakoffer.jl/dev).

## Authors
Klimakoffer is maintained by
[Gregor Gassner](https://www.mi.uni-koeln.de/NumSim/gassner),
[Johannes Markert](https://www.mi.uni-koeln.de/NumSim/markert),
[Christof Czernik](https://www.mi.uni-koeln.de/NumSim/christof-czernik),
[Andrés Rueda-Ramírez](https://www.mi.uni-koeln.de/NumSim/dr-andres-rueda-ramirez),
and
[Michael Schlottke-Lakemper](https://www.mi.uni-koeln.de/NumSim/schlottke-lakemper)
(all University of Cologne, Germany),
who are also the principal developers.

## License and contributing
Klimakoffer is licensed under the MIT license (see [LICENSE.md](LICENSE.md)).
