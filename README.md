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
julia> NT = 48 # this is a good default
48
```
Now you can create the mesh and model with
```julia
julia> mesh = Mesh()
Mesh() with 128×65 degrees of freedom

julia> model = Model(mesh, NT)
Model() with 128×65 degrees of freedom
```
and combine everything into the discretization:
```julia
julia> discretization = Discretization(mesh, model, NT)
Discretization() with 128×65 degrees of freedom
```
Finally, you can solve for the equilibrium temperature with
```julia
julia> GlobTemp = compute_equilibrium!(discretization)
year  Average Temperature
0  5.000000000000189
1  9.004451135208686
2  9.79480397404041
3  10.891438923925572
4  11.790964037289493
5  12.481170727332831
6  12.99908587812284
7  13.384303292306818
8  13.669812313384279
9  13.88115978138856
10  14.037575162852074
11  14.153361080313475
12  14.239103641567818
13  14.30262469707876
14  14.349702021948694
15  14.384604941755951
16  14.410489957304135
17  14.429692258370507
18  14.443940447343182
19  14.454514784995716
20  14.462363914213213
21  14.468191036721114
22  14.472517593123692
23  14.47573035613406
24  14.478116281193246
25  14.479888312227843
26  14.481204499489236
27  14.482182168723126
28  14.482908426404402
29  14.483447950855698
30  14.483848771972168
31  14.484146559783808
32  14.484367807186084
33  14.484532192378012
34  14.484654332676165
35  14.48474508666562
36  14.484812521111287
37  14.484862629006399
38  14.484899862855066
39  14.484927530761858
40  14.484948090641815
41  14.484963368768746
EQUILIBRIUM REACHED!
14.484963368768746
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
