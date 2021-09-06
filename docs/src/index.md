# Klimakoffer.jl

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
Then, obtain the Answer to the Ultimate Question of Life, The Universe, and Everything by executing
```julia
julia> answer()
42
```

## Authors
Klimakoffer is maintained by
[Gregor Gassner](https://www.mi.uni-koeln.de/NumSim/gassner),
[Johannes Markert](https://www.mi.uni-koeln.de/NumSim/markert),
[Andrés Rueda-Ramírez](https://www.mi.uni-koeln.de/NumSim/dr-andres-rueda-ramirez),
and
[Michael Schlottke-Lakemper](https://www.mi.uni-koeln.de/NumSim/schlottke-lakemper)
(all University of Cologne, Germany),
who are also the principal developers.

## License and contributing
Klimakoffer is licensed under the MIT license (see [License](@ref)).
