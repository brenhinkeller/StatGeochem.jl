[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/brenhinkeller/StatGeochem.jl/master?filepath=examples/BootstrapResamplingDemo.ipynb)

# StatGeochem.jl
Some computational tools for geochemistry

## Installation

In the Julia package manager (type `]` in the REPL to enter)
```
(v1.0) pkg> add "https://github.com/brenhinkeller/StatGeochem.jl"
```
or for previous versions of Julia, in the REPL
```
julia> Pkg.clone("https://github.com/brenhinkeller/StatGeochem.jl")
```

## Usage

This package can be used in the [Julia](https://julialang.org) REPL, in scripts or functions in Julia `.jl` files, in the [Juno/Atom IDE](http://junolab.org/), or in a Jupyter notebook. There aren't examples yet for most of the code in this repository, but for a quick demonstration, try the interactive Jupyter notebooks (it may take a few minutes for these to launch)
* Weighted bootstrap resampling: [BootstrapResamplingDemo.ipynb](https://mybinder.org/v2/gh/brenhinkeller/StatGeochem.jl/master?filepath=examples/BootstrapResamplingDemo.ipynb)
* Julia-alphaMELTS interface demo: [MeltsExamples.ipynb](https://mybinder.org/v2/gh/brenhinkeller/StatGeochem.jl/master?filepath=examples%2FMeltsExamples.ipynb)


The above links run notebooks from the [examples/](examples/) folder on a [JupyterHub](https://github.com/jupyterhub/jupyterhub) server hosted by the [Binder](https://mybinder.org) project. If you make changes to the online notebook, you can save them with `File` > `Download as` > `Notebook (.ipynb)` To run a downloaded notebook locally, use [IJulia](https://github.com/JuliaLang/IJulia.jl)

```Julia
julia> using IJulia
julia> notebook()
```
