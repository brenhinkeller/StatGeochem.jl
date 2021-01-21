# StatGeochem.jl
> Some computational tools for geochemistry

[![Build Status](https://github.com/brenhinkeller/StatGeochem.jl/workflows/CI/badge.svg)](https://github.com/brenhinkeller/StatGeochem.jl/actions?query=workflow%3ACI)

## Installation

In the Julia package manager (type `]` in the REPL to enter)
```
pkg> add "https://github.com/brenhinkeller/StatGeochem.jl"
```

## Usage

This package can be used in the [Julia](https://julialang.org) REPL, in scripts or functions in Julia `.jl` files, in the [Juno/Atom IDE](http://junolab.org/), or in a Jupyter notebook. There aren't examples yet for most of the code in this repository, but for a quick demonstration, try the interactive Jupyter notebooks (it may take a few minutes for these to launch)
* [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/brenhinkeller/StatGeochem.jl/master?filepath=examples/BootstrapResamplingDemo.ipynb) Weighted bootstrap resampling: [BootstrapResamplingDemo.ipynb](https://mybinder.org/v2/gh/brenhinkeller/StatGeochem.jl/master?filepath=examples/BootstrapResamplingDemo.ipynb)
* [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/brenhinkeller/StatGeochem.jl/master?filepath=examples/ConstantSilicaReferenceModel.ipynb) Constant-silica reference model: [ConstantSilicaReferenceModel.ipynb](https://mybinder.org/v2/gh/brenhinkeller/StatGeochem.jl/master?filepath=examples/ConstantSilicaReferenceModel.ipynb)
* [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/brenhinkeller/StatGeochem.jl/master?filepath=examples/MeltsExamples.ipynb) Julia-alphaMELTS interface demo: [MeltsExamples.ipynb](https://mybinder.org/v2/gh/brenhinkeller/StatGeochem.jl/master?filepath=examples%2FMeltsExamples.ipynb)
* [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/brenhinkeller/StatGeochem.jl/master?filepath=examples/PerplexExamples.ipynb) Julia-Perple_X interface demo: [PerplexExamples.ipynb](https://mybinder.org/v2/gh/brenhinkeller/StatGeochem.jl/master?filepath=examples%2FPerplexExamples.ipynb)

The above links run notebooks from the [examples/](examples/) folder on a [JupyterHub](https://github.com/jupyterhub/jupyterhub) server hosted by the [Binder](https://mybinder.org) project. If you make changes to the online notebook, you can save them with `File` > `Download as` > `Notebook (.ipynb)` To run a downloaded notebook locally, use [IJulia](https://github.com/JuliaLang/IJulia.jl)

```Julia
julia> using IJulia
julia> notebook()
```
