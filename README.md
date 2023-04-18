# StatGeochem.jl
[![DOI](osf_io_TJHMW.svg)](https://doi.org/10.17605/OSF.IO/TJHMW)
[![Docs][docs-dev-img]][docs-dev-url]
[![CI][ci-img]][ci-url]
[![CI (Julia nightly)][ci-nightly-img]][ci-nightly-url]
[![Coverage][codecov-img]][codecov-url]

_Computational tools for statistical geochemistry and petrology_

In addition to functions exported by StatGeochem directly, StatGeochem also reexports (and depends upon internally) both [StatGeochemBase.jl](https://github.com/brenhinkeller/StatGeochemBase.jl) and [NaNStatistics.jl](https://github.com/brenhinkeller/NaNStatistics.jl)

## Installation

StatGeochem.jl is written in the [Julia programming language](https://julialang.org/), and is registered on the General registry. To install, enter the Julia package manager (type `]` in the REPL) and type:
```
pkg> add StatGeochem
```
If you are trying to use a script written prior to ~2021, you may want to use the oldest registered version of the package, which you can install with (e.g.)
```Julia
julia> add StatGeochem@v0.1
```

## Usage

This package can be used in the [Julia](https://julialang.org) REPL, in scripts or functions in Julia `.jl` files, in the [Juno/Atom IDE](http://junolab.org/), or in a Jupyter notebook. There aren't examples yet for most of the code in this repository, but for a quick demonstration, try the interactive Jupyter notebooks (it may take a few minutes for these to launch)
* [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/brenhinkeller/StatGeochem.jl/main?filepath=examples/BootstrapResamplingDemo.ipynb) Weighted bootstrap resampling: [BootstrapResamplingDemo.ipynb](https://mybinder.org/v2/gh/brenhinkeller/StatGeochem.jl/main?filepath=examples/BootstrapResamplingDemo.ipynb)
* [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/brenhinkeller/StatGeochem.jl/main?filepath=examples/ConstantSilicaReferenceModel.ipynb) Constant-silica reference model: [ConstantSilicaReferenceModel.ipynb](https://mybinder.org/v2/gh/brenhinkeller/StatGeochem.jl/main?filepath=examples/ConstantSilicaReferenceModel.ipynb)
* [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/brenhinkeller/StatGeochem.jl/main?filepath=examples/MeltsExamples.ipynb) Julia-alphaMELTS interface demo: [MeltsExamples.ipynb](https://mybinder.org/v2/gh/brenhinkeller/StatGeochem.jl/main?filepath=examples%2FMeltsExamples.ipynb)
* [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/brenhinkeller/StatGeochem.jl/main?filepath=examples/PerplexExamples.ipynb) Julia-Perple_X interface demo: [PerplexExamples.ipynb](https://mybinder.org/v2/gh/brenhinkeller/StatGeochem.jl/main?filepath=examples%2FPerplexExamples.ipynb)

The above links run notebooks from the [examples/](examples/) folder on a [JupyterHub](https://github.com/jupyterhub/jupyterhub) server hosted by the [Binder](https://mybinder.org) project. If you make changes to the online notebook, you can save them with `File` > `Download as` > `Notebook (.ipynb)` To run a downloaded notebook locally, use [IJulia](https://github.com/JuliaLang/IJulia.jl)

```Julia
julia> using IJulia
julia> notebook()
```

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://brenhinkeller.github.io/StatGeochem.jl/dev/
[ci-img]: https://github.com/brenhinkeller/StatGeochem.jl/workflows/CI/badge.svg
[ci-url]: https://github.com/brenhinkeller/StatGeochem.jl/actions/workflows/CI.yml
[ci-nightly-img]:https://github.com/brenhinkeller/StatGeochem.jl/workflows/CI%20(Julia%20nightly)/badge.svg
[ci-nightly-url]:https://github.com/brenhinkeller/StatGeochem.jl/actions/workflows/CI-julia-nightly.yml
[codecov-img]: http://codecov.io/github/brenhinkeller/StatGeochem.jl/coverage.svg?branch=main
[codecov-url]: http://app.codecov.io/github/brenhinkeller/StatGeochem.jl?branch=main
