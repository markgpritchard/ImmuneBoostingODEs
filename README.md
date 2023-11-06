# ImmuneBoostingODEs

This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> ImmuneBoostingODEs

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

All scripts start with the commands:
```julia
using DrWatson
@quickactivate "ImmuneBoostingODEs"
```
which auto-activate the project and enable local path handling from DrWatson.

## Scripts 

* `figure1.jl` creates the model flowchart (Figure 1 in the manuscript)
* `immuneduration.jl` plots expected immune duration without immune boosting for supplementary figure
* `equilibriumvalues.jl` calculates equilibrium values and limit cycles
* `fouriertransforms.jl` perfoms simulations of long duration and calculates discrete Fourier transforms (included in supplementary figure in manuscript)
* `npisimulation.jl` runs simulations with a period of reduced transmission representing use of non-pharmaceutical interventions
* `rsvanalysis.jl` plots Scottish respiratory syncytial virus data and simulations with similar characteristics

## Data 

Processed data to reproduce all analysis are provided in `data/exp_pro`. If you prefer to use your own raw data, load the raw data into `exp_raw` and change the filename of the processed data files so they are not called in error.
