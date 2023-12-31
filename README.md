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

* `immuneduration.jl` calculates expected immune duration without immune boosting for supplementary figure
* `equilibriumvalues.jl` calculates equilibrium values and limit cycles
* `fouriertransforms.jl` perfoms simulations of long duration and calculates discrete Fourier transforms (included in supplementary figure in manuscript)
* `npisimulation.jl` runs simulations with a period of reduced transmission representing use of non-pharmaceutical interventions
* `rsvanalysis.jl` calculates Scottish respiratory syncytial virus data and simulations with similar characteristics
* `figuresformanuscript.jl` produces figures formatted for a manuscript 
* `figuresforposter.jl` produces figures formatted for a poster

## Data 

Processed data to reproduce all analysis are provided in `data/exp_pro`. If you prefer to use your own raw data, load the raw data into `exp_raw` and change the filename of the processed data files so they are not called in error.

Data used in this analysis were released by Public Health Scotland under the [Open Government Licence for public sector information](https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/).

## Funding 

This study is funded by the National Institute for Health Research (NIHR) Health Protection Research Unit in Emerging and Zoonotic Infections (200907), a partnership between the United Kingdom Health Security Agency (UKHSA), the University of Liverpool, the University of Oxford and the Liverpool School of Tropical Medicine. The views expressed are those of the authors and not necessarily those of the NIHR, the UKHSA or the Department of Health and Social Care.

