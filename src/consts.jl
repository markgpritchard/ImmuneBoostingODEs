
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Days before first day of each month 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# NB assuming 28 days in February
const MONTHDAYS = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# const used in equilibria.jl 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

const TRIALPSIS = [10^i for i ∈ collect(0:1:18)]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Order of model compartments 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

const COMPARTMENTINDICES = Dict(:S => 1, :I => 2, :R1 => 3, :R2 => 4, :R3 => 5, :cc => 8)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Constants used in plotting 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Consistent colour scheme across plots 

const COLOURVECTOR = [
    RGBf(33 / 255, 145 / 255, 140 / 255),
    RGBf(68 / 255, 57 / 255, 131 / 255),
    RGBf(253 / 255, 231 / 255, 37 / 255),
    RGBf(53 / 255, 183 / 255, 121 / 255),
    RGBf(49 / 255, 104 / 255, 142 / 255),
    RGBf(68 / 255, 1 / 255, 84 / 255),
    RGBf(144 / 255, 215 / 255, 67 / 255),
]

const COLOUR_S = COLOURVECTOR[1]
const COLOUR_I = COLOURVECTOR[2]
const COLOUR_R = COLOURVECTOR[3]

# outputs from MCMC that are not plotted 
const _NOPLOTNAMES = [ 
    "iteration", 
    "chain", 
    "lp", 
    "n_steps", 
    "is_accept", 
    "acceptance_rate", 
    "log_density", 
    "hamiltonian_energy", 
    "hamiltonian_energy_error", 
    "max_hamiltonian_energy_error", 
    "tree_depth", 
    "numerical_error", 
    "step_size", 
    "nom_step_size",
]
