
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Days before first day of each month 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# NB assuming 28 days in February
const MONTHDAYS = [ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 ]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# const used in equilibria.jl 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

const TRIALPSIS = [ 10^i for i ∈ collect(0:1:18) ]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Order of model compartments 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

const COMPARTMENTINDICES = Dict(:S => 1, :I => 2, :R1 => 3, :R2 => 4, :R3 => 5, :cc => 8)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Constants used in plotting 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Consistent colour scheme across plots 

const COLOUR_S = :blue
const COLOUR_I = :darkgoldenrod1
const COLOUR_R = :seagreen4
const COLOURVECTOR = [ 
    COLOUR_S, COLOUR_I, COLOUR_R, :plum, :brown2, :dodgerblue3, :skyblue2, :lightgray 
]

# outputs from MCMC that are not plotted 
const _NOPLOTNAMES = [ 
    "iteration", "chain", "lp", "n_steps", "is_accept", "acceptance_rate", "log_density", 
    "hamiltonian_energy", "hamiltonian_energy_error", "max_hamiltonian_energy_error", 
    "tree_depth", "numerical_error", "step_size", "nom_step_size"
]
