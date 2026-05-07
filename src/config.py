"""Configuration and constants for the microkinetics system."""

# Reaction network entities
ENTITIES = ['bispyr', 'roh', 'me2pyr', 'map', 'bis']

# Species classification (Phase 0 decisions)
# Intermediates: perturbed by A only
INTERMEDIATES = ['map', 'bis', 'a1', 'a2']

# Transition states: perturbed by B = A * x
TRANSITION_STATES = ['ts1', 'ts2']
TRANSITION_STATES_LABELS = TRANSITION_STATES

# Default initial concentrations (in mol/L)                          <<<< CHANGE HERE IF NECESSARY >>>>
INITIAL_CONCENTRATIONS = {
    'bispyr': 0.1,
    'roh': 0.2,
    'me2pyr': 0.0,
    'map': 0.0,
    'bis': 0.0,
}

# ODE solver options
ODE_SOLVER_OPTIONS = {
    'method': 'BDF',
    'rtol': 1e-12,
    'atol': 1e-12,
}

# Perturbation mode defaults                                         <<<< CHANGE HERE IF NECESSARY >>>>
PERTURBATION_SIGMA_DEFAULT = 3.0
PERTURBATION_SEED_DEFAULT = 42
PERTURBATION_N_SAMPLES_DEFAULT = 1000
PERTURBATION_PERCENTILES_DEFAULT = (5, 25, 50, 75, 95)

# Barrier to rate constant calculation                               <<<< CHANGE HERE IF NECESSARY >>>>
# Temperature in Celsius
TEMPERATURE_DEFAULT = 70

# Rate threshold for accessibility check (s^-1)
RATE_ACCESSIBILITY_THRESHOLD = 1e-10

# Simulation defaults                                                <<<< CHANGE HERE IF NECESSARY >>>>
SIMULATION_TIME_DEFAULT = 3600 * 12  # n hours in seconds
TRAJECTORY_POINTS_DEFAULT = 1000
TRAJECTORY_EXPONENT_DEFAULT = 3

# Column mappings for energy to rate conversion
BARRIER_COLUMN_MAPPING = {
    'b1': 'k1d',  # barrier 1 direct -> k1d
    'b2': 'k2d',  # barrier 2 direct -> k2d
    'b3': 'k1r',  # barrier 3 reverse -> k1r
    'b4': 'k2r',  # barrier 4 reverse -> k2r
}

# Rate column order (matches entity order for ODE system)
RATE_COLUMN_ORDER = ['k1d', 'k1r', 'k2d', 'k2r']
