"""ODE system definition for the reaction mechanism."""

from . import config


def kinetic_equations(t, c, *rates):
    """
    Define the ODEs for the reaction mechanism kinetics.
    
    Mass-action kinetics for a multi-step reaction network:
        Step 1: bispyr + roh ⇌ map + me2pyr  (k1d, k1r)
        Step 2: map + roh ⇌ bis + me2pyr     (k2d, k2r)
    
    Args:
        t (float): Time (unused, but required for scipy.integrate.solve_ivp)
        c (list or array): Concentration vector [bispyr, roh, me2pyr, map, bis]
        *rates: Variable rate constants (k1d, k1r, k2d, k2r)
    
    Returns:
        list: Time derivatives [d[bispyr]/dt, d[roh]/dt, d[me2pyr]/dt, d[map]/dt, d[bis]/dt]
    """
    # Convert concentration list to dict for readability
    entities = config.ENTITIES  # ['bispyr', 'roh', 'me2pyr', 'map', 'bis']
    conc_dict = {k: v for k, v in zip(entities, c)}
    
    # Unpack rate constants
    k1d = rates[0]
    k1r = rates[1]
    k2d = rates[2]
    k2r = rates[3]
    
    # Define ODEs based on mechanism
    # Step 1: bispyr + roh <-> map + me2pyr
    # Step 2: map + roh <-> bis + me2pyr
    
    d_bispyr_dt = (
        - (k1d * conc_dict['bispyr'] * conc_dict['roh'])
        + (k1r * conc_dict['map'] * conc_dict['me2pyr'])
    )
    
    d_roh_dt = (
        - (k1d * conc_dict['bispyr'] * conc_dict['roh'])
        + (k1r * conc_dict['map'] * conc_dict['me2pyr'])
        - (k2d * conc_dict['map'] * conc_dict['roh'])
        + (k2r * conc_dict['bis'] * conc_dict['me2pyr'])
    )
    
    d_me2pyr_dt = (
        + (k1d * conc_dict['bispyr'] * conc_dict['roh'])
        - (k1r * conc_dict['map'] * conc_dict['me2pyr'])
        + (k2d * conc_dict['map'] * conc_dict['roh'])
        - (k2r * conc_dict['bis'] * conc_dict['me2pyr'])
    )
    
    d_map_dt = (
        + (k1d * conc_dict['bispyr'] * conc_dict['roh'])
        - (k1r * conc_dict['map'] * conc_dict['me2pyr'])
        - (k2d * conc_dict['map'] * conc_dict['roh'])
        + (k2r * conc_dict['bis'] * conc_dict['me2pyr'])
    )
    
    d_bis_dt = (
        + (k2d * conc_dict['map'] * conc_dict['roh'])
        - (k2r * conc_dict['bis'] * conc_dict['me2pyr'])
    )
    
    return [d_bispyr_dt, d_roh_dt, d_me2pyr_dt, d_map_dt, d_bis_dt]


def concentrations_list_to_dict(concentrations_list):
    """Convert concentration list to dictionary."""
    return {k: v for k, v in zip(config.ENTITIES, concentrations_list)}


def concentrations_dict_to_list(concentrations_dict):
    """Convert concentration dictionary to list in entity order."""
    return [concentrations_dict[entity] for entity in config.ENTITIES]
