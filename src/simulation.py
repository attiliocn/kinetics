"""Numerical integration and trajectory utilities."""

import numpy as np
from scipy.integrate import solve_ivp

from . import config
from .ode_system import kinetic_equations


def run_simulation(rates, concentrations, times):
    """
    Run the ODE simulation for the reaction mechanism.
    
    Integrates the kinetic ODEs over the specified time points using the
    configured ODE solver (BDF with tight tolerances for stiff systems).
    
    Args:
        rates (dict): Rate constants {k1d, k1r, k2d, k2r}
        concentrations (dict): Initial concentrations {bispyr, roh, ...}
        times (array-like): Time points to evaluate the solution
    
    Returns:
        scipy.integrate.OdeResult: Solution object with attributes:
            - t: Time points
            - y: Concentration profiles (5 x len(times))
            - success: Boolean flag
            - message: Status message
    """
    # Convert concentrations dict to list in ENTITIES order
    c0 = [concentrations[entity] for entity in config.ENTITIES]
    
    # Convert rates dict to tuple in expected order
    rates_tuple = tuple(rates[key] for key in config.RATE_COLUMN_ORDER)
    
    # Run the integration
    simulation = solve_ivp(
        fun=kinetic_equations,
        t_span=[0, times[-1]],
        y0=c0,
        args=rates_tuple,
        t_eval=times,
        method=config.ODE_SOLVER_OPTIONS['method'],
        rtol=config.ODE_SOLVER_OPTIONS['rtol'],
        atol=config.ODE_SOLVER_OPTIONS['atol'],
    )
    
    return simulation


def exponential_trajectory_points(x_end, num_points=100, exponent=3):
    """
    Generate time points clustered more densely at the beginning.
    
    Useful for capturing fast initial transients in kinetic simulations
    while still covering the full time span.
    
    Args:
        x_end (float): Final time value (starts from 0)
        num_points (int): Total number of points to generate. Default 100.
        exponent (float): Controls clustering: >1 compresses early times,
                          <1 compresses late times. Default 3.
    
    Returns:
        np.ndarray: Array of num_points time values in [0, x_end]
    
    Example:
        >>> points = exponential_trajectory_points(3600, num_points=500, exponent=3)
        >>> # Returns 500 points in [0, 3600], denser at t=0
    """
    s = np.linspace(0, 1, num_points)
    stretched = s ** exponent
    return stretched * x_end
