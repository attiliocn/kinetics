"""Kinetic theory calculations for barrier energies to rate constants."""

import numpy as np


def calculate_rate(reaction_barrier, temperature):
    """
    Converts a barrier energy to a reaction rate constant.
    
    Uses transition state theory (Eyring-Polanyi equation variant).
    
    Args:
        reaction_barrier (float): Barrier energy in kcal/mol
        temperature (float): Temperature in Celsius
    
    Returns:
        float: Rate constant in s^-1
    
    Notes:
        - Uses fundamental physical constants (Boltzmann, Planck)
        - Temperature is internally converted from Celsius to Kelvin
    """
    R = 8.314  # Gas constant (J/mol·K)
    T = temperature + 273.15  # Convert Celsius to Kelvin
    RT = (R * T) / 1000  # Convert to kJ/mol
    
    kB = 1.380649e-23  # Boltzmann constant (J/K)
    h = 6.62607015e-34  # Planck constant (J·s)
    kBTh = kB * T / h  # Pre-exponential factor (1/s)
    
    # Exponent uses barrier in kcal/mol, conversion factor 4.184 kJ/kcal
    rate = kBTh * np.exp(-reaction_barrier / (RT / 4.184))
    return rate
