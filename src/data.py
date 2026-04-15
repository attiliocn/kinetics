"""Data loading, parsing, and filtering utilities."""

import pandas as pd
from . import config
from .kinetics import calculate_rate


def load_dataset(filepath, sheet_name='barriers'):
    """
    Load the reaction energy dataset from Excel.
    
    Args:
        filepath (str): Path to the Excel file
        sheet_name (str): Name of the sheet to load. Default 'barriers'.
    
    Returns:
        pd.DataFrame: Energy data with species as rows and reactions as columns
    """
    dataset = pd.read_excel(filepath, index_col=0, sheet_name=sheet_name)
    return dataset


def parse_reaction_id(reaction_label):
    """
    Parse a reaction ID into components.
    
    Args:
        reaction_label (str): Label like "im1-w1m"
    
    Returns:
        dict: Components {imido, roh}
    
    Example:
        >>> parse_reaction_id("im1-w1m")
        {'imido': 'im1', 'roh': 'w1m'}
    """
    parts = reaction_label.split('-')
    if len(parts) != 2:
        raise ValueError(f"Cannot parse reaction label: {reaction_label}")
    return {'imido': parts[0], 'roh': parts[1]}


def select_reactions(df, filters=None):
    """
    Select a subset of reactions based on filters.
    
    Args:
        df (pd.DataFrame): Full reaction dataset
        filters (dict): Boolean filters, e.g., {'imido': 'im1', 'roh': ['w1m', 'd1m']}
                       If None, returns all reactions.
    
    Returns:
        pd.DataFrame: Filtered subset
    """
    if not filters:
        return df
    
    selected = df.copy()
    
    # Parse reaction IDs to create filterable columns
    selected['imido'] = [parse_reaction_id(label)['imido'] for label in selected.index]
    selected['roh'] = [parse_reaction_id(label)['roh'] for label in selected.index]
    
    # Apply filters
    for column, value in filters.items():
        if isinstance(value, (list, tuple)):
            selected = selected.loc[selected[column].isin(value)]
        else:
            selected = selected.loc[selected[column] == value]
    
    # Drop the helper columns
    selected = selected.drop(columns=['imido', 'roh'], errors='ignore')
    
    return selected


def compute_rate_constants(barriers_df, temperature=None):
    """
    Compute rate constants from barrier energies.
    
    Applies calculate_rate() to all barrier columns and renames them
    according to the column mapping.
    
    Args:
        barriers_df (pd.DataFrame): DataFrame with barrier columns (b1, b2, b3, b4)
        temperature (float): Temperature in Celsius. If None, uses default.
    
    Returns:
        pd.DataFrame: DataFrame with rate constant columns (k1d, k2d, k1r, k2r)
    """
    if temperature is None:
        temperature = config.TEMPERATURE_DEFAULT
    
    # Apply calculate_rate to all rows
    rates_df = barriers_df.apply(lambda row: calculate_rate(row, temperature), axis=0)
    
    # Rename columns according to mapping
    rates_df = rates_df.rename(columns=config.BARRIER_COLUMN_MAPPING)
    
    # Select only the rate constants in the correct order
    rates_df = rates_df.loc[:, config.RATE_COLUMN_ORDER]
    
    return rates_df
