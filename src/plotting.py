"""Visualization utilities for kinetic simulations."""

import matplotlib.pyplot as plt
import seaborn as sns


def plot_kinetics_roh(times, concentrations, reaction_label, colors=None, 
                      set_ylim=False, set_xlim=False, save_path=None):
    """
    Plot concentration vs. time for the reaction species.
    
    Creates a figure with concentration profiles for bispyr, roh, map, and bis.
    
    Args:
        times (np.ndarray): Time array in seconds
        concentrations (dict): Species concentrations {entity_name: concentration_array}
        reaction_label (str): Label for the reaction (used in title and filename)
        colors (list or palette): Color palette for the species. If None, uses muted.
        set_ylim (bool): Whether to set fixed y-axis limits. Default False.
        set_xlim (bool): Whether to set fixed x-axis limits. Default False.
        save_path (str): Directory to save plots. If None, plots are not saved.
    
    Returns:
        tuple: (fig, ax)
    """
    if colors is None:
        colors = sns.color_palette('muted')
    
    plt.close('all')
    fig, ax = plt.subplots(figsize=(8, 5))
    
    # Convert times from seconds to minutes
    times_min = times / 60
    
    # Bispyrrolide (initial reactant)
    ax.plot(
        times_min,
        concentrations['bispyr'] * 1e3,
        color=colors[0],
        label='Bispyrrolide',
        linestyle='-',
        linewidth=2,
    )
    
    # ROH (reactant)
    ax.plot(
        times_min,
        concentrations['roh'] * 1e3,
        color=colors[3],
        label='ROH',
        linestyle='--',
        linewidth=2,
    )
    
    # MAP (intermediate product)
    ax.plot(
        times_min,
        concentrations['map'] * 1e3,
        color=colors[1],
        label='MAP',
        linestyle=':',
        linewidth=2,
    )
    
    # BIS (final product)
    ax.plot(
        times_min,
        concentrations['bis'] * 1e3,
        color=colors[2],
        label='BIS',
        linestyle='-.',
        linewidth=2,
    )
    
    # Axis labels and title
    ax.set_xlabel('Time (min)', fontsize=11)
    ax.set_ylabel('Concentration (mmol/L)', fontsize=11)
    ax.set_title(f'{reaction_label}', fontsize=12, fontweight='bold')
    
    # Optional axis limits
    if set_ylim:
        ax.set_ylim(-0.05, 2.5)
    if set_xlim:
        ax.set_xlim(-5, 310)
    
    ax.legend(loc='best', fontsize=10)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    
    # Optional saving
    if save_path:
        fig.savefig(f'{save_path}/{reaction_label}.svg', dpi=150, bbox_inches='tight')
        fig.savefig(f'{save_path}/{reaction_label}.png', dpi=150, bbox_inches='tight')
    
    return fig, ax
