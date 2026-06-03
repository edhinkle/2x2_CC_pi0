""" 
Defines Matplotlib styles, colors, and fonts for 2x2 CC1pi0 analysis.
"""

"""
Styling and theme utilities for consistent plot appearance.

Provides configuration for colors, fonts, and matplotlib styles.
"""

from typing import Dict, Tuple, Optional, List
import matplotlib as mpl
import matplotlib.pyplot as plt


# Color palettes
COLORS = {
    # Primary analysis colors
    "reco": '#0072B2', # Blue
    "truth": '#F0E442', # Yellow
    "truthmatched": '#009E73', # Bluish Green
    
    # Signal vs Background
    "signal": '#D55E00', # Vermillion
    "background": '#0072B2', # Blue
    "data": "#000000",          # Black
    
    # Interaction modes
    "ccqe": '#E69F00', # Orange
    "ccmec": '#56B4E9', # Sky Blue
    "ccdis": '#009E73', # Bluish Green
    "ccres": '#F0E442', # Yellow
    "cccoh": '#CC79A7',  # Reddish Purple
    
    # Particle types
    "muon": '#009E73', # Bluish Green
    "electron":'#D55E00', # Vermillion
    "photon": '#F0E442', # Yellow
    "pion": '#CC79A7',  # Reddish Purple
    "proton": '#0072B2', # Blue
    
    # Multiplicity/Count colors
    "primary": '#E69F00', # Orange
    "secondary":  '#56B4E9', # Sky Blue

}

# Extended palette for categorical data
CATEGORICAL_COLORS = [
    '#E69F00', # Orange
    '#56B4E9', # Sky Blue
    '#009E73', # Bluish Green
    '#F0E442', # Yellow
    '#0072B2', # Blue
    '#D55E00', # Vermillion
    '#CC79A7',  # Reddish Purple
    '#999999',  # Grey
]

# Matplotlib style configuration
FONT_CONFIG = {
    "family": "sans-serif",
    "sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "size": 11,
}

TICK_CONFIG = {
    "labelsize": 10,
    "length": 5,
    "width": 1.2,
}

LEGEND_CONFIG = {
    "fontsize": 10,
    "framealpha": 0.95,
    "edgecolor": "black",
    "fancybox": True,
    "shadow": True,
}

LINE_CONFIG = {
    "linewidth": 1.5,
    "markersize": 6,
}

GRID_CONFIG = {
    "visible": True,
    "which": "major",
    "alpha": 0.3,
    "linestyle": "-",
    "linewidth": 0.8,
}


def apply_theme(theme: str = "default") -> None:
    """
    Apply a matplotlib theme globally.
    
    Args:
        theme: Theme name ("default", "dark", "publication")
    """
    if theme == "default":
        _apply_default_theme()
    elif theme == "dark":
        _apply_dark_theme()
    elif theme == "publication":
        _apply_publication_theme()
    else:
        raise ValueError(f"Unknown theme: {theme}")


def _apply_default_theme() -> None:
    """Apply default analysis theme."""
    plt.style.use("seaborn-v0_8-darkgrid")
    
    mpl.rcParams.update({
        # Font
        "font.family": FONT_CONFIG["family"],
        "font.size": FONT_CONFIG["size"],
        
        # Axes
        "axes.labelsize": 12,
        "axes.titlesize": 13,
        "axes.grid": True,
        "axes.axisbelow": True,
        
        # Ticks
        "xtick.labelsize": TICK_CONFIG["labelsize"],
        "ytick.labelsize": TICK_CONFIG["labelsize"],
        "xtick.major.size": TICK_CONFIG["length"],
        "ytick.major.size": TICK_CONFIG["length"],
        "xtick.major.width": TICK_CONFIG["width"],
        "ytick.major.width": TICK_CONFIG["width"],
        
        # Legend
        "legend.fontsize": LEGEND_CONFIG["fontsize"],
        "legend.framealpha": LEGEND_CONFIG["framealpha"],
        "legend.edgecolor": LEGEND_CONFIG["edgecolor"],
        
        # Grid
        "grid.alpha": GRID_CONFIG["alpha"],
        "grid.linestyle": GRID_CONFIG["linestyle"],
        "grid.linewidth": GRID_CONFIG["linewidth"],
        
        # Figure
        "figure.figsize": (10, 6),
        "figure.dpi": 100,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
    })


def _apply_dark_theme() -> None:
    """Apply dark theme for presentations."""
    plt.style.use("dark_background")
    
    mpl.rcParams.update({
        # Font
        "font.family": FONT_CONFIG["family"],
        "font.size": FONT_CONFIG["size"],
        
        # Axes
        "axes.labelsize": 12,
        "axes.titlesize": 13,
        "axes.facecolor": "#1a1a1a",
        
        # Ticks
        "xtick.labelsize": TICK_CONFIG["labelsize"],
        "ytick.labelsize": TICK_CONFIG["labelsize"],
        "xtick.color": "white",
        "ytick.color": "white",
        
        # Legend
        "legend.fontsize": LEGEND_CONFIG["fontsize"],
        "legend.framealpha": 0.9,
        "legend.facecolor": "#2a2a2a",
        
        # Grid
        "grid.alpha": 0.15,
        "grid.linestyle": "--",
        
        # Figure
        "figure.figsize": (10, 6),
        "figure.facecolor": "#1a1a1a",
        "savefig.facecolor": "#1a1a1a",
    })


def _apply_publication_theme() -> None:
    """Apply publication-quality theme."""
    mpl.rcParams.update({
        # Font - serif for publications
        "font.family": "serif",
        "font.serif": ["Times", "Palatino", "Computer Modern"],
        "font.size": 10,
        
        # Axes
        "axes.labelsize": 11,
        "axes.titlesize": 12,
        "axes.grid": False,
        "axes.spines.left": True,
        "axes.spines.bottom": True,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 1.0,
        
        # Ticks
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "xtick.major.size": 4,
        "ytick.major.size": 4,
        "xtick.major.width": 1.0,
        "ytick.major.width": 1.0,
        "xtick.direction": "in",
        "ytick.direction": "in",
        
        # Legend
        "legend.fontsize": 9,
        "legend.frameon": True,
        "legend.framealpha": 1.0,
        "legend.edgecolor": "black",
        "legend.fancybox": False,
        
        # Figure
        "figure.figsize": (8, 6),
        "figure.dpi": 100,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.1,
        
        # Lines
        "lines.linewidth": 1.5,
        "lines.markersize": 5,
    })


def get_colors(n_colors: Optional[int] = None, 
               category: Optional[str] = None) -> Dict[str, str] | List[str]:
    """
    Get color palette.
    
    Args:
        n_colors: If specified, return first n colors from categorical palette
        category: If specified, return colors dict for that category
                 (e.g., "interaction_modes", "particle_types")
    
    Returns:
        Dictionary of colors, or list of n colors
    """
    if category == "interaction_modes":
        return {
            "CCQE": COLORS["ccqe"],
            "CCMEC": COLORS["ccmec"],
            "CCDIS": COLORS["ccdis"],
            "CCRES": COLORS["ccres"],
            "CCCOH": COLORS["cccoh"],
        }
    elif category == "particle_types":
        return {
            "muon": COLORS["muon"],
            "electron": COLORS["electron"],
            "photon": COLORS["photon"],
            "pion": COLORS["pion"],
            "proton": COLORS["proton"],
        }
    elif category == "reco_truth_matched":
        return {
            "reco": COLORS["reco"],
            "truth": COLORS["truth"],
            "truthmatched": COLORS["truthmatched"],
        }
    elif n_colors is not None:
        return CATEGORICAL_COLORS[:n_colors]
    else:
        return COLORS


def get_color(key: str, default: str = "#000000") -> str:
    """
    Get a single color by key.
    
    Args:
        key: Color key (e.g., "reco", "truth", "signal")
        default: Default color if key not found
    
    Returns:
        Hex color string
    """
    return COLORS.get(key, default)


def get_font_config() -> Dict:
    """
    Get font configuration dictionary.
    
    Returns:
        Font config dict suitable for plt.rcParams
    """
    return {
        "font.family": FONT_CONFIG["family"],
        "font.sans-serif": FONT_CONFIG["sans-serif"],
        "font.size": FONT_CONFIG["size"],
    }


def get_tick_config() -> Dict:
    """
    Get tick configuration dictionary.
    
    Returns:
        Tick config dict suitable for plt.rcParams
    """
    return {
        "xtick.labelsize": TICK_CONFIG["labelsize"],
        "ytick.labelsize": TICK_CONFIG["labelsize"],
        "xtick.major.size": TICK_CONFIG["length"],
        "ytick.major.size": TICK_CONFIG["length"],
        "xtick.major.width": TICK_CONFIG["width"],
        "ytick.major.width": TICK_CONFIG["width"],
    }


def get_legend_config() -> Dict:
    """
    Get legend configuration dictionary.
    
    Returns:
        Legend config dict suitable for ax.legend()
    """
    return {
        "fontsize": LEGEND_CONFIG["fontsize"],
        "framealpha": LEGEND_CONFIG["framealpha"],
        "edgecolor": LEGEND_CONFIG["edgecolor"],
        "fancybox": LEGEND_CONFIG["fancybox"],
        "shadow": LEGEND_CONFIG["shadow"],
    }


def get_grid_config() -> Dict:
    """
    Get grid configuration dictionary.
    
    Returns:
        Grid config dict suitable for ax.grid()
    """
    return {
        "visible": GRID_CONFIG["visible"],
        "which": GRID_CONFIG["which"],
        "alpha": GRID_CONFIG["alpha"],
        "linestyle": GRID_CONFIG["linestyle"],
        "linewidth": GRID_CONFIG["linewidth"],
    }


def configure_axes(ax, 
                   xlabel: Optional[str] = None,
                   ylabel: Optional[str] = None,
                   title: Optional[str] = None,
                   grid: bool = True,
                   legend: bool = False) -> None:
    """
    Configure axis labels, title, and grid for a single axes.
    
    Args:
        ax: Matplotlib axes object
        xlabel: X-axis label
        ylabel: Y-axis label
        title: Plot title
        grid: Whether to show grid
        legend: Whether to show legend
    """
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    
    if grid:
        ax.grid(**get_grid_config())
    else:
        ax.grid(False)
    
    if legend:
        ax.legend(**get_legend_config())


def set_color_cycle(n_colors: int) -> None:
    """
    Set matplotlib color cycle to use first n colors from categorical palette.
    
    Args:
        n_colors: Number of colors to use
    """
    colors = CATEGORICAL_COLORS[:n_colors]
    mpl.rcParams['axes.prop_cycle'] = mpl.cycler('color', colors)
