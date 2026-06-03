"""
Formatting utilities for histogram plots.

Provides helper functions for common formatting tasks like axis labels,
legends, colorbars, and other plot decorations.
"""

from typing import Optional, List, Tuple, Dict, Any
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import FuncFormatter, MaxNLocator
import numpy as np


# Common axis labels for physics quantities
PHYSICS_LABELS = {
    "cosL": r"$\cos(\theta)$",
    "cos_theta": r"$\cos(\theta)$",
    "cos_angle": r"$\cos(\theta)$",
    "energy": "Energy [MeV]",
    "E": "Energy [MeV]",
    "frac": "Fraction of Events",
    "momentum": "Momentum [MeV/c]",
    "p": "Momentum [MeV/c]",
    "angle": r"Angle [$^\circ$]",
    "theta": r"$\theta$ [$^\circ$]",
    "phi": r"$\phi$ [$^\circ$]",
    "mass": "Invariant Mass [MeV/c$^2$]",
    "multiplicity": "Multiplicity",
    "count": "Count",
    "entries": "Number of Events",
    "vertex_x": "Vertex X [cm]",
    "vertex_y": "Vertex Y [cm]",
    "vertex_z": "Vertex Z [cm]",
    "position_x": "X Position [cm]",
    "position_y": "Y Position [cm]",
    "position_z": "Z Position [cm]",
    "pot": "POT",
    "pot_count": "Protons on Target",
}

# Common histogram titles/descriptions
HIST_DESCRIPTIONS = {
    "reco": "Reconstructed",
    "truth": "Truth",
    "matched": "Truth-Matched",
    "signal": "Signal",
    "background": "Background",
    "data": "Data",
    "mc": "Monte Carlo",
}


def get_axis_label(key: str, unit: Optional[str] = None, 
                   default: Optional[str] = None) -> str:
    """
    Get formatted axis label for common physics quantities.
    
    Args:
        key: Quantity key (e.g., "cosL", "energy", "multiplicity")
        unit: Optional unit override
        default: Default label if key not found
    
    Returns:
        Formatted axis label string
        
    Example:
        >>> get_axis_label("cosL")
        '$\\cos(\\theta)$'
        >>> get_axis_label("energy", unit="GeV")
        'Energy (GeV)'
    """
    label = PHYSICS_LABELS.get(key, default or key)
    
    if unit and "[" not in label:
        label = f"{label} [{unit}]"
    
    return label


def get_hist_title(hist_name: str, category: Optional[str] = None) -> str:
    """
    Generate a formatted title for a histogram.
    
    Args:
        hist_name: Histogram name (e.g., "h_reco_VertexX")
        category: Optional category ("reco", "truth", "matched", etc.)
    
    Returns:
        Formatted title string
        
    Example:
        >>> get_hist_title("h_reco_VertexX", category="reco")
        'Reconstructed Vertex X'
    """
    # Remove 'h_' prefix and category prefix
    title = hist_name.replace("h_", "")
    
    # Try to extract parts
    parts = title.split("_")
    
    # If category is specified, use it
    if category and category in HIST_DESCRIPTIONS:
        prefix = HIST_DESCRIPTIONS[category]
        # Remove category from parts if present
        if parts and parts[0].lower() in HIST_DESCRIPTIONS:
            parts = parts[1:]
        title = " ".join(parts)
        return f"{prefix} {title}"
    
    # Convert camelCase or snake_case to Title Case
    title = title.replace("_", " ")
    title = title.replace("Reco", "Reconstructed")
    title = title.replace("Truth", "Truth")
    title = title.replace("truthMatch", "Truth Matched")
    
    return title.title()


def setup_legend(ax, loc: str = "best", ncol: int = 1,
                framealpha: float = 0.95, fontsize: Optional[int] = None,
                **kwargs) -> None:
    """
    Configure and display legend with consistent styling.
    
    Args:
        ax: Matplotlib axes object
        loc: Legend location
        ncol: Number of columns
        framealpha: Frame transparency
        fontsize: Font size (uses theme default if None)
        **kwargs: Additional arguments passed to ax.legend()
    """
    legend_kwargs = {
        "loc": loc,
        "ncol": ncol,
        "framealpha": framealpha,
        "edgecolor": "black",
        "fancybox": True,
        "shadow": True,
    }
    
    if fontsize:
        legend_kwargs["fontsize"] = fontsize
    
    legend_kwargs.update(kwargs)
    ax.legend(**legend_kwargs)


def add_grid(ax, which: str = "major", visible: bool = True,
            alpha: float = 0.3, linestyle: str = "-",
            linewidth: float = 0.8, **kwargs) -> None:
    """
    Add grid to axes with consistent styling.
    
    Args:
        ax: Matplotlib axes object
        which: "major", "minor", or "both"
        visible: Whether grid is visible
        alpha: Grid transparency
        linestyle: Grid line style
        linewidth: Grid line width
        **kwargs: Additional arguments passed to ax.grid()
    """
    grid_kwargs = {
        "which": which,
        "alpha": alpha,
        "linestyle": linestyle,
        "linewidth": linewidth,
        "visible": visible,
    }
    
    grid_kwargs.update(kwargs)
    ax.grid(**grid_kwargs)


def set_axis_limits(ax, xlim: Optional[Tuple[float, float]] = None,
                   ylim: Optional[Tuple[float, float]] = None) -> None:
    """
    Set axis limits with validation.
    
    Args:
        ax: Matplotlib axes object
        xlim: X-axis limits (xmin, xmax)
        ylim: Y-axis limits (ymin, ymax)
        
    Raises:
        ValueError: If limits are invalid (e.g., xmin >= xmax)
    """
    if xlim is not None:
        if xlim[0] >= xlim[1]:
            raise ValueError(f"Invalid xlim: {xlim}")
        ax.set_xlim(xlim)
    
    if ylim is not None:
        if ylim[0] >= ylim[1]:
            raise ValueError(f"Invalid ylim: {ylim}")
        ax.set_ylim(ylim)


def format_axis_ticks(ax, axis: str = "both", 
                     format_str: Optional[str] = None,
                     n_ticks: Optional[int] = None,
                     rotation: int = 0) -> None:
    """
    Format tick labels on axis.
    
    Args:
        ax: Matplotlib axes object
        axis: "x", "y", or "both"
        format_str: Format string (e.g., "%.2f" for 2 decimals, "%.0e" for scientific)
        n_ticks: Approximate number of ticks to display
        rotation: Tick label rotation angle
        
    Example:
        >>> format_axis_ticks(ax, axis="y", format_str="%.2f", n_ticks=5)
    """
    if axis in ["x", "both"]:
        if n_ticks:
            ax.xaxis.set_major_locator(MaxNLocator(nbins=n_ticks))
        if format_str:
            ax.xaxis.set_major_formatter(FuncFormatter(lambda x, p: format_str % x))
        ax.tick_params(axis="x", rotation=rotation)
    
    if axis in ["y", "both"]:
        if n_ticks:
            ax.yaxis.set_major_locator(MaxNLocator(nbins=n_ticks))
        if format_str:
            ax.yaxis.set_major_formatter(FuncFormatter(lambda x, p: format_str % x))
        ax.tick_params(axis="y", rotation=rotation)


def format_scientific_notation(ax, axis: str = "both", 
                               offset_notation: bool = True) -> None:
    """
    Format axis ticks in scientific notation.
    
    Args:
        ax: Matplotlib axes object
        axis: "x", "y", or "both"
        offset_notation: Use offset notation (e.g., 1e5) vs full notation
    """
    from matplotlib.ticker import ScalarFormatter
    
    formatter = ScalarFormatter(useMathText=True)
    if not offset_notation:
        formatter.set_useOffset(False)
    
    if axis in ["x", "both"]:
        ax.xaxis.set_major_formatter(formatter)
    if axis in ["y", "both"]:
        ax.yaxis.set_major_formatter(formatter)


def format_log_scale(ax, axis: str = "y", base: int = 10) -> None:
    """
    Set logarithmic scale on axis.
    
    Args:
        ax: Matplotlib axes object
        axis: "x", "y", or "both"
        base: Logarithm base (10, 2, etc.)
    """
    if axis in ["x", "both"]:
        ax.set_xscale("log", base=base)
    if axis in ["y", "both"]:
        ax.set_yscale("log", base=base)


def add_text_box(ax, text: str, loc: str = "upper right",
                bbox_props: Optional[Dict] = None,
                fontsize: int = 10, **kwargs) -> None:
    """
    Add text box to axes.
    
    Args:
        ax: Matplotlib axes object
        text: Text to display
        loc: Box location (e.g., "upper right", "lower left")
        bbox_props: Dictionary of bbox properties
        fontsize: Font size
        **kwargs: Additional arguments passed to ax.text()
        
    Example:
        >>> add_text_box(ax, "Signal region", loc="upper left")
    """
    if bbox_props is None:
        bbox_props = {
            "boxstyle": "round,pad=0.5",
            "facecolor": "white",
            "alpha": 0.8,
            "edgecolor": "black",
        }
    
    # Map location names to coordinates
    loc_map = {
        "upper right": (0.96, 0.96),
        "upper left": (0.04, 0.96),
        "lower right": (0.96, 0.04),
        "lower left": (0.04, 0.04),
        "center": (0.5, 0.5),
    }
    
    x, y = loc_map.get(loc, (0.96, 0.96))
    
    ax.text(x, y, text, transform=ax.transAxes, fontsize=fontsize,
           verticalalignment="top" if y > 0.5 else "bottom",
           horizontalalignment="right" if x > 0.5 else "left",
           bbox=bbox_props, **kwargs)


def add_colorbar(im, ax: Optional[plt.Axes] = None,
                label: Optional[str] = None,
                orientation: str = "vertical") -> mpl.colorbar.Colorbar:
    """
    Add colorbar to figure with consistent styling.
    
    Args:
        im: Image object (returned from imshow, pcolormesh, etc.)
        ax: Axes object (required for constrained layout)
        label: Colorbar label
        orientation: "vertical" or "horizontal"
    
    Returns:
        Colorbar object
    """
    cbar = plt.colorbar(im, ax=ax, orientation=orientation)
    
    if label:
        if orientation == "vertical":
            cbar.set_label(label, rotation=270, labelpad=20)
        else:
            cbar.set_label(label, rotation=0, labelpad=15)
    
    return cbar


def add_stat_box(ax, hist_data: np.ndarray, 
                loc: str = "upper right",
                show_mean: bool = True, show_std: bool = True,
                show_entries: bool = True) -> None:
    """
    Add statistics box to histogram plot.
    
    Args:
        ax: Matplotlib axes object
        hist_data: Histogram data (bin contents)
        loc: Box location
        show_mean: Whether to show mean
        show_std: Whether to show standard deviation
        show_entries: Whether to show number of entries
    """
    # Filter out zero/negative entries for statistics
    data = hist_data[hist_data > 0]
    
    if len(data) == 0:
        return
    
    # Calculate statistics
    entries = len(data)
    mean = np.mean(data)
    std = np.std(data)
    
    # Build text
    text_lines = []
    if show_entries:
        text_lines.append(f"Entries: {entries}")
    if show_mean:
        text_lines.append(f"Mean: {mean:.2f}")
    if show_std:
        text_lines.append(f"Std Dev: {std:.2f}")
    
    text = "\n".join(text_lines)
    add_text_box(ax, text, loc=loc, fontsize=9)


def format_title_with_info(ax, title: str, info: Optional[Dict[str, str]] = None) -> None:
    """
    Set title with additional info in subtitle format.
    
    Args:
        ax: Matplotlib axes object
        title: Main title
        info: Dictionary of info to display below title
              (e.g., {"Region": "Signal", "Cuts": "Mx2"})
    
    Example:
        >>> format_title_with_info(ax, "Muon Kinematics",
        ...                        info={"Region": "Signal", "Cuts": "Mx2"})
    """
    ax.set_title(title, fontsize=13, fontweight="bold")
    
    if info:
        subtitle = ", ".join([f"{k}: {v}" for k, v in info.items()])
        ax.text(0.5, 0.98, subtitle, transform=ax.transAxes,
               fontsize=9, ha="center", va="top",
               style="italic", color="gray",
               transform_offset=(0, -15))


def create_ratio_plot(fig, gs, hist1: np.ndarray, hist2: np.ndarray,
                     bins: Optional[np.ndarray] = None,
                     xlabel: Optional[str] = None,
                     ylabel: Optional[str] = None,
                     title: Optional[str] = None,
                     label1: str = "Numerator",
                     label2: str = "Denominator") -> Tuple[plt.Axes, plt.Axes]:
    """
    Create a ratio plot with main histogram on top and ratio below.
    
    Args:
        fig: Matplotlib figure object
        gs: GridSpec object for layout
        hist1: Numerator histogram data
        hist2: Denominator histogram data
        bins: Bin edges (computed from data length if None)
        xlabel: X-axis label
        ylabel: Y-axis label for ratio
        title: Plot title
        label1: Label for first histogram
        label2: Label for second histogram
    
    Returns:
        Tuple of (main_ax, ratio_ax)
    """
    if bins is None:
        bins = np.arange(len(hist1) + 1)
    
    # Create subplots with shared x-axis
    ax_main = fig.add_subplot(gs[0])
    ax_ratio = fig.add_subplot(gs[1], sharex=ax_main)
    
    # Plot histograms
    bin_centers = (bins[:-1] + bins[1:]) / 2
    ax_main.hist(bin_centers, bins=bins, weights=hist1, 
                label=label1, alpha=0.7, edgecolor="black")
    ax_main.hist(bin_centers, bins=bins, weights=hist2,
                label=label2, alpha=0.7, edgecolor="black")
    
    # Plot ratio
    ratio = np.divide(hist1, hist2, where=hist2!=0, 
                     out=np.zeros_like(hist1, dtype=float))
    ax_ratio.plot(bin_centers, ratio, "ko-", markersize=5)
    ax_ratio.axhline(y=1, color="r", linestyle="--", linewidth=1.5, alpha=0.7)
    
    # Formatting
    ax_main.set_ylabel(ylabel or "Count")
    ax_main.legend()
    ax_main.set_title(title or "")
    
    if xlabel:
        ax_ratio.set_xlabel(xlabel)
    ax_ratio.set_ylabel("Ratio")
    ax_ratio.grid(True, alpha=0.3)
    ax_ratio.set_ylim([0.5, 1.5])
    
    # Hide x labels on main plot
    ax_main.tick_params(labelbottom=False)
    
    return ax_main, ax_ratio
