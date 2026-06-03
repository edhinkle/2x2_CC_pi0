"""
Styling and theme utilities for plots.
"""

from .theme import (
    apply_theme,
    get_colors,
    get_color,
    get_font_config,
    get_tick_config,
    get_legend_config,
    get_grid_config,
    configure_axes,
    set_color_cycle,
    COLORS,
    CATEGORICAL_COLORS,
    FONT_CONFIG,
    LEGEND_CONFIG,
    GRID_CONFIG,
)

__all__ = [
    # Functions
    "apply_theme",
    "get_colors",
    "get_color",
    "get_font_config",
    "get_tick_config",
    "get_legend_config",
    "get_grid_config",
    "configure_axes",
    "set_color_cycle",
    # Constants (color dicts, config dicts)
    "COLORS",
    "CATEGORICAL_COLORS",
    "FONT_CONFIG",
    "LEGEND_CONFIG",
    "GRID_CONFIG",
]