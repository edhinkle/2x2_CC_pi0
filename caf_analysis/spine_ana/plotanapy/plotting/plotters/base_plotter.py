from abc import ABC, abstractmethod
from typing import Dict, List, Optional, Tuple, Any
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import uproot


class HistogramPlotter(ABC):
    """
    Abstract base class for histogram plotting for CC1pi0 analysis.
    
    Provides common functionality for 1D and 2D histogram visualization
    with consistent styling and layout management.
    
    Attributes:
        hist_dict (Dict): Dictionary mapping histogram names to histogram objects
        style_config (Dict): Configuration for plot styling
        sel_config (Dict): Selection configuration parameters
        beam_config (Dict): Beam configuration parameters
        det_config (Dict): Detector configuration parameters
    """
    
    def __init__(self, hist_dict: Optional[Dict[str, Any]] = None, 
                 style_config: Optional[Dict[str, Any]] = None,
                 sel_config: Optional[Dict[str, Any]] = None,
                 beam_config: Optional[Dict[str, Any]] = None,
                 det_config: Optional[Dict[str, Any]] = None):
        """
        Initialize the plotter.
        
        Args:
            hist_dict: Dictionary of histograms (name -> histogram object)
            style_config: Configuration dictionary for styling options
            sel_config: Selection configuration parameters
            beam_config: Beam configuration parameters
            det_config: Detector configuration parameters
        """
        self.hist_dict = hist_dict or {}
        self.style_config = style_config or {}
        self.sel_config = sel_config or {}
        self.beam_config = beam_config or {}
        self.det_config = det_config or {}

    
    def load_histograms(self, hist_dict: Dict[str, Any]) -> None:
        """
        Load histograms into the plotter.
        
        Args:
            hist_dict: Dictionary mapping histogram names to histogram objects
        """
        self.hist_dict.update(hist_dict)
    
    def add_histogram(self, name: str, hist: Any) -> None:
        """
        Add a single histogram to the plotter.
        
        Args:
            name: Histogram name/identifier
            hist: Histogram object
        """
        self.hist_dict[name] = hist
    
    def get_histogram(self, name: str) -> Any:
        """
        Retrieve a histogram by name.
        
        Args:
            name: Histogram name
            
        Returns:
            The histogram object, or None if not found
            
        Raises:
            KeyError: If histogram not found
        """
        if name not in self.hist_dict:
            raise KeyError(f"Histogram '{name}' not found. Available: {list(self.hist_dict.keys())}")
        return self.hist_dict[name]
    
    def plot_1d(self, hist_name: str, normalized: bool = True, 
                ax: Optional[plt.Axes] = None,
                label: Optional[str] = None, color: Optional[str] = None,
                linestyle: str = '-', linewidth: float = 1.5,
                alpha: float = 0.7, **kwargs) -> int:
        """
        Plot a 1D histogram.
        
        Args:
            hist_name: Name of the histogram to plot
            normalized: Whether to normalize the histogram to unit area
            ax: Matplotlib axes to plot on (creates new if None)
            label: Label for the histogram (uses hist_name if None)
            color: Line/fill color
            linestyle: Line style
            linewidth: Line width
            alpha: Transparency (0-1)
            **kwargs: Additional arguments passed to matplotlib
            
        Returns:
            Number of entries in the histogram
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 6))
        else:
            fig = ax.get_figure()
        
        ax = ax
        hist = self.get_histogram(hist_name)
        
        # Extract histogram data (handles both TH1D and numpy arrays)
        if hasattr(hist, 'to_numpy'):  # ROOT TH1D (loaded with uproot)
            data, bin_edges = hist.to_numpy()
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        else:  # numpy array or similar
            data = np.asarray(hist)
            bin_edges = np.arange(len(data) + 1)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        
        # Plot as histogram or line
        if 'drawstyle' not in kwargs:
            ax.hist(bin_centers, bins=bin_edges, 
                    weights=data/data.sum() if normalized else data, 
                    label=label or hist_name, color=color, 
                    alpha=alpha, linewidth=linewidth, **kwargs)
        else:
            normalized = False  # Line plots should not be normalized
            ax.scatter(bin_centers, data, label=label or hist_name, 
                    color=color, linestyle=linestyle, linewidth=linewidth, 
                    alpha=alpha, **kwargs)
        
        return data.sum()
    
    def plot_cutflow(self, hist_name: str, truth_events_total: Optional[int] = None,
                     axes: Optional[plt.Axes] = None,
                     label: Optional[str] = None, color: Optional[str] = None,
                     linestyle: str = '-', linewidth: float = 1.5,
                     alpha: float = 0.7, **kwargs) -> Tuple[plt.figure, plt.axes]:
        """
        Plot a cutflow histogram.
        
        Args:
            hist_name: Name of the histogram to plot
            truth_events_total: Number of true signal events for efficiency/purity
            axes: Matplotlib axes to plot on (creates new if None)
            label: Label for the histogram (uses hist_name if None)
            color: Line/fill color
            linestyle: Line style
            linewidth: Line width
            alpha: Transparency (0-1)
            **kwargs: Additional arguments passed to matplotlib
            
        Returns:
            Tuple of figure, axes
        """
        if axes is None:
            fig, axes = plt.subplots(2,1,figsize=(10, 8), sharex=True, 
                                           gridspec_kw={"height_ratios": [2, 1]})
        else:
            fig = axes[0].get_figure()
        
        hist = self.get_histogram(hist_name)
        
        # Extract histogram data (handles both TH1D and numpy arrays)
        if hasattr(hist, 'to_numpy'):  # ROOT TH1D (loaded with uproot)
            values = hist.values()
            cut_labels = list(hist.axis().labels())
            bin_edges = np.arange(len(values) + 1)
            bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
            print(values, bin_centers)
        else:  # numpy array or similar
            data = np.asarray(hist)
            bin_edges = np.arange(len(data) + 1)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        
        # Plot as histogram
        axes[0].stairs(values, bin_edges, linewidth=2,label=label)
        axes[0].set_yscale("log")
        axes[0].set_xticks(np.arange(len(values)) + 0.5)
        axes[0].set_xticklabels(cut_labels, rotation=45, ha="right")

        # Event counts above each cut
        for i, y in enumerate(values):

            axes[0].text(i+0.5, y * 1.15, f"{int(y)}", ha="center", fontsize=12, rotation=0)

        return fig, axes, cut_labels, values


    def plot_2d(self, hist_name: str = None, normalized: bool = True, 
                ax: Optional[plt.Axes] = None,
                cmap: str = 'cividis', **kwargs) -> Tuple[mpl.collections.QuadMesh, int]:
        """
        Plot a 2D histogram as a heatmap.
        
        Args:
            hist_name: Name of the 2D histogram to plot
            normalized: Whether to normalize the histogram to unit area
            ax: Matplotlib axes to plot on (creates new if None)
            cmap: Colormap name
            **kwargs: Additional arguments passed to matplotlib imshow
            
        Returns:
            Tuple of (mesh, number of hist entries)
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 8))
        else:
            fig = ax.get_figure()
        
        hist = self.get_histogram(hist_name)
        
        # Extract 2D histogram data
        if hasattr(hist, 'to_numpy'):  # ROOT TH1D (loaded with uproot)
            data, x_edges, y_edges = hist.to_numpy()
        else:
            data = np.asarray(hist)
        
        # Plot heatmap
        mesh = ax.pcolormesh(
            x_edges,
            y_edges,
            (data / np.sum(data) if normalized else data).T, # data = (x,y) but pcolormesh wants (y,x), so transpose
            cmap=cmap,
            **kwargs,
        )
        #ax.set_title(hist_name)
        #plt.colorbar(mesh, ax=ax)
        
        return mesh, np.sum(data)
    
    
    @abstractmethod
    def plot_all(self, output_dir: str) -> None:
        """
        Generate all plots for this histogram set.
        
        Subclasses should implement this to create all relevant plots
        and save them to output_dir.
        
        Args:
            output_dir: Directory to save plots to
        """
        pass
