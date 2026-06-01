from abc import ABC, abstractmethod
from typing import Dict, List, Optional, Tuple, Any
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import uproot


class HistogramPlotter(ABC):
    """
    Abstract base class for histogram plotting for CC1pi0 analysis.
    
    Provides common functionality for 1D and 2D histogram visualization
    with consistent styling and layout management.
    
    Attributes:
        hist_dict (Dict): Dictionary mapping histogram names to histogram objects
        style_config (Dict): Configuration for plot styling
        fig (plt.Figure): Current figure object
        ax (plt.Axes or np.ndarray): Current axes object(s)
    """
    
    def __init__(self, hist_dict: Optional[Dict[str, Any]] = None, 
                 style_config: Optional[Dict[str, Any]] = None):
        """
        Initialize the plotter.
        
        Args:
            hist_dict: Dictionary of histograms (name -> histogram object)
            style_config: Configuration dictionary for styling options
        """
        self.hist_dict = hist_dict or {}
        self.style_config = style_config or {}
        self.fig = None
        self.ax = None
    
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
    
    def plot_1d(self, hist_name: str, ax: Optional[plt.Axes] = None,
                label: Optional[str] = None, color: Optional[str] = None,
                linestyle: str = '-', linewidth: float = 1.5,
                alpha: float = 0.7, **kwargs) -> Tuple[plt.Figure, plt.Axes]:
        """
        Plot a 1D histogram.
        
        Args:
            hist_name: Name of the histogram to plot
            ax: Matplotlib axes to plot on (creates new if None)
            label: Label for the histogram (uses hist_name if None)
            color: Line/fill color
            linestyle: Line style
            linewidth: Line width
            alpha: Transparency (0-1)
            **kwargs: Additional arguments passed to matplotlib
            
        Returns:
            Tuple of (figure, axes)
        """
        if ax is None:
            self.fig, ax = plt.subplots(figsize=(10, 6))
        else:
            self.fig = ax.get_figure()
        
        self.ax = ax
        hist = self.get_histogram(hist_name)
        
        # Extract histogram data (handles both TH1D and numpy arrays)
        if hasattr(hist, 'GetArray'):  # ROOT TH1D
            data = hist.GetArray()
            bin_edges = self._get_bin_edges_root(hist)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        else:  # numpy array or similar
            data = np.asarray(hist)
            bin_edges = np.arange(len(data) + 1)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        
        # Plot as histogram or line
        if 'drawstyle' not in kwargs:
            ax.hist(bin_centers, bins=bin_edges, weights=data[:-1], 
                   label=label or hist_name, color=color, 
                   alpha=alpha, linewidth=linewidth, **kwargs)
        else:
            ax.plot(bin_centers, data[:-1], label=label or hist_name, 
                   color=color, linestyle=linestyle, linewidth=linewidth, 
                   alpha=alpha, **kwargs)
        
        return self.fig, ax
    
    def plot_2d(self, hist_name: str, ax: Optional[plt.Axes] = None,
                cmap: str = 'viridis', **kwargs) -> Tuple[plt.Figure, plt.Axes]:
        """
        Plot a 2D histogram as a heatmap.
        
        Args:
            hist_name: Name of the 2D histogram to plot
            ax: Matplotlib axes to plot on (creates new if None)
            cmap: Colormap name
            **kwargs: Additional arguments passed to matplotlib imshow
            
        Returns:
            Tuple of (figure, axes)
        """
        if ax is None:
            self.fig, ax = plt.subplots(figsize=(10, 8))
        else:
            self.fig = ax.get_figure()
        
        self.ax = ax
        hist = self.get_histogram(hist_name)
        
        # Extract 2D histogram data
        if hasattr(hist, 'GetArray'):  # ROOT TH2D
            data = self._get_2d_array_root(hist)
        else:
            data = np.asarray(hist)
        
        # Plot heatmap
        im = ax.imshow(data, cmap=cmap, aspect='auto', origin='lower', **kwargs)
        ax.set_title(hist_name)
        plt.colorbar(im, ax=ax)
        
        return self.fig, ax
    
    def create_figure(self, n_rows: int = 1, n_cols: int = 1,
                     figsize: Tuple[float, float] = (10, 6),
                     **kwargs) -> Tuple[plt.Figure, Any]:
        """
        Create a new figure with subplots.
        
        Args:
            n_rows: Number of subplot rows
            n_cols: Number of subplot columns
            figsize: Figure size (width, height) in inches
            **kwargs: Additional arguments passed to plt.subplots()
            
        Returns:
            Tuple of (figure, axes)
        """
        self.fig, self.ax = plt.subplots(n_rows, n_cols, figsize=figsize, **kwargs)
        return self.fig, self.ax
    
    def save(self, output_path: str, dpi: int = 300, 
             bbox_inches: str = 'tight', **kwargs) -> None:
        """
        Save the current figure.
        
        Args:
            output_path: Path to save the figure (e.g., 'plot.png')
            dpi: Resolution in dots per inch
            bbox_inches: Bounding box setting
            **kwargs: Additional arguments passed to fig.savefig()
        """
        if self.fig is None:
            raise RuntimeError("No figure to save. Create a plot first.")
        
        self.fig.savefig(output_path, dpi=dpi, bbox_inches=bbox_inches, **kwargs)
        print(f"Figure saved to {output_path}")
    
    def show(self) -> None:
        """Display the current figure."""
        if self.fig is None:
            raise RuntimeError("No figure to display.")
        plt.show()
    
    def close(self) -> None:
        """Close the current figure."""
        if self.fig is not None:
            plt.close(self.fig)
            self.fig = None
            self.ax = None
    
    @staticmethod
    def _get_bin_edges_root(hist) -> np.ndarray:
        """
        Extract bin edges from ROOT TH1D histogram.
        
        Args:
            hist: ROOT TH1D histogram
            
        Returns:
            Array of bin edges
        """
        n_bins = hist.GetNbinsX()
        edges = []
        for i in range(n_bins + 1):
            edges.append(hist.GetBinLowEdge(i + 1))
        return np.array(edges)
    
    @staticmethod
    def _get_2d_array_root(hist) -> np.ndarray:
        """
        Extract 2D array from ROOT TH2D histogram.
        
        Args:
            hist: ROOT TH2D histogram
            
        Returns:
            2D numpy array
        """
        nx = hist.GetNbinsX()
        ny = hist.GetNbinsY()
        data = np.zeros((ny, nx))
        
        for i in range(1, nx + 1):
            for j in range(1, ny + 1):
                data[j - 1, i - 1] = hist.GetBinContent(i, j)
        
        return data
    
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
