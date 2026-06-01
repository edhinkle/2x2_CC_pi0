# python/plotting/plotters/reco_plotter.py
from .base_plotter import HistogramPlotter
import matplotlib.pyplot as plt

class RecoPlotter(HistogramPlotter):

    """
    Plotter for reco histograms for CC1pi0 analysis.
    """

    def plot_vertex_diagnostics(self, output_dir):
        """Plot vertex XZ, X, Y, Z with/without cuts"""
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        self.plot_1d('h_reco_VertexXNoCuts', ax=axes[0, 0], label='No Cuts')
        self.plot_1d('h_reco_VertexXWithCuts', ax=axes[0, 0], label='With Cuts')
        # ... etc
        
        fig.savefig(f'{output_dir}/reco_vertices.png')
    
    def plot_shower_multiplicity(self, output_dir):
        """Plot primary/secondary shower, electron, photon multiplicities"""
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        # ... organized by particle type
        fig.savefig(f'{output_dir}/reco_multiplicities.png')