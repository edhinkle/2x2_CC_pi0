# plotdqmpy/plotting/plotters/beam_plotter.py
"""
Plotter for Beam histograms in CC1pi0 analysis.

Organizes plotting of all BeamHists histograms into logical groupings:
- Detector basics (spills, POT)
"""

from typing import Optional
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from .base_plotter import HistogramPlotter
from ..styles import get_color, get_colors, configure_axes
from ..utils import add_colorbar, get_axis_label, add_text_box


class BeamPlotter(HistogramPlotter):
    """
    Plotter for Beam histograms in DQM

    Provides organized methods to plot:
    - Detector basics (spills, POT)
    - POT by spill time
    """

    def __init__(self, hist_dict, sample_name: Optional[str] = None):
        
        self.hist_dict = hist_dict
        self.sample_name = sample_name

    def plot_total_data_counts(self, output_dir: Optional[str] = None) -> plt.Figure:
        """
        Plot counts of data collected: spills and POT processed.
        Histograms:
            - h_beam_TotalSpillsProcessed
            - h_beam_TotalPOT
            - h_beam_POTperSpill
        """

        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        
        total_spills = self.plot_1d('h_beam_TotalSpillsProcessed', normalized=False, ax=axes[0], color=get_colors(1))
        configure_axes(axes[0], ylabel=get_axis_label('count'), title=f"Total Spills Processed (Beam) for {self.sample_name}")
        add_text_box(axes[0], f"{total_spills:.0f} Spills", loc="upper right")
        
        total_pot = self.plot_1d('h_beam_TotalPOT', normalized=False, ax=axes[1], color=get_colors(1))
        configure_axes(axes[1], ylabel=get_axis_label('count'), title=f"Total POT Processed (Beam) for {self.sample_name}")
        total_pot_rescaled = total_pot * 1e13  # Saved as POT/1e13 for ease in histograms
        add_text_box(axes[1], f"{total_pot_rescaled:.2e} POT", loc="upper right")
        
        pot_per_spill = self.plot_1d('h_beam_POTperSpill', normalized=False, ax=axes[2], color=get_colors(1))
        configure_axes(axes[2], ylabel='POT/1e13 /Spill', title=f"POT by Spill (Beam) for {self.sample_name}")
        add_text_box(axes[2], f"{pot_per_spill:.2e} POT/Spill", loc="upper right")

        plt.tight_layout()
        if output_dir:
            fig.savefig(f'{output_dir}/pngs/{self.sample_name}_beam_total_data_counts.png', dpi=300, bbox_inches='tight')
            with PdfPages(f'{output_dir}/pdfs/{self.sample_name}_beam_total_data_counts.pdf') as pdf:
                pdf.savefig(fig, dpi=300, bbox_inches='tight')
        return fig
    
    def plot_POT_by_spill_time(self, output_dir: Optional[str] = None) -> plt.Figure:
        """
        Plot POT by spill time.
        TTrees:
            - t_beam_SpillPOTAndStartTime (branches : pot, start_time_sec, start_time_nsec)
        """
        fig, ax = plt.subplots(1, 1, figsize=(18, 5))

        spill_data = self.get_ttree('t_beam_SpillPOTAndStartTime')
        print(f"Retrieved spill data with {len(spill_data)} entries for plotting POT by spill time.")
        print(f"Spill data columns: {spill_data.columns.tolist()}")
        spill_pot = spill_data['pot'].to_numpy()*1e13  # Convert to POT
        print(f"Spill POT data (first 5 entries): {spill_pot[:5]}")
        spill_time = spill_data['start_time_sec'].to_numpy() + spill_data['start_time_nsec'].to_numpy()*1e-9  # Convert to seconds
        print(f"Spill time data (first 5 entries): {spill_time[:5]}")

        fig = self.plot_scatter(spill_time, spill_pot, ax=ax, xlabel='Spill Time [s]', ylabel='POT',\
                                 title=f"POT by Spill Time (Beam) for {self.sample_name}", color='k', alpha=0.4, s=1)
        ax.set_yscale('log')

        plt.tight_layout()
        if output_dir:
            fig.savefig(f'{output_dir}/pngs/{self.sample_name}_beam_pot_by_spill_time.png', dpi=300, bbox_inches='tight')
            with PdfPages(f'{output_dir}/pdfs/{self.sample_name}_beam_pot_by_spill_time.pdf') as pdf:
                pdf.savefig(fig, dpi=300, bbox_inches='tight')
        return fig

    
    def plot_all(self, output_dir: str = "/outputs/beam/") -> None:
        """Generate all BeamPlotter plots and save to output directory."""
        from pathlib import Path
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        Path(output_dir+"/pngs").mkdir(parents=True, exist_ok=True)
        Path(output_dir+"/pdfs").mkdir(parents=True, exist_ok=True)
        
        print(f"Generating Beam plots in {output_dir}...")
        fig1 = self.plot_total_data_counts(output_dir)
        plt.close(fig1)  # Close figure to save memory
        print("  ✓ Total data counts")
        fig2 = self.plot_POT_by_spill_time(output_dir)
        plt.close(fig2)  # Close figure to save memory
        print("  ✓ POT by spill time")
        print(f"Done! All plots saved to {output_dir}")