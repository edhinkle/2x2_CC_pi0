# python/plotting/plotters/reco_plotter.py
"""
Plotter for reconstructed (Reco) histograms in CC1pi0 analysis.

Organizes plotting of all RecoHists histograms into logical groupings:
- Detector basics (spills, POT)
- Vertex diagnostics (with/without cuts)
- Kinematics (muon angle)
- Particle multiplicities
"""

from typing import Optional
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from .base_plotter import HistogramPlotter
from ..styles import get_color, get_colors, configure_axes
from ..utils import add_colorbar, get_axis_label, add_text_box


class RecoPlotter(HistogramPlotter):
    """
    Plotter for reconstructed histograms from RecoHists.
    
    Provides organized methods to plot:
    - Detector basics (spills, POT)
    - Vertex with/without cuts
    - Vertex comparisons
    - Muon kinematics
    - Shower multiplicities
    """

    def __init__(self, hist_dict):
        self.hist_dict = hist_dict

    
    def plot_total_data_counts(self, output_dir: Optional[str] = None) -> plt.Figure:
        """
        Plot counts of data collected: spills and POT processed.
        Histograms:
            - h_reco_TotalSpillsProcessed
            - h_reco_TotalPOT
        """

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        total_spills = self.plot_1d('h_reco_TotalSpillsProcessed', normalized=False, ax=axes[0], color=get_colors(1))
        configure_axes(axes[0], ylabel=get_axis_label('count'), title="Total Spills Processed (Reco)")
        add_text_box(axes[0], f"{total_spills:.0f} Spills", loc="upper right")
        
        total_pot = self.plot_1d('h_reco_TotalPOT', normalized=False, ax=axes[1], color=get_colors(1))
        configure_axes(axes[1], ylabel=get_axis_label('count'), title="Total POT Processed (Reco)")
        total_pot_rescaled = total_pot * 1e13  # Saved as POT/1e13 for ease in histograms
        add_text_box(axes[1], f"{total_pot_rescaled:.2e} POT", loc="upper right")
        
        plt.tight_layout()
        if output_dir:
            fig.savefig(f'{output_dir}/pngs/reco_total_data_counts.png', dpi=300, bbox_inches='tight')
            with PdfPages(f'{output_dir}/pdfs/reco_total_data_counts.pdf') as pdf:
                pdf.savefig(fig, dpi=300, bbox_inches='tight')
        return fig
    
    def plot_vertex_no_cuts(self, output_dir: Optional[str] = None) -> plt.Figure:
        """
        Plot vertex with no cuts: XZ 2D + X,Y,Z 1D projections.
        Histograms:
            - h_reco_VertexZXNoCuts 
            - h_reco_VertexXNoCuts
            - h_reco_VertexYNoCut
            - h_reco_VertexZNoCuts
        """
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # 2D XZ
        mesh, vtx_xz_entries = self.plot_2d('h_reco_VertexZXNoCuts', normalized=True, ax=axes[0, 0], cmap='viridis')
        configure_axes(axes[0, 0], xlabel=get_axis_label("vertex_z"), \
                       ylabel=get_axis_label("vertex_x"), title="Reco Vertex XZ (No Cuts)")
        add_colorbar(mesh, axes[0, 0], label=get_axis_label('frac'))
        add_text_box(axes[0, 0], f"{vtx_xz_entries:.0f} Events", loc="upper right")
        
        # 1D projections
        vtx_x_entries = self.plot_1d('h_reco_VertexXNoCuts', normalized=True, ax=axes[0, 1], color=get_colors(1))
        configure_axes(axes[0, 1], xlabel=get_axis_label("vertex_x"), \
                       ylabel=get_axis_label('frac'), title="Reco Vertex X (No Cuts)")
        add_text_box(axes[0, 1], f"{vtx_x_entries:.0f} Events", loc="upper right")
        
        vtx_y_entries = self.plot_1d('h_reco_VertexYNoCuts', normalized=True, ax=axes[1, 0], color=get_colors(1))
        configure_axes(axes[1, 0], xlabel=get_axis_label("vertex_y"), \
                       ylabel=get_axis_label('frac'), title="Reco Vertex Y (No Cuts)")
        add_text_box(axes[1, 0], f"{vtx_y_entries:.0f} Events", loc="upper right")
        
        vtx_z_entries = self.plot_1d('h_reco_VertexZNoCuts', normalized=True, ax=axes[1, 1], color=get_colors(1))
        configure_axes(axes[1, 1], xlabel=get_axis_label("vertex_z"), \
                       ylabel=get_axis_label('frac'), title="Reco Vertex Z (No Cuts)")
        add_text_box(axes[1, 1], f"{vtx_z_entries:.0f} Events", loc="upper right")
        
        plt.tight_layout()
        if output_dir:
            fig.savefig(f'{output_dir}/pngs/reco_vertex_no_cuts.png', dpi=300, bbox_inches='tight')
            with PdfPages(f'{output_dir}/pdfs/reco_vertex_no_cuts.pdf') as pdf:
                pdf.savefig(fig, dpi=300, bbox_inches='tight')
        return fig
    
    def plot_vertex_with_cuts(self, output_dir: Optional[str] = None) -> plt.Figure:
        """
        Plot vertex with cuts: XZ 2D + X,Y,Z 1D projections.
        Histograms:
            - h_reco_VertexZXWithCuts 
            - h_reco_VertexXWithCuts
            - h_reco_VertexYWithCuts
            - h_reco_VertexZWithCuts
        """
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # 2D XZ
        mesh, vtx_xz_entries = self.plot_2d('h_reco_VertexZXWithCuts', normalized=True, ax=axes[0, 0], cmap='viridis')
        configure_axes(axes[0, 0], xlabel=get_axis_label("vertex_z"), \
                       ylabel=get_axis_label("vertex_x"), title="Reco Vertex XZ (With Cuts)")
        add_colorbar(mesh, axes[0, 0], label=get_axis_label('frac'))
        add_text_box(axes[0, 0], f"{vtx_xz_entries:.0f} Events", loc="upper right")
        
        # 1D projections
        vtx_x_entries = self.plot_1d('h_reco_VertexXWithCuts', normalized=True, ax=axes[0, 1], color=get_colors(1))
        configure_axes(axes[0, 1], xlabel=get_axis_label("vertex_x"), \
                       ylabel=get_axis_label('frac'), title="Reco Vertex X (With Cuts)")
        add_text_box(axes[0, 1], f"{vtx_x_entries:.0f} Events", loc="upper right")
        
        vtx_y_entries = self.plot_1d('h_reco_VertexYWithCuts', normalized=True, ax=axes[1, 0], color=get_colors(1))
        configure_axes(axes[1, 0], xlabel=get_axis_label("vertex_y"), \
                       ylabel=get_axis_label('frac'), title="Reco Vertex Y (With Cuts)")
        add_text_box(axes[1, 0], f"{vtx_y_entries:.0f} Events", loc="upper right")
        
        vtx_z_entries = self.plot_1d('h_reco_VertexZWithCuts', normalized=True, ax=axes[1, 1], color=get_colors(1))
        configure_axes(axes[1, 1], xlabel=get_axis_label("vertex_z"), \
                       ylabel=get_axis_label('frac'), title="Reco Vertex Z (With Cuts)")
        add_text_box(axes[1, 1], f"{vtx_z_entries:.0f} Events", loc="upper right")
        
        plt.tight_layout()
        if output_dir:
            fig.savefig(f'{output_dir}/pngs/reco_vertex_with_cuts.png', dpi=300, bbox_inches='tight')
            with PdfPages(f'{output_dir}/pdfs/reco_vertex_with_cuts.pdf') as pdf:
                pdf.savefig(fig, dpi=300, bbox_inches='tight')
        return fig
    
    def plot_muon_kinematics(self, output_dir: Optional[str] = None) -> plt.Figure:
        """
        Plot muon cos(angle with beam) with different zoom levels.
        Histograms:
            - h_reco_CosL
            - h_reco_CosL_zoomOut
        """
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        num_ev_truth_sig = self.plot_1d('h_reco_CosL', normalized=True, ax=axes[0], color=get_color("muon"))
        configure_axes(axes[0], xlabel=get_axis_label("cosL"), \
                       ylabel=get_axis_label('frac'), title="Muon Angle with Beam (Truth Sel Mx2 Match Cut)")
        add_text_box(axes[0], f"{num_ev_truth_sig:.0f} Events", loc="upper left")
        
        num_ev_total = self.plot_1d('h_reco_CosL_zoomOut', normalized=True, ax=axes[1], color=get_color("muon"))
        configure_axes(axes[1], xlabel=get_axis_label("cosL"), \
                       ylabel=get_axis_label('frac'), title="Muon Angle with Beam (All Reco Events)")
        add_text_box(axes[1], f"{num_ev_total:.0f} Events", loc="upper left")
        #add_text_box(axes[1], "All Reco Selected Events", loc="upper right")
        
        plt.tight_layout()
        if output_dir:
            fig.savefig(f'{output_dir}/pngs/reco_muon_kinematics.png', dpi=300, bbox_inches='tight')
            with PdfPages(f'{output_dir}/pdfs/reco_muon_kinematics.pdf') as pdf:
                pdf.savefig(fig, dpi=300, bbox_inches='tight')
        return fig
    
    def plot_shower_multiplicity(self, output_dir: Optional[str] = None) -> plt.Figure:
        """
        Plot shower multiplicities: primary/secondary electrons and photons.
        Histograms:
            - h_reco_PrimElectronMultiplicity
            - h_reco_PrimPhotonMultiplicity
            - h_reco_PrimShowerMultiplicity
            - h_reco_SecElectronMultiplicity
            - h_reco_SecPhotonMultiplicity
            - h_reco_SecShowerMultiplicity
        """
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        
        # Primary
        num_ev_prim_electron = self.plot_1d('h_reco_PrimElectronMultiplicity', normalized=True, ax=axes[0, 0], color=get_color("electron"))
        configure_axes(axes[0, 0], xlabel=get_axis_label('multiplicity'),\
                       ylabel=get_axis_label('frac'), title="Reco Primary Electrons")
        add_text_box(axes[0, 0], f"{num_ev_prim_electron:.0f} Events", loc="upper right")
        
        num_ev_prim_photon = self.plot_1d('h_reco_PrimPhotonMultiplicity', normalized=True, ax=axes[0, 1], color=get_color("photon"))
        configure_axes(axes[0, 1], xlabel=get_axis_label('multiplicity'),\
                       ylabel=get_axis_label('frac'), title="Reco Primary Photons")
        add_text_box(axes[0, 1], f"{num_ev_prim_photon:.0f} Events", loc="upper right")
        
        num_ev_prim_shower = self.plot_1d('h_reco_PrimShowerMultiplicity', normalized=True, ax=axes[0, 2], color=get_color("primary"))
        configure_axes(axes[0, 2], xlabel=get_axis_label('multiplicity'),\
                       ylabel=get_axis_label('frac'), title="Reco Primary Showers (e + γ)")
        add_text_box(axes[0, 2], f"{num_ev_prim_shower:.0f} Events", loc="upper right")
        
        # Secondary
        num_ev_sec_electron = self.plot_1d('h_reco_SecElectronMultiplicity', normalized=True, ax=axes[1, 0], color=get_color("electron"))
        configure_axes(axes[1, 0], xlabel=get_axis_label('multiplicity'),\
                       ylabel=get_axis_label('frac'),  title="Reco Secondary Electrons")
        add_text_box(axes[1, 0], f"{num_ev_sec_electron:.0f} Events", loc="upper right")
        
        num_ev_sec_photon = self.plot_1d('h_reco_SecPhotonMultiplicity', normalized=True, ax=axes[1, 1], color=get_color("photon"))
        configure_axes(axes[1, 1], xlabel=get_axis_label('multiplicity'),\
                       ylabel=get_axis_label('frac'), title="Reco Secondary Photons")
        add_text_box(axes[1, 1], f"{num_ev_sec_photon:.0f} Events", loc="upper right")
        
        num_ev_sec_shower = self.plot_1d('h_reco_SecShowerMultiplicity', normalized=True, ax=axes[1, 2], color=get_color("secondary"))
        configure_axes(axes[1, 2], xlabel=get_axis_label('multiplicity'),\
                       ylabel=get_axis_label('frac'), title="Reco Secondary Showers (e + γ)")
        add_text_box(axes[1, 2], f"{num_ev_sec_shower:.0f} Events", loc="upper right")

        plt.tight_layout()
        if output_dir:
            fig.savefig(f'{output_dir}/pngs/reco_shower_multiplicity.png', dpi=300, bbox_inches='tight')
            with PdfPages(f'{output_dir}/pdfs/reco_shower_multiplicity.pdf') as pdf:
                pdf.savefig(fig, dpi=300, bbox_inches='tight')
        return fig
    
    def plot_all(self, output_dir: str = "/outputs/reco/") -> None:
        """Generate all RecoPlotter plots and save to output directory."""
        from pathlib import Path
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        Path(output_dir+"/pngs").mkdir(parents=True, exist_ok=True)
        Path(output_dir+"/pdfs").mkdir(parents=True, exist_ok=True)
        
        print(f"Generating Reco plots in {output_dir}...")
        fig1 = self.plot_total_data_counts(output_dir)
        plt.close(fig1)  # Close figure to save memory
        print("  ✓ Total data counts")
        fig2 = self.plot_vertex_no_cuts(output_dir)
        plt.close(fig2)  # Close figure to save memory
        print("  ✓ Vertex (no cuts)")
        fig3 = self.plot_vertex_with_cuts(output_dir)
        plt.close(fig3)  # Close figure to save memory
        print("  ✓ Vertex (with cuts)")
        fig4 = self.plot_muon_kinematics(output_dir)
        plt.close(fig4)  # Close figure to save memory
        print("  ✓ Muon kinematics")
        fig5 = self.plot_shower_multiplicity(output_dir)
        plt.close(fig5)  # Close figure to save memory
        print("  ✓ Shower multiplicities")
        print(f"Done! All plots saved to {output_dir}")