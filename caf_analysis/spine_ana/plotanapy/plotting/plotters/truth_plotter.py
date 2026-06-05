"""
Plotter for truth histograms in CC1pi0 analysis.

Organizes plotting of all TruthHists histograms into logical groupings:
- Interactions above KE threshold per spill
- Muon kinematics (cos(angle with beam)) and energy)
- Neutrino energy
- Interaction mode
- Primary pi0 multiplicity
- Secondary pi0 multiplicity pre- and post- "Mx2" truth cuts
- Primary and secondary e, gamma, and shower multiplicities pre- "Mx2" truth cuts
- Primary and secondary e, gamma, and shower multiplicities post- "Mx2" truth cuts
"""

from typing import Optional
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from .base_plotter import HistogramPlotter
from ..styles import get_color, get_colors, configure_axes
from ..utils import add_colorbar, get_axis_label, add_text_box, set_axis_limits


class TruthPlotter(HistogramPlotter):
    """
    Plotter for truth histograms in CC1pi0 analysis.

    Organizes plotting of all TruthHists histograms into logical groupings:
    - Interactions above KE threshold per spill
    - Muon kinematics (cos(angle with beam)) and energy)
    - Neutrino energy
    - Interaction mode
    - Primary pi0 multiplicity
    - Secondary pi0 multiplicity pre- and post- "Mx2" truth cuts
    - Primary and secondary e, gamma, and shower multiplicities pre- "Mx2" truth cuts
    - Primary and secondary e, gamma, and shower multiplicities post- "Mx2" truth cuts
    """

    def __init__(self, hist_dict, sel_config, beam_config, det_config):
        
        self.hist_dict = hist_dict
        self.sel_config = sel_config
        self.beam_config = beam_config
        self.det_config = det_config

    def plot_ixn_above_ke_threshold_per_spill(self, output_dir: Optional[str] = None) -> plt.Figure:
        """
        Plot counts of data collected: spills and POT processed.
        Histograms:
            - h_true_ixnsAboveKEThresholdPerSpill
        """

        fig, ax = plt.subplots(1, 1, figsize=(6, 5))
        
        thresh = self.det_config["keThreshold"]*1000 # Given in GeV in config, want MeV for plot title

        total_spills = self.plot_1d('h_true_ixnsAboveKEThresholdPerSpill', normalized=True, ax=ax, color=get_colors(1))
        configure_axes(ax, ylabel=get_axis_label('frac_sp'), title=f"True Interaction above KE Threshold of {thresh:.0f} MeV per Spill")
        add_text_box(ax, f"{total_spills:.0f} Spills", loc="upper right")
        

        plt.tight_layout()
        if output_dir:
            fig.savefig(f'{output_dir}/pngs/truth_ixn_above_ke_threshold_{thresh}MeV.png', dpi=300, bbox_inches='tight')
            with PdfPages(f'{output_dir}/pdfs/truth_ixn_above_ke_threshold_{thresh}MeV.pdf') as pdf:
                pdf.savefig(fig, dpi=300, bbox_inches='tight')
        return fig
    
    def plot_pre_mx2_muon_angle_numu_numubar(self, output_dir: Optional[str] = None) -> plt.Figure:
        """
        Plot muon cos(angle with beam) with different pre- Mx2 Truth Match "cheat" cut for numu/numubar events.
        Histograms:
            - h_true_CosLNumubar_zoomOut
            - h_true_CosLNumu_zoomOut
        """
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        muonECut = self.sel_config["muonEnergyCut"]
        cosThetaCut = self.sel_config["cosThetaCut"]
        
        num_ev_numu = self.plot_1d('h_true_CosLNumu_zoomOut', normalized=True, ax=axes[0], color=get_color("muon"))
        configure_axes(axes[0], xlabel=get_axis_label("cosL"), \
                       ylabel=get_axis_label('frac'), title=r"$\nu_\mu$ Events")
        add_text_box(axes[0], f"{num_ev_numu:.0f} Events", loc="upper left")
        
        num_ev_numubar = self.plot_1d('h_true_CosLNumubar_zoomOut', normalized=True, ax=axes[1], color=get_color("muon"))
        configure_axes(axes[1], xlabel=get_axis_label("cosL"), \
                       ylabel=get_axis_label('frac'), title=r"$\bar{\nu}_\mu$ Events")
        add_text_box(axes[1], f"{num_ev_numubar:.0f} Events", loc="upper left")

        fig.suptitle(fr"True $\mu$ Angle w/ Beam (Pre-Mx2 Cuts: >{muonECut:.1f} GeV, $\cos(\theta)$>{cosThetaCut:.2f}")
        
        plt.tight_layout()
        if output_dir:
            fig.savefig(f'{output_dir}/pngs/truth_muon_angle_pre_mx2_cheat_cut_numu_numubar.png', dpi=300, bbox_inches='tight')
            with PdfPages(f'{output_dir}/pdfs/truth_muon_angle_pre_mx2_cheat_cut_numu_numubar.pdf') as pdf:
                pdf.savefig(fig, dpi=300, bbox_inches='tight')
        return fig

    def plot_pre_post_mx2_muon_angle(self, output_dir: Optional[str] = None) -> plt.Figure:
        """
        Plot muon cos(angle with beam) with different pre- and post Mx2 Truth Match "cheat" cut. 
        Histograms:
            - h_true_CosL
            - h_true_CosL_zoomOut
        """
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        muonECut = self.sel_config["muonEnergyCut"]
        cosThetaCut = self.sel_config["cosThetaCut"]
        
        num_ev_truth_sig = self.plot_1d('h_true_CosL_zoomOut', normalized=True, ax=axes[0], color=get_color("muon"))
        configure_axes(axes[0], xlabel=get_axis_label("cosL"), \
                       ylabel=get_axis_label('frac'), title=fr"True $\mu$ Angle w/ Beam (Pre-Mx2 Cuts: >{muonECut:.1f} GeV, $\cos(\theta)$>{cosThetaCut:.2f})")
        add_text_box(axes[0], f"{num_ev_truth_sig:.0f} Events", loc="upper left")
        
        num_ev_total = self.plot_1d('h_true_CosL', normalized=True, ax=axes[1], color=get_color("muon"))
        configure_axes(axes[1], xlabel=get_axis_label("cosL"), \
                       ylabel=get_axis_label('frac'), title=fr"True $\mu$ Angle w/ Beam (Post-Mx2 Cuts: >{muonECut:.1f} GeV, $\cos(\theta)$>{cosThetaCut:.2f})")
        add_text_box(axes[1], f"{num_ev_total:.0f} Events", loc="upper left")
        
        plt.tight_layout()
        if output_dir:
            fig.savefig(f'{output_dir}/pngs/truth_muon_angle_pre_mx2_cheat_cut.png', dpi=300, bbox_inches='tight')
            with PdfPages(f'{output_dir}/pdfs/truth_muon_angle_pre_mx2_cheat_cut.pdf') as pdf:
                pdf.savefig(fig, dpi=300, bbox_inches='tight')
        return fig
    
    def plot_post_mx2_muon_energy(self, output_dir: Optional[str] = None) -> plt.Figure:
        """
        Plot muon energy.
        Plots filled using Mx2 Truth Match "cheat" cut. 
        Histograms:
            - h_true_Elep
        """
        fig, ax = plt.subplots(1, 1, figsize=(6, 5))

        muonECut = self.sel_config["muonEnergyCut"]
        cosThetaCut = self.sel_config["cosThetaCut"]
        
        num_ev_truth = self.plot_1d('h_true_Elep', normalized=True, ax=ax, color=get_color("muon"))
        configure_axes(ax, xlabel=get_axis_label("E", unit="GeV"), \
                       ylabel=get_axis_label('frac'), title=fr"True Muon Energy (Post-Mx2 Cuts: >{muonECut:.1f} GeV, $\cos(\theta)$>{cosThetaCut:.2f})")
        add_text_box(ax, f"{num_ev_truth:.0f} Events", loc="upper right")
        
        plt.tight_layout()
        if output_dir:
            fig.savefig(f'{output_dir}/pngs/truth_muon_energy_post_mx2_cheat_cut.png', dpi=300, bbox_inches='tight')
            with PdfPages(f'{output_dir}/pdfs/truth_muon_energy_post_mx2_cheat_cut.pdf') as pdf:
                pdf.savefig(fig, dpi=300, bbox_inches='tight')
        return fig
    
    def plot_selected_ixn_nu_energy(self, output_dir: Optional[str] = None) -> plt.Figure:
        """
        Plot neutrino energy.
        Plots filled AFTER Mx2 Truth Match "cheat" cut. 
        Histograms:
            - h_true_Enu
        """
        fig, ax = plt.subplots(1, 1, figsize=(6, 5))
        
        num_ev_truth = self.plot_1d('h_true_Enu', normalized=True, ax=ax, color=get_color("signal"))
        configure_axes(ax, xlabel=get_axis_label("E", unit="GeV"), \
                       ylabel=get_axis_label('frac'), title=fr"True Neutrino Energy")
        add_text_box(ax, f"{num_ev_truth:.0f} Events", loc="upper right")
        
        plt.tight_layout()
        if output_dir:
            fig.savefig(f'{output_dir}/pngs/truth_nu_energy.png', dpi=300, bbox_inches='tight')
            with PdfPages(f'{output_dir}/pdfs/truth_nu_energy.pdf') as pdf:
                pdf.savefig(fig, dpi=300, bbox_inches='tight')
        return fig
    
    def plot_pre_mx2_pre1pi0_primary_pi0_mult(self, output_dir: Optional[str] = None) -> plt.Figure:
        """
        Plot primary pi0 multiplicity before cut on one primary pi0.
        Plots filled prior to Mx2 Truth Match "cheat" cut as well.
        Histograms:
            - h_true_nPrimPi0
        """
        fig, ax = plt.subplots(1, 1, figsize=(6, 5))

        muonECut = self.sel_config["muonEnergyCut"]
        cosThetaCut = self.sel_config["cosThetaCut"]
        
        num_ev_truth = self.plot_1d('h_true_nPrimPi0', normalized=True, ax=ax, color=get_color("pion"))
        configure_axes(ax, xlabel=get_axis_label("multiplicity"), \
                       ylabel=get_axis_label('frac'), title=fr"True Primary $\pi^0$ Multiplicity (Pre 1 $\pi^0$ and Mx2 Cuts: >{muonECut:.1f} GeV, $\cos(\theta)$>{cosThetaCut:.2f})")
        add_text_box(ax, f"{num_ev_truth:.0f} Events", loc="upper right")
        
        plt.tight_layout()
        if output_dir:
            fig.savefig(f'{output_dir}/pngs/truth_primary_pi0_multiplicity_pre_pi0_pre_mx2.png', dpi=300, bbox_inches='tight')
            with PdfPages(f'{output_dir}/pdfs/truth_primary_pi0_multiplicity_pre_pi0_pre_mx2.pdf') as pdf:
                pdf.savefig(fig, dpi=300, bbox_inches='tight')
        return fig
    
    def plot_pre_post_mx2_sec_pi0_mult(self, output_dir: Optional[str] = None) -> plt.Figure:
        """
        Plot secondary pi0 multiplicity pre- and post-Mx2 Truth Match "cheat" cut. 
        Histograms:
            - h_true_nSecPi0_preMx2
            - h_true_nSecPi0_postMx2
        """
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        muonECut = self.sel_config["muonEnergyCut"]
        cosThetaCut = self.sel_config["cosThetaCut"]
        
        num_ev_truth_pre_mx2 = self.plot_1d('h_true_nSecPi0_preMx2', normalized=True, ax=axes[0], color=get_color("muon"))
        configure_axes(axes[0], xlabel=get_axis_label("multiplicity"), \
                       ylabel=get_axis_label('frac'), title=fr"True Sec $\pi^0$ Multiplicity (Pre-Mx2 Cuts: >{muonECut:.1f} GeV, $\cos(\theta)$>{cosThetaCut:.2f})")
        add_text_box(axes[0], f"{num_ev_truth_pre_mx2:.0f} Events", loc="upper right")
        
        num_ev_truth_post_mx2 = self.plot_1d('h_true_nSecPi0_postMx2', normalized=True, ax=axes[1], color=get_color("muon"))
        configure_axes(axes[1], xlabel=get_axis_label("multiplicity"), \
                       ylabel=get_axis_label('frac'), title=fr"True Sec $\pi^0$ Multiplicity (Post-Mx2 Cuts: >{muonECut:.1f} GeV, $\cos(\theta)$>{cosThetaCut:.2f})")
        add_text_box(axes[1], f"{num_ev_truth_post_mx2:.0f} Events", loc="upper right")
        
        plt.tight_layout()
        if output_dir:
            fig.savefig(f'{output_dir}/pngs/truth_sec_pi0_mult_pre_and_post_mx2_cheat_cut.png', dpi=300, bbox_inches='tight')
            with PdfPages(f'{output_dir}/pdfs/truth_sec_pi0_mult_pre_and_post_mx2_cheat_cut.pdf') as pdf:
                pdf.savefig(fig, dpi=300, bbox_inches='tight')
        return fig
    
    def plot_shower_multiplicity_pre_mx2(self, output_dir: Optional[str] = None) -> plt.Figure:
        """
        Plot shower multiplicities: primary/secondary electrons and photons pre-Mx2 Truth Cheat Cut
        Histograms:
            - h_true_nPrimElectron_preMx2
            - h_true_nSecElectron_preMx2
            - h_true_nPrimPhoton_preMx2
            - h_true_nSecPhoton_preMx2
            - h_true_nPrimShower_preMx2
            - h_true_nSecShower_preMx2
        """
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))

        muonECut = self.sel_config["muonEnergyCut"]
        cosThetaCut = self.sel_config["cosThetaCut"]
        
        # Primary
        num_ev_prim_electron = self.plot_1d('h_true_nPrimElectron_preMx2', normalized=True, ax=axes[0, 0], color=get_color("electron"))
        configure_axes(axes[0, 0], xlabel=get_axis_label('multiplicity'),\
                       ylabel=get_axis_label('frac'), title="True Primary Electrons")
        add_text_box(axes[0, 0], f"{num_ev_prim_electron:.0f} Events", loc="upper right")
        
        num_ev_prim_photon = self.plot_1d('h_true_nPrimPhoton_preMx2', normalized=True, ax=axes[0, 1], color=get_color("photon"))
        configure_axes(axes[0, 1], xlabel=get_axis_label('multiplicity'),\
                       ylabel=get_axis_label('frac'), title="True Primary Photons")
        add_text_box(axes[0, 1], f"{num_ev_prim_photon:.0f} Events", loc="upper right")
        
        num_ev_prim_shower = self.plot_1d('h_true_nPrimShower_preMx2', normalized=True, ax=axes[0, 2], color=get_color("primary"))
        configure_axes(axes[0, 2], xlabel=get_axis_label('multiplicity'),\
                       ylabel=get_axis_label('frac'), title="True Primary Showers (e + γ)")
        add_text_box(axes[0, 2], f"{num_ev_prim_shower:.0f} Events", loc="upper right")
        
        # Secondary
        num_ev_sec_electron = self.plot_1d('h_true_nSecElectron_preMx2', normalized=True, ax=axes[1, 0], color=get_color("electron"))
        configure_axes(axes[1, 0], xlabel=get_axis_label('multiplicity'),\
                       ylabel=get_axis_label('frac'),  title="True Secondary Electrons")
        add_text_box(axes[1, 0], f"{num_ev_sec_electron:.0f} Events", loc="upper right")
        
        num_ev_sec_photon = self.plot_1d('h_true_nSecPhoton_preMx2', normalized=True, ax=axes[1, 1], color=get_color("photon"))
        configure_axes(axes[1, 1], xlabel=get_axis_label('multiplicity'),\
                       ylabel=get_axis_label('frac'), title="True Secondary Photons")
        add_text_box(axes[1, 1], f"{num_ev_sec_photon:.0f} Events", loc="upper right")
        
        num_ev_sec_shower = self.plot_1d('h_true_nSecShower_preMx2', normalized=True, ax=axes[1, 2], color=get_color("secondary"))
        configure_axes(axes[1, 2], xlabel=get_axis_label('multiplicity'),\
                       ylabel=get_axis_label('frac'), title="True Secondary Showers (e + γ)")
        add_text_box(axes[1, 2], f"{num_ev_sec_shower:.0f} Events", loc="upper right")

        fig.suptitle(fr"Pre-Mx2 Cuts: >{muonECut:.1f} GeV, $\cos(\theta)$>{cosThetaCut:.2f}")

        plt.tight_layout()
        if output_dir:
            fig.savefig(f'{output_dir}/pngs/true_shower_multiplicity_pre_mx2.png', dpi=300, bbox_inches='tight')
            with PdfPages(f'{output_dir}/pdfs/true_shower_multiplicity_pre_mx2.pdf') as pdf:
                pdf.savefig(fig, dpi=300, bbox_inches='tight')
        return fig
    
    def plot_shower_multiplicity_post_mx2(self, output_dir: Optional[str] = None) -> plt.Figure:
        """
        Plot shower multiplicities: primary/secondary electrons and photons post-Mx2 Truth Cheat Cut
        Histograms:
            - h_true_nPrimElectron_postMx2
            - h_true_nSecElectron_postMx2
            - h_true_nPrimPhoton_postMx2
            - h_true_nSecPhoton_postMx2
            - h_true_nPrimShower_postMx2
            - h_true_nSecShower_postMx2
        """
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))

        muonECut = self.sel_config["muonEnergyCut"]
        cosThetaCut = self.sel_config["cosThetaCut"]
        
        # Primary
        num_ev_prim_electron = self.plot_1d('h_true_nPrimElectron_postMx2', normalized=True, ax=axes[0, 0], color=get_color("electron"))
        configure_axes(axes[0, 0], xlabel=get_axis_label('multiplicity'),\
                       ylabel=get_axis_label('frac'), title="True Primary Electrons")
        add_text_box(axes[0, 0], f"{num_ev_prim_electron:.0f} Events", loc="upper right")
        
        num_ev_prim_photon = self.plot_1d('h_true_nPrimPhoton_postMx2', normalized=True, ax=axes[0, 1], color=get_color("photon"))
        configure_axes(axes[0, 1], xlabel=get_axis_label('multiplicity'),\
                       ylabel=get_axis_label('frac'), title="True Primary Photons")
        add_text_box(axes[0, 1], f"{num_ev_prim_photon:.0f} Events", loc="upper right")
        
        num_ev_prim_shower = self.plot_1d('h_true_nPrimShower_postMx2', normalized=True, ax=axes[0, 2], color=get_color("primary"))
        configure_axes(axes[0, 2], xlabel=get_axis_label('multiplicity'),\
                       ylabel=get_axis_label('frac'), title="True Primary Showers (e + γ)")
        add_text_box(axes[0, 2], f"{num_ev_prim_shower:.0f} Events", loc="upper right")
        
        # Secondary
        num_ev_sec_electron = self.plot_1d('h_true_nSecElectron_postMx2', normalized=True, ax=axes[1, 0], color=get_color("electron"))
        configure_axes(axes[1, 0], xlabel=get_axis_label('multiplicity'),\
                       ylabel=get_axis_label('frac'),  title="True Secondary Electrons")
        add_text_box(axes[1, 0], f"{num_ev_sec_electron:.0f} Events", loc="upper right")
        
        num_ev_sec_photon = self.plot_1d('h_true_nSecPhoton_postMx2', normalized=True, ax=axes[1, 1], color=get_color("photon"))
        configure_axes(axes[1, 1], xlabel=get_axis_label('multiplicity'),\
                       ylabel=get_axis_label('frac'), title="True Secondary Photons")
        add_text_box(axes[1, 1], f"{num_ev_sec_photon:.0f} Events", loc="upper right")
        
        num_ev_sec_shower = self.plot_1d('h_true_nSecShower_postMx2', normalized=True, ax=axes[1, 2], color=get_color("secondary"))
        configure_axes(axes[1, 2], xlabel=get_axis_label('multiplicity'),\
                       ylabel=get_axis_label('frac'), title="True Secondary Showers (e + γ)")
        add_text_box(axes[1, 2], f"{num_ev_sec_shower:.0f} Events", loc="upper right")

        fig.suptitle(fr"Post-Mx2 Cuts: >{muonECut:.1f} GeV, $\cos(\theta)$>{cosThetaCut:.2f}")

        plt.tight_layout()
        if output_dir:
            fig.savefig(f'{output_dir}/pngs/true_shower_multiplicity_post_mx2.png', dpi=300, bbox_inches='tight')
            with PdfPages(f'{output_dir}/pdfs/true_shower_multiplicity_post_mx2.pdf') as pdf:
                pdf.savefig(fig, dpi=300, bbox_inches='tight')
        return fig
    
    def plot_true_signal_ixn_mode(self, output_dir: Optional[str] = None) -> plt.Figure:
        """
        Plot GENIE interaction mode for truth signal events.
        Plots filled AFTER Mx2 Truth Match "cheat" cut. 
        Histograms:
            - h_true_IxnMode
        """
        fig, ax = plt.subplots(1, 1, figsize=(6, 5))
        
        num_ev_truth = self.plot_1d('h_true_IxnMode', normalized=True, ax=ax, color=get_color("signal"))
        configure_axes(ax, xlabel="GENIE Interaction Modes", \
                       ylabel=get_axis_label('frac'), title=fr"True Signal Interaction Modes")
        set_axis_limits(ax,xlim=(0,11))
        add_text_box(ax, f"{num_ev_truth:.0f} Events", loc="upper left")
        add_text_box(ax, f"1/1001: QE\n10: MEC\n3:DIS\n4:RES\n5:COH", loc="upper right")
        
        plt.tight_layout()
        if output_dir:
            fig.savefig(f'{output_dir}/pngs/truth_signal_ixn_mode.png', dpi=300, bbox_inches='tight')
            with PdfPages(f'{output_dir}/pdfs/truth_signal_ixn_mode.pdf') as pdf:
                pdf.savefig(fig, dpi=300, bbox_inches='tight')
        return fig

    def plot_all(self, output_dir: str = "/outputs/truth/") -> None:
        """Generate all TruthPlotter plots and save to output directory."""
        from pathlib import Path
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        Path(output_dir+"/pngs").mkdir(parents=True, exist_ok=True)
        Path(output_dir+"/pdfs").mkdir(parents=True, exist_ok=True)
        
        print(f"Generating Truth plots in {output_dir}...")
        fig1 = self.plot_ixn_above_ke_threshold_per_spill(output_dir)
        plt.close(fig1)  # Close figure to save memory
        print("  ✓ Interactions above KE threshold per spill")
        fig2 = self.plot_pre_post_mx2_muon_angle(output_dir)
        plt.close(fig2)  # Close figure to save memory
        print("  ✓ Muon angle pre- and post-Mx2 Truth Cut")
        fig3 = self.plot_post_mx2_muon_energy(output_dir)
        plt.close(fig3)  # Close figure to save memory
        print("  ✓ Muon energy post-Mx2 Truth Cut")
        fig4 = self.plot_selected_ixn_nu_energy(output_dir)
        plt.close(fig4)  # Close figure to save memory
        print("  ✓ Neutrino energy")
        fig5 = self.plot_pre_mx2_pre1pi0_primary_pi0_mult(output_dir)
        plt.close(fig5)  # Close figure to save memory
        print("  ✓ Primary pi0 multiplicity pre-pi0 cut, pre-Mx2 Truth Cut")
        fig6 = self.plot_pre_post_mx2_sec_pi0_mult(output_dir)
        plt.close(fig6)  # Close figure to save memory
        print("  ✓ Secondary pi0 multiplicity pre- and post-Mx2 Truth Cut")
        fig7 = self.plot_shower_multiplicity_pre_mx2(output_dir)
        plt.close(fig7)  # Close figure to save memory
        print("  ✓ Shower multiplicity pre-Mx2 Truth Cut")
        fig8 = self.plot_shower_multiplicity_post_mx2(output_dir)
        plt.close(fig8)  # Close figure to save memory
        print("  ✓ Shower multiplicity post-Mx2 Truth Cut")
        fig9 = self.plot_true_signal_ixn_mode(output_dir)
        plt.close(fig9)  # Close figure to save memory
        print("  ✓ Interaction Mode")
        fig10 = self.plot_pre_mx2_muon_angle_numu_numubar(output_dir)
        plt.close(fig10)  # Close figure to save memory
        print("  ✓ Muon angle with beam pre-Mx2 Cheat Cut (numu vs. numbar)")

