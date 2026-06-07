#plotanapy/plotting/plotters/truth_match_plotter.py
"""
Plotter for truth-matched histograms in CC1pi0 analysis.
 
For each reconstructed event, these histograms reflect the best truth-matched
interaction. Organizes plotting into logical groupings:
- Track PDG & kinematics (mx2 primary/secondary track PDG, energy, angle)
- Muon kinematics (cos(angle), energy, response matrices)
- Vertex resolution (reco-truth differences in x, y, z)
- Shower & particle multiplicities (pi0, showers, electrons, photons)
- Multiplicity response matrices (truth vs reco)
- Signal/background breakdown by interaction mode (CCQE, CCMEC, CCDIS, CCRES, CCCOH, NC, Rock, CC bkg)
- Interaction mode & origin flags
"""
 
from typing import Optional, Dict
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from .base_plotter import HistogramPlotter
from ..styles import get_color, get_colors, configure_axes
from ..utils import add_colorbar, get_axis_label, add_text_box, set_axis_limits
 
 
class TruthMatchedPlotter(HistogramPlotter):
    """
    Plotter for truth-matched histograms from truth-matching step of CC1pi0 analysis.
 
    For each reconstructed selected event, a best-matched truth interaction is found.
    These histograms show distributions of truth and reco quantities for those matched events.
 
    Grouped plotting methods:
    - Mx2 match Track PDG & kinematics       : plot_mx2_track_kinematics
    - Muon in ixn kinematics                 : plot_muon_kinematics
    - Response matrices (ixn muon/mx2 match) : plot_response_matrices
    - Vertex resolution                      : plot_vertex_resolution
    - Particle (pi0/shower) multiplicities   : plot_multiplicities
    - Multiplicity response (true/reco)      : plot_multiplicity_response
    - Signal event pi0/shower breakdown      : plot_signal_multiplicity_by_mode
    - Background breakdown                   : plot_background_multiplicity_by_mode
    - Interaction mode & origin              : plot_interaction_mode_and_origin
    """
 
    def __init__(self, hist_dict, sel_config, beam_config, det_config):
        self.hist_dict = hist_dict
        self.sel_config = sel_config
        self.beam_config = beam_config
        self.det_config = det_config
 
 
    # -------------------------------------------------------------------------
    # Track PDG & Kinematics (Matched Mx2 track)
    # -------------------------------------------------------------------------
 
    def plot_mx2_track_kinematics_and_pdg(self, output_dir: Optional[str] = None) -> plt.Figure:
        """
        Plot PDG code, energy, and cos(angle) for primary and secondary mx2 tracks.
        Histograms:
            - h_truthMatchPrim_mx2TrackPDG
            - h_truthMatchSec_mx2TrackPDG
            - h_truthMatchPrim_mx2TrackE
            - h_truthMatchSec_mx2TrackE
            - h_truthMatchPrim_mx2TrackCosL
            - h_truthMatchSec_mx2TrackCosL
        """
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
 
        # PDG
        n = self.plot_1d('h_truthMatchPrim_mx2TrackPDG', normalized=True, ax=axes[0, 0], color=get_color("primary"))
        configure_axes(axes[0, 0], xlabel=get_axis_label('pdg'), ylabel=get_axis_label('frac'), title="PDG (True Primary)")
        add_text_box(axes[0, 0], f"{n:.0f} Events", loc="upper right")
        set_axis_limits(axes[0, 0],xlim=(-250,250))
 
        n = self.plot_1d('h_truthMatchSec_mx2TrackPDG', normalized=True, ax=axes[1, 0], color=get_color("secondary"))
        configure_axes(axes[1, 0], xlabel=get_axis_label('pdg'), ylabel=get_axis_label('frac'), title="PDG (True Secondary)")
        add_text_box(axes[1, 0], f"{n:.0f} Events", loc="upper right")
 
        # Energy
        n = self.plot_1d('h_truthMatchPrim_mx2TrackE', normalized=True, ax=axes[0, 1], color=get_color("primary"))
        configure_axes(axes[0, 1], xlabel=get_axis_label('energy', unit="GeV"), ylabel=get_axis_label('frac'), title="Energy (True Primary)")
        add_text_box(axes[0, 1], f"{n:.0f} Events", loc="upper right")
 
        n = self.plot_1d('h_truthMatchSec_mx2TrackE', normalized=True, ax=axes[1, 1], color=get_color("secondary"))
        configure_axes(axes[1, 1], xlabel=get_axis_label('energy', unit="GeV"), ylabel=get_axis_label('frac'), title="Energy (True Secondary)")
        add_text_box(axes[1, 1], f"{n:.0f} Events", loc="upper right")
 
        # CosL
        n = self.plot_1d('h_truthMatchPrim_mx2TrackCosL', normalized=True, ax=axes[0, 2], color=get_color("primary"))
        configure_axes(axes[0, 2], xlabel=get_axis_label('cosL'), ylabel=get_axis_label('frac'), title="cos(θ), Beam Angle (True Primary)")
        add_text_box(axes[0, 2], f"{n:.0f} Events", loc="upper left")
 
        n = self.plot_1d('h_truthMatchSec_mx2TrackCosL', normalized=True, ax=axes[1, 2], color=get_color("secondary"))
        configure_axes(axes[1, 2], xlabel=get_axis_label('cosL'), ylabel=get_axis_label('frac'), title="cos(θ), Beam Angle (True Secondary)")
        add_text_box(axes[1, 2], f"{n:.0f} Events", loc="upper left")

        fig.suptitle("Mx2 Track Match Truth Match PDG and Kinematics")
 
        plt.tight_layout()
        if output_dir:
            fig.savefig(f'{output_dir}/pngs/truthmatch_mx2_match_track_pdg_and_kinematics.png', dpi=300, bbox_inches='tight')
            with PdfPages(f'{output_dir}/pdfs/truthmatch_mx2_match_track_pdg_and_kinematics.pdf') as pdf:
                pdf.savefig(fig, dpi=300, bbox_inches='tight')
        return fig
 
 
    # -------------------------------------------------------------------------
    # Muon Kinematics for muon in truth match ixn
    # -------------------------------------------------------------------------
 
    def plot_truth_match_muon_kinematics(self, output_dir: Optional[str] = None) -> plt.Figure:
        """
        Plot truth-matched muon cos(angle) and energy distributions.
        Histograms:
            - h_truthMatchIxn_MuonCosL
            - h_truthMatchIxn_MuonCosL_unbinned
            - h_truthMatchIxn_MuonCosL_zoomOut
            - h_truthMatchIxn_MuonE
        """
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        cosThetaCut = self.sel_config["cosThetaCut"]
 
        n = self.plot_1d('h_truthMatchIxn_MuonCosL', normalized=True, ax=axes[0, 0], color=get_color("muon"))
        configure_axes(axes[0, 0], xlabel=get_axis_label('cosL'), ylabel=get_axis_label('frac'), title=f"cos(θ), Beam Angle (>{cosThetaCut})")
        add_text_box(axes[0, 0], f"{n:.0f} Events", loc="upper left")
 
        n = self.plot_1d('h_truthMatchIxn_MuonCosL_unbinned', normalized=True, ax=axes[0, 1], color=get_color("muon"))
        configure_axes(axes[0, 1], xlabel=get_axis_label('cosL'), ylabel=get_axis_label('frac'), title="cos(θ), Beam Angle (Unbinned)")
        add_text_box(axes[0, 1], f"{n:.0f} Events", loc="upper left")
 
        n = self.plot_1d('h_truthMatchIxn_MuonCosL_zoomOut', normalized=True, ax=axes[1, 0], color=get_color("muon"))
        configure_axes(axes[1, 0], xlabel=get_axis_label('cosL'), ylabel=get_axis_label('frac'), title="cos(θ), Beam Angle (>0.8)")
        add_text_box(axes[1, 0], f"{n:.0f} Events", loc="upper left")
 
        n = self.plot_1d('h_truthMatchIxn_MuonE', normalized=True, ax=axes[1, 1], color=get_color("muon"))
        configure_axes(axes[1, 1], xlabel=get_axis_label('energy', unit="GeV"), ylabel=get_axis_label('frac'), title="Energy")
        add_text_box(axes[1, 1], f"{n:.0f} Events", loc="upper right")

        fig.suptitle("Truth Matched Interaction Muon Kinematics")
 
        plt.tight_layout()
        if output_dir:
            fig.savefig(f'{output_dir}/pngs/truthmatch_muon_kinematics.png', dpi=300, bbox_inches='tight')
            with PdfPages(f'{output_dir}/pdfs/truthmatch_muon_kinematics.pdf') as pdf:
                pdf.savefig(fig, dpi=300, bbox_inches='tight')
        return fig
 
 
    # -------------------------------------------------------------------------
    # Response Matrices
    # -------------------------------------------------------------------------
 
    def plot_muon_response_matrices(self, output_dir: Optional[str] = None) -> plt.Figure:
        """
        Plot reco vs. truth response matrices for matched ixn muon and mx2 match track cos(angle).
        Histograms:
            - h_truthMatchIxn_Reco_responseMuonCosL
            - h_truthMatchIxn_Reco_responseMx2MatchPrimCosL
            - h_truthMatchIxn_Reco_responseMx2MatchSecCosL
        """
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
 
        mesh, n = self.plot_2d('h_truthMatchIxn_Reco_responseMuonCosL', normalized=True, ax=axes[0], cmap='viridis')
        configure_axes(axes[0], xlabel="Reco "+get_axis_label('cosL'), ylabel="Truth Matched"+get_axis_label('cosL'), title="Truth Match Ixn Muon vs. Reco")
        add_colorbar(mesh, axes[0], label=get_axis_label('frac'))
        add_text_box(axes[0], f"{n:.0f} Events", loc="upper left")
 
        mesh, n = self.plot_2d('h_truthMatchIxn_Reco_responseMx2MatchPrimCosL', normalized=True, ax=axes[1], cmap='viridis')
        configure_axes(axes[1], xlabel="Reco "+get_axis_label('cosL'), ylabel="Truth Matched"+get_axis_label('cosL'), title="Mx2 Match Truth Match (Primary) vs. Reco")
        add_colorbar(mesh, axes[1], label=get_axis_label('frac'))
        add_text_box(axes[1], f"{n:.0f} Events", loc="upper left")
 
        mesh, n = self.plot_2d('h_truthMatchIxn_Reco_responseMx2MatchSecCosL', normalized=True, ax=axes[2], cmap='viridis')
        configure_axes(axes[2], xlabel="Reco "+get_axis_label('cosL'), ylabel="Truth Matched"+get_axis_label('cosL'), title="Mx2 Match Truth Match (Secondary) vs. Reco")
        add_colorbar(mesh, axes[2], label=get_axis_label('frac'))
        add_text_box(axes[2], f"{n:.0f} Events", loc="upper left")

        fig.suptitle(r"cos($\theta$) (Beam/muon angle) Reco/Truth Match Response")
 
        plt.tight_layout()
        if output_dir:
            fig.savefig(f'{output_dir}/pngs/truthmatch_ixnmuon_mx2match_kinematics_response_matrices_with_reco.png', dpi=300, bbox_inches='tight')
            with PdfPages(f'{output_dir}/pdfs/truthmatch_ixnmuon_mx2match_kinematics_response_matrices_with_reco.pdf') as pdf:
                pdf.savefig(fig, dpi=300, bbox_inches='tight')
        return fig
 
 
    # -------------------------------------------------------------------------
    # Vertex and Mx2 Track Match Start Resolution
    # -------------------------------------------------------------------------
 
    def plot_vertex_mx2_track_match_start_resolution(self, output_dir: Optional[str] = None) -> plt.Figure:
        """
        Plot reco-truth vertex differences for primary and secondary tracks, plus overall.
        Histograms:
            - h_truthMatchPrim_mx2TrackRecoDiffLArStartX/Y/Z
            - h_truthMatchSec_mx2TrackRecoDiffLArStartX/Y/Z
            - h_truthMatch_diffTruthRecoVertex
        """
        fig, axes = plt.subplots(3, 3, figsize=(15, 13))
 
        for col, coord in enumerate(['X', 'Y', 'Z']):
            n = self.plot_1d(f'h_truthMatchPrim_mx2TrackRecoDiffLArStart{coord}', normalized=True,
                             ax=axes[0, col], color=get_color("primary"))
            configure_axes(axes[0, col], xlabel=get_axis_label(f'vertex_{coord.lower()}'),
                           ylabel=get_axis_label('frac'), title=f"Mx2 Match Track Δ{coord} (Truth Primary-Reco)")
            add_text_box(axes[0, col], f"{n:.0f} Events", loc="upper right")
 
            n = self.plot_1d(f'h_truthMatchSec_mx2TrackRecoDiffLArStart{coord}', normalized=True,
                             ax=axes[1, col], color=get_color("secondary"))
            configure_axes(axes[1, col], xlabel=get_axis_label(f'vertex_{coord.lower()}'),
                           ylabel=get_axis_label('frac'), title=f"Mx2 Match Track Δ{coord} (Truth Secondary-Reco)")
            add_text_box(axes[1, col], f"{n:.0f} Events", loc="upper right")
 
        # Overall vertex diff (spans bottom row center)
        n = self.plot_1d('h_truthMatch_diffTruthRecoVertex', normalized=True, ax=axes[2, 1], color=get_colors(1))
        configure_axes(axes[2, 1], xlabel=get_axis_label('vertex_diff'), ylabel=get_axis_label('frac'),
                       title="Overall Reco−Truth Vertex Δ")
        add_text_box(axes[2, 1], f"{n:.0f} Events", loc="upper right")
 
        # Hide unused bottom-row axes
        axes[2, 0].set_visible(False)
        axes[2, 2].set_visible(False)
 
        plt.tight_layout()
        if output_dir:
            fig.savefig(f'{output_dir}/pngs/truthmatch_mx2trackstart_and_vertex_resolution.png', dpi=300, bbox_inches='tight')
            with PdfPages(f'{output_dir}/pdfs/truthmatch_mx2trackstart_and_vertex_resolution.pdf') as pdf:
                pdf.savefig(fig, dpi=300, bbox_inches='tight')
        return fig
 
 
    # -------------------------------------------------------------------------
    # Particle (pi0, shower) Multiplicities (Primary and Secondary)
    # -------------------------------------------------------------------------
 
    def plot_truth_matched_ixn_pi0_shower_multiplicities(self, output_dir: Optional[str] = None) -> plt.Figure:
        """
        Plot pi0, shower, electron, and photon multiplicities (primary and secondary).
        Histograms:
            - h_truthMatchIxn_PrimPi0Multiplicity / SecPi0Multiplicity
            - h_truthMatchIxn_Prim/SecShowerMultiplicity
            - h_truthMatchIxn_Prim/SecElectronMultiplicity
            - h_truthMatchIxn_Prim/SecPhotonMultiplicity
        """
        fig, axes = plt.subplots(4, 2, figsize=(12, 18))
 
        particle_rows = [
            ('Pi0',      get_color("pion")),
            ('Shower',   get_color("primary")),
            ('Electron', get_color("electron")),
            ('Photon',   get_color("photon")),
        ]
 
        for row, (particle, color) in enumerate(particle_rows):
            for col, level in enumerate(['Prim', 'Sec']):
                key = f'h_truthMatchIxn_{level}{particle}Multiplicity'
                n = self.plot_1d(key, normalized=True, ax=axes[row, col], color=color)
                if (particle == 'Pi0'):
                    configure_axes(axes[row, col], xlabel=get_axis_label('multiplicity'),
                               ylabel=get_axis_label('frac'),
                               title=fr"Truth-Matched Ixn {'Primary' if level == 'Prim' else 'Secondary'} $\pi^0$ Mult")
                else: 
                    configure_axes(axes[row, col], xlabel=get_axis_label('multiplicity'),
                               ylabel=get_axis_label('frac'),
                               title=fr"Truth-Matched Ixn {'Primary' if level == 'Prim' else 'Secondary'} {particle} Mult")
                add_text_box(axes[row, col], f"{n:.0f} Events", loc="upper right")
 
        plt.tight_layout()
        if output_dir:
            fig.savefig(f'{output_dir}/pngs/truthmatch_ixn_pi0_shower_multiplicities.png', dpi=300, bbox_inches='tight')
            with PdfPages(f'{output_dir}/pdfs/truthmatch_ixn_pi0_shower_multiplicities.pdf') as pdf:
                pdf.savefig(fig, dpi=300, bbox_inches='tight')
        return fig
 
 
    # -------------------------------------------------------------------------
    # Pi0/Shower Multiplicity Response Matrices
    # -------------------------------------------------------------------------
 
    def plot_reco_truthmatch_shower_multiplicity_response(self, output_dir: Optional[str] = None) -> plt.Figure:
        """
        Plot reco vs. truth multiplicity response matrices for showers, electrons, photons.
        Histograms:
            - h_truthMatchIxn_Reco_Prim/SecShowerResponseMult
            - h_truthMatchIxn_Reco_Prim/SecElectronResponseMult
            - h_truthMatchIxn_Reco_Prim/SecPhotonResponseMult
        """
        fig, axes = plt.subplots(3, 2, figsize=(12, 15))
 
        particles = [
            ('Shower',   get_color("primary")),
            ('Electron', get_color("electron")),
            ('Photon',   get_color("photon")),
        ]
 
        for row, (particle, _) in enumerate(particles):
            for col, level in enumerate(['Prim', 'Sec']):
                key = f'h_truthMatchIxn_Reco_{level}{particle}ResponseMult'
                mesh, n = self.plot_2d(key, normalized=True, ax=axes[row, col], cmap='viridis')
                configure_axes(axes[row, col],
                               xlabel="Reco "+get_axis_label('multiplicity'),
                               ylabel="Truth Match "+get_axis_label('multiplicity'),
                               title=f"Multiplicity Response: {'Primary' if level == 'Prim' else 'Secondary'} {particle}s")
                add_colorbar(mesh, axes[row, col], label=get_axis_label('frac'))
                add_text_box(axes[row, col], f"{n:.0f} Events", loc="upper left")
 
        plt.tight_layout()
        if output_dir:
            fig.savefig(f'{output_dir}/pngs/truthmatch_reco_truthmatch_shower_multiplicity_response.png', dpi=300, bbox_inches='tight')
            with PdfPages(f'{output_dir}/pdfs/truthmatch_reco_truthmatch_shower_multiplicity_response.pdf') as pdf:
                pdf.savefig(fig, dpi=300, bbox_inches='tight')
        return fig
 
 
    # -------------------------------------------------------------------------
    # Signal Multiplicity Breakdown by Interaction Mode
    # -------------------------------------------------------------------------
 
    def plot_signal_reco_multiplicity_by_truth_match_ixn_mode(self, output_dir: Optional[str] = None) -> Dict[str, plt.Figure]:
        """
        Plot reco shower/electron/photon multiplicities for reco ixns matched to true signal events,
        broken down by interaction mode: inclusive signal, CCQE, CCMEC, CCDIS, CCRES, CCCOH.
        Histograms (pattern):
            h_truthMatchIxn_TruthSignal{Mode}_Reco{Level}{Particle}Mult
        where Mode ∈ {'' (inclusive), CCQE, CCMEC, CCDIS, CCRES, CCCOH}
              Level ∈ {Prim, Sec}
              Particle ∈ {Shower, Electron, Photon}
        """
        modes = [
            ('',      'Signal (All)'),
            ('CCQE',  'Signal CCQE'),
            ('CCMEC', 'Signal CCMEC'),
            ('CCDIS', 'Signal CCDIS'),
            ('CCRES', 'Signal CCRES'),
            ('CCCOH', 'Signal CCCOH'),
        ]
        particles = ['Shower', 'Electron', 'Photon']
        levels    = ['Prim', 'Sec']
 
        # One figure per particle type: rows = modes, cols = Prim/Sec
        figs = {}
        for particle in particles:
            fig, axes = plt.subplots(len(modes), 2, figsize=(12, 4 * len(modes)))
            figs[particle] = fig
 
            for row, (mode_key, mode_label) in enumerate(modes):
                for col, level in enumerate(levels):
                    key = f'h_truthMatchIxn_TruthSignal{mode_key}_Reco{level}{particle}Mult'
                    n = self.plot_1d(key, normalized=True, ax=axes[row, col],
                                     color=get_color(particle.lower()))
                    configure_axes(axes[row, col],
                                   xlabel="Reco "+get_axis_label('multiplicity'),
                                   ylabel=get_axis_label('frac'),
                                   title=f"{mode_label}: {'Prim' if level == 'Prim' else 'Sec'} {particle}s")
                    add_text_box(axes[row, col], f"{n:.0f} Events", loc="upper right")

            fig.suptitle(f"Reco {particle} Multiplicity for Truth-Matched Signal by Ixn Mode")
 
            plt.tight_layout()
            if output_dir:
                fname = f'truthmatch_signal_{particle.lower()}_reco_mult_by_ixn_mode'
                fig.savefig(f'{output_dir}/pngs/{fname}.png', dpi=300, bbox_inches='tight')
                with PdfPages(f'{output_dir}/pdfs/{fname}.pdf') as pdf:
                    pdf.savefig(fig, dpi=300, bbox_inches='tight')
 
        return figs   # dict keyed by particle name
 
 
    # -------------------------------------------------------------------------
    # Background Multiplicity Breakdown by Interaction Mode
    # -------------------------------------------------------------------------
 
    def plot_background_reco_multiplicity_by_truth_match_ixn_mode(self, output_dir: Optional[str] = None) -> Dict[str, plt.Figure]:
        """
        Plot reco shower/electron/photon multiplicities for truth-background events,
        broken down by background type: all bkg, NC, Rock, CC bkg.
        Histograms (pattern):
            h_truthMatchIxn_TruthBkg{Type}_Reco{Level}{Particle}Mult
        where Type ∈ {'' (inclusive), NC, ROCK, CC}
              Level ∈ {Prim, Sec}
              Particle ∈ {Shower, Electron, Photon}
        """
        bkg_types = [
            ('',     'All Background'),
            ('NC',   'NC Background'),
            ('ROCK', 'Rock Background'),
            ('CC',   'CC Background'),
        ]
        particles = ['Shower', 'Electron', 'Photon']
        levels    = ['Prim', 'Sec']
 
        figs = {}
        for particle in particles:
            fig, axes = plt.subplots(len(bkg_types), 2, figsize=(12, 4 * len(bkg_types)))
            figs[particle] = fig
 
            for row, (bkg_key, bkg_label) in enumerate(bkg_types):
                for col, level in enumerate(levels):
                    key = f'h_truthMatchIxn_TruthBkg{bkg_key}_Reco{level}{particle}Mult'
                    n = self.plot_1d(key, normalized=True, ax=axes[row, col],
                                     color=get_color(particle.lower()))
                    configure_axes(axes[row, col],
                                   xlabel="Reco "+get_axis_label('multiplicity'),
                                   ylabel=get_axis_label('frac'),
                                   title=f"{bkg_label}: {'Prim' if level == 'Prim' else 'Sec'} {particle}s")
                    add_text_box(axes[row, col], f"{n:.0f} Events", loc="upper right")

            fig.suptitle(f"Reco {particle} Multiplicity for Truth-Matched Backgroung by Ixn Type")
 
            plt.tight_layout()
            if output_dir:
                fname = f'truthmatch_bkg_{particle.lower()}_reco_mult_by_mode'
                fig.savefig(f'{output_dir}/pngs/{fname}.png', dpi=300, bbox_inches='tight')
                with PdfPages(f'{output_dir}/pdfs/{fname}.pdf') as pdf:
                    pdf.savefig(fig, dpi=300, bbox_inches='tight')
 
        return figs
 
 
    # -------------------------------------------------------------------------
    # Interaction Mode & Origin Flags (All, Non-Rock and Rock Bkgs)
    # -------------------------------------------------------------------------
 
    def plot_truth_match_interaction_mode_and_origin(self, output_dir: Optional[str] = None) -> plt.Figure:
        """
        Plot interaction mode, rock flag, and neutrino PDG for signal and background subsets.
        Histograms:
            - h_truthMatchIxn_IxnMode
            - h_truthMatchIxn_IsRock
            - h_truthMatchIxn_TruthBkgROCK_isCC
            - h_truthMatchIxn_TruthBkgROCK_nuPDG
            - h_truthMatchIxn_TruthBkgNONROCK_isCC
            - h_truthMatchIxn_TruthBkgNONROCK_nuPDG
        """
        fig, axes = plt.subplots(3, 2, figsize=(10, 15))
        colors = get_colors(4)
 
        n = self.plot_1d('h_truthMatchIxn_IxnMode', normalized=True, ax=axes[0, 0], color=colors[0])
        configure_axes(axes[0, 0], xlabel="GENIE Interaction Mode", ylabel=get_axis_label('frac'),
                       title="Truth-Matched Interaction Mode")
        add_text_box(axes[0, 0], f"{n:.0f} Events", loc="upper left")
        add_text_box(axes[0, 0], f"1/1001: QE\n10: MEC\n3:DIS\n4:RES\n5:COH", loc="upper right")
        set_axis_limits(axes[0, 0],xlim=(0,11))
 
        n = self.plot_1d('h_truthMatchIxn_IsRock', normalized=True, ax=axes[0, 1], color=colors[0])
        configure_axes(axes[0, 1], xlabel="Rock = 1, Non-Rock=0", ylabel=get_axis_label('frac'),
                       title="Truth-Matched Rock Interaction?")
        add_text_box(axes[0, 1], f"{n:.0f} Events", loc="upper right")
 
 
        n = self.plot_1d('h_truthMatchIxn_TruthBkgROCK_isCC', normalized=True, ax=axes[1, 0], color=colors[1])
        configure_axes(axes[1, 0], xlabel="CC=1, NC=0", ylabel=get_axis_label('frac'),
                       title="Truth-Matched Background Rock CC vs. NC")
        add_text_box(axes[1, 0], f"{n:.0f} Events", loc="upper right")
 
        n = self.plot_1d('h_truthMatchIxn_TruthBkgROCK_nuPDG', normalized=True, ax=axes[1, 1], color=colors[1])
        configure_axes(axes[1, 1], xlabel=get_axis_label('pdg'), ylabel=get_axis_label('frac'),
                       title=r"Truth-Matched Background Rock $\nu$ PDG")
        add_text_box(axes[1, 1], f"{n:.0f} Events", loc="upper right")
 
        n = self.plot_1d('h_truthMatchIxn_TruthBkgNONROCK_isCC', normalized=True, ax=axes[2, 0], color=colors[3])
        configure_axes(axes[2, 0], xlabel="CC=1, NC=0", ylabel=get_axis_label('frac'),
                       title="Truth-Matched Background Non-Rock CC vs. NC")
        add_text_box(axes[2, 0], f"{n:.0f} Events", loc="upper right")

        n = self.plot_1d('h_truthMatchIxn_TruthBkgNONROCK_nuPDG', normalized=True, ax=axes[2, 1], color=colors[3])
        configure_axes(axes[2, 1], xlabel=get_axis_label('pdg'), ylabel=get_axis_label('frac'),
                       title=r"Truth-Matched Background Non-Rock $\nu$ PDG")
        add_text_box(axes[2, 1], f"{n:.0f} Events", loc="upper right")

 
        plt.tight_layout()
        if output_dir:
            fig.savefig(f'{output_dir}/pngs/truthmatch_rock_interaction_mode_origin.png', dpi=300, bbox_inches='tight')
            with PdfPages(f'{output_dir}/pdfs/truthmatch_rock_interaction_mode_origin.pdf') as pdf:
                pdf.savefig(fig, dpi=300, bbox_inches='tight')
        return fig
 
 
    # -------------------------------------------------------------------------
    # plot_all
    # -------------------------------------------------------------------------
 
    def plot_all(self, output_dir: str = "/outputs/truth_matched_plots/") -> None:
        """Generate all TruthMatchedPlotter plots and save to output directory."""
        from pathlib import Path
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        Path(output_dir + "/pngs").mkdir(parents=True, exist_ok=True)
        Path(output_dir + "/pdfs").mkdir(parents=True, exist_ok=True)
 
        print(f"Generating TruthMatched plots in {output_dir}...")
 
        fig = self.plot_mx2_track_kinematics_and_pdg(output_dir)
        plt.close(fig)
        print("  ✓ Mx2 track match kinematics and PDG")
 
        fig = self.plot_truth_match_muon_kinematics(output_dir)
        plt.close(fig)
        print("  ✓ Muon kinematics in truth matched interaction")
 
        fig = self.plot_muon_response_matrices(output_dir)
        plt.close(fig)
        print("  ✓ Muon cos(theta) Response matrices")
 
        fig = self.plot_vertex_mx2_track_match_start_resolution(output_dir)
        plt.close(fig)
        print("  ✓ Vertex and Mx2 Track match start point resolution")
 
        fig = self.plot_truth_matched_ixn_pi0_shower_multiplicities(output_dir)
        plt.close(fig)
        print("  ✓ Truth matched interaction pi0 and shower multiplicities")
 
        fig = self.plot_reco_truthmatch_shower_multiplicity_response(output_dir)
        plt.close(fig)
        print("  ✓ Shower reco vs. truth match multiplicity response matrices")
 
        figs = self.plot_signal_reco_multiplicity_by_truth_match_ixn_mode(output_dir)
        for fig in figs.values():
            plt.close(fig)
        print("  ✓ Signal truth matched reco multiplicities by interaction mode")
 
        figs = self.plot_background_reco_multiplicity_by_truth_match_ixn_mode(output_dir)
        for fig in figs.values():
            plt.close(fig)
        print("  ✓ Background truth matched reco multiplicities by type of background")
 
        fig = self.plot_truth_match_interaction_mode_and_origin(output_dir)
        plt.close(fig)
        print("  ✓ Interaction mode & origin for truth matched ixns including Rock, Non-Rock bkg categories")
 
        print(f"Done! All truth matched plots saved to {output_dir}")
 

