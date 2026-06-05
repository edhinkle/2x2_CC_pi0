# plotanapy/plotting/plotters/selection_plotter.py

from typing import Dict, List, Optional, Tuple, Any
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

from .base_plotter import HistogramPlotter
from ..styles import configure_axes, get_colors
from ..utils import add_text_box


class SelectionPlotter(HistogramPlotter):
    """
    Plotter for cutflow, efficiency and purity studies.

    Histograms expected:
        - Reco
        - Truth
        - RecoMatchedSignal
        - RecoMatchedSignalwithMx2
        - RecoMatchedValidVtx
    """

    def __init__(
        self,
        hist_dict,
        sel_config,
        beam_config,
        det_config,
        mcOnly
    ):
        self.hist_dict = hist_dict
        self.sel_config = sel_config
        self.beam_config = beam_config
        self.det_config = det_config
        self.mcOnly = True if mcOnly=='1' else False

    def plot_sel_truth(self, output_dir: Optional[str] = None) -> Tuple[plt.Figure, int]:
        """ Plot truth selection cutflow (Histogram == "Truth")"""

        fig, axes, true_cut_labels, true_cut_values = self.plot_cutflow("Truth")
        configure_axes(axes[0], ylabel="Events", title="CC1π⁰ Truth Selection")
        true_signal_events_total = np.full_like(true_cut_values, true_cut_values[-1])

        purity = np.divide(true_signal_events_total, true_cut_values,
                           out=np.zeros_like(true_cut_values, dtype=float),
                           where=true_cut_values > 0)

        axes[1].set_ylim(0, 1.05)
        axes[1].grid(alpha=0.3)

        true_cut_label_xvals = np.arange(len(true_cut_labels))+0.5
        axes[1].plot(true_cut_label_xvals, purity,marker="s",linewidth=2,label="Purity")
        configure_axes(axes[1], xlabel="Selection Cut", ylabel="Fraction")

        axes[1].legend()
        axes[0].set_xlim(0, len(true_cut_labels))
        axes[1].set_xticks(true_cut_label_xvals)
        axes[0].set_xticklabels(true_cut_labels, rotation=45, ha="right")

        # Plot pur, eff
        for i in range(len(purity)):

            axes[1].text(i+0.5, purity[i] + 0.03, f"{purity[i]:.2f}", ha="center", fontsize=12)

        plt.tight_layout()
        if output_dir:
            fig.savefig(f'{output_dir}/pngs/truth_selection.png', dpi=300, bbox_inches='tight')
            with PdfPages(f'{output_dir}/pdfs/truth_selection.pdf') as pdf:
                pdf.savefig(fig, dpi=300, bbox_inches='tight')

        return fig, true_cut_values[-1]

    def plot_sel_reco(self, truth_total_events: Optional[int] = None,
                      output_dir: Optional[str] = None) -> plt.Figure:
        """
        Plot truth selection cutflow. Histograms:
            - Reco
            - RecoMatchedSignal (if with MC)
            - RecoMatchedSignalwithMx2 (if with MC) -- not currently using
            - RecoMatchedValidVtx (if with MC) -- not currently using
        """
        colors = get_colors(6)
        fig, axes, reco_cut_labels, reco_cut_values = self.plot_cutflow("Reco", label="Reco Selected", colors=colors[5])
        configure_axes(axes[0], ylabel="Events", title="CC1π⁰ Reco Selection")


        axes[0].set_xlim(0, len(reco_cut_labels))

        if (truth_total_events != 0): 

            fig, axes, reco_matched_signal_cut_labels, reco_matched_signal_cut_values = self.plot_cutflow("RecoMatchedSignal", \
                                                                                                     label="Matched Signal Selected", \
                                                                                                      axes=axes, colors=colors[1])
            
            axes[0].legend()
            purity = np.divide(reco_matched_signal_cut_values, reco_cut_values,
                               out=np.zeros_like(reco_cut_values, dtype=float),
                               where=reco_cut_values > 0)
            
            truth_signal_arr = np.full_like(reco_matched_signal_cut_values, truth_total_events)
            efficiency = np.divide(reco_matched_signal_cut_values, truth_signal_arr,
                               out=np.zeros_like(reco_matched_signal_cut_values, dtype=float))


            axes[1].set_ylim(0, 1.05)
            axes[1].grid(alpha=0.3)

            reco_cut_label_xvals = np.arange(len(reco_cut_labels))+0.5
            axes[1].plot(reco_cut_label_xvals, purity,marker="s",linewidth=2,label="Purity", color=colors[2])
            axes[1].plot(reco_cut_label_xvals, efficiency,marker="o",linewidth=2,label="Efficiency", color=colors[4])
            axes[1].set_xticks(reco_cut_label_xvals)
            axes[0].set_xticklabels(reco_cut_labels, rotation=45, ha="right")
            configure_axes(axes[1], ylabel="Fraction")

            axes[1].legend(loc="upper right")

            # Plot pur, eff
            for i in range(len(purity)):

                axes[1].text(i+0.5, purity[i] + 0.03, f"{purity[i]:.2f}", ha="center", fontsize=8)
                axes[1].text(i+0.5, efficiency[i] + 0.03, f"{efficiency[i]:.2f}", ha="center", fontsize=8)

            add_text_box(axes[0],(f"Truth: {truth_total_events:.0f}\n"
                                  f"Final Reco: {reco_cut_values[-1]:.0f}\n"
                                  f"Final Eff: {100*efficiency[-1]:.1f}%\n"
                                  f"Final Pur: {100*purity[-1]:.1f}%"), loc="upper right")

        plt.tight_layout()
        if output_dir:
            fig.savefig(f'{output_dir}/pngs/reco_selection.png', dpi=300, bbox_inches='tight')
            with PdfPages(f'{output_dir}/pdfs/reco_selection.pdf') as pdf:
                pdf.savefig(fig, dpi=300, bbox_inches='tight')

        return fig 






    def plot_all(self, output_dir: str = "/outputs/cutflow/") -> None:
        """Generate all SelectionPlotter plots and save to output directory."""
        from pathlib import Path
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        Path(output_dir+"/pngs").mkdir(parents=True, exist_ok=True)
        Path(output_dir+"/pdfs").mkdir(parents=True, exist_ok=True)
        
        print(f"Generating Selection plots in {output_dir}...")

        if self.mcOnly == True:
            fig1, true_signal_total = self.plot_sel_truth(output_dir)
            plt.close(fig1)  # Close figure to save memory
            print("  ✓ Truth Selection")
        else:
            print("   No MC Selection...")
            true_signal_total = 0

        fig2 = self.plot_sel_reco(true_signal_total, output_dir)
        plt.close(fig2)  # Close figure to save memory
        print("  ✓ Reco Selection")
        