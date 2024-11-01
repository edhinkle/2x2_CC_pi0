################################################################################
##                                                                            ##
##    CONTAINS: Script to create plots describing pi0 in signal               ##
##              events using a pi0 dictionary created using methods in        ##
##              /truth_studies/selections/dictionary_defs.py                  ##
##              and a scale factor for scaling event counts to those expected ##
##              with 1.2e19 POT. Adapted from https://github.com/edhinkle/    ##
##              mesonless_numubarCC/blob/main/truth_kinematics/plotting/      ##
##              plot_signal_hadrons.py                                        ##
##                                                                            ##
################################################################################

import matplotlib.pyplot as plt
import numpy as np
import sys
import math
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.axes import Axes
sys.path.append('../../')
import truth_studies.util.geometry_defs as geo_defs
import truth_studies.util.particlePDG_defs as pdg_defs

# Rasterize plots 
_old_axes_init = Axes.__init__
def _new_axes_init(self, *a, **kw):
    _old_axes_init(self, *a, **kw)
    # https://matplotlib.org/stable/gallery/misc/zorder_demo.html
    # 3 => leave text and legends vectorized
    self.set_rasterization_zorder(3)
def rasterize_plots():
    Axes.__init__ = _new_axes_init
def vectorize_plots():
    Axes.__init__ = _old_axes_init


# PLOT: pi0 kinematics
#       sig_bkg is an int such that 0 == signal, 1 == 'dirt' backgrounds, 2 == 'beam' backgrounds
def plot_pi0s(d, scale_factor, sig_bkg = 0):
    
    rasterize_plots()
    # DEFINE: Plotting pi0 kinematics for signal or background events
    sample_type = ''
    sample_title = ''
    if sig_bkg == 0: 
        sample_type = 'single_pi0_signal'
        sample_title = 'Single '+r'$\pi^0$ Signal'
    elif sig_bkg == 1:
        sample_type = 'dirt_bkg'
        sample_title = 'Dirt Background'
    elif sig_bkg == 2:
        sample_type = 'beam_bkg'
        sample_title = 'Beam Background'
    else: 
        return "Error: plot_pi0 function given undefined signal/background definition"
    output_pdf_name = sample_type+"_events_pi0_kinematics.pdf"
    print("Output PDF name: ", output_pdf_name)
    # put file in this directory for now
    with PdfPages(output_pdf_name, keep_empty=False) as output:   

        # PLOT: truth-level pi0 start momentum 
        fig, ax = plt.subplots(figsize=(6,4))
        data = np.array([np.sqrt(np.sum(np.array(d[key]['pi0_start_mom'])**2)) for key in d.keys()])/1000.
        counts, bins =np.histogram(data, bins=np.linspace(0,8,33))
        ax.hist(bins[:-1], bins=bins, weights = counts*scale_factor, histtype='stepfilled', alpha=0.5, color='blue')
        ax.set_xlabel(r"$\pi^0$ Initial Momentum [GeV/c]")
        ax.set_ylabel("Count / 0.25 GeV/c")
        output.savefig(fig)
        plt.close(fig)
        
        ## PLOT: total visible energy + contained visible energy
        #fig0, ax0 = plt.subplots(figsize=(8,4))
        #data0tot = np.array([d[key]['total_edep'] for key in d.keys()])
        ##data0cont = np.array([d[key]['contained_edep'] for key in d.keys()])
        #counts0tot, bins0tot = np.histogram(data0tot, bins=np.linspace(0,400,20))
        ##counts0cont, bins0cont = np.histogram(data0cont, bins=np.linspace(0,400,20))
        #ax0.hist(bins0tot[:-1], bins=bins0tot, weights = counts0tot*scale_factor, label='Total', histtype='step')
        ##ax0.hist(bins0cont[:-1], bins=bins0cont, weights = counts0cont*scale_factor, label='Contained',histtype='step', linestyle='--')
        #ax0.set_xlabel('Total Visible Muon Energy [MeV]')
        #ax0.set_ylabel('Count / 20 MeV')
        ##ax0.set_title(r'Muon Energy')
        ##ax0.legend()s
        #ax0.set_yscale('log')
        #ax0.grid(True)
        #output.savefig(fig0)
        #plt.close(fig0)

    return