################################################################################
##                                                                            ##
##    CONTAINS: Script to create plots describing pi0 in signal and bkg       ##
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
import pandas as pd
import itertools
import sys
import math
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.axes import Axes
import json
import argparse
sys.path.append('../../')
import truth_studies.util.geometry_defs as geo_defs
import truth_studies.util.particlePDG_defs as pdg_defs
import cmasher as cmr

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
def plot_sig_and_bkg_CC_1pi0(df, scale_factor):
    
    #rasterize_plots()
    # Get signal and background dataframes
    hist_labels = df['sig_bkg_label'].value_counts().index
    hist_labels = np.sort(hist_labels)
    list_of_sig_bkg_dfs = []
    for i in range(len(hist_labels)):

        sig_bkg_label = hist_labels[i]
        sig_bkg_df_i = df[df['sig_bkg_label'] == sig_bkg_label]
        list_of_sig_bkg_dfs.append(sig_bkg_df_i)
    #df_signal = df[df['sig']==True]
    #print("Number of signal events: ", len(df_signal)*scale_factor)
    #df_all_background = df[df['sig']==False]
    #print("Number of background events: ", len(df_all_background)*scale_factor)
    #df_bkg_1pi0 = df_all_background[df_all_background['pi0_mult']==1]
    #df_bkg_1pi0_with_charged_pions = df_bkg_1pi0[df_bkg_1pi0['primary_charged_pion_mult']>0]
    ##df_bkg_1pi0_with_charged_kaons = df_bkg_1pi0[df_bkg_1pi0['primary_charged_kaon_mult']>0]
    ##df_bkg_1pi0_no_charged_mesons = df_bkg_1pi0[df_bkg_1pi0['primary_charged_pion_mult']==0]
    #df_bkg_1pi0_no_charged_pions = df_bkg_1pi0[df_bkg_1pi0['primary_charged_pion_mult']==0]
    #df_bkg_other = df_all_background[df_all_background['pi0_mult']!=1]


    output_pdf_name = "signal_and_background_plots.pdf"
    print("Output PDF name: ", output_pdf_name)
    # put file in this directory for now
    with PdfPages(output_pdf_name, keep_empty=False) as output:   

        ## PLOT: truth-level pi0 multiplicity
        fig, ax = plt.subplots()
        hist_data = []
        hist_weights = []
        for i in range(len(list_of_sig_bkg_dfs)):
            data = np.array(list_of_sig_bkg_dfs[i]['pi0_mult'])
            counts, bins = np.histogram(data, bins=np.linspace(0,10,11) )
            hist_data.append(bins[:-1])
            hist_weights.append(counts*scale_factor)
        hist_colors = cmr.take_cmap_colors('nipy_spectral', len(hist_labels))
        ax.hist(hist_data, bins=bins, weights = hist_weights, histtype='stepfilled', alpha=0.9, color=hist_colors, label=hist_labels, stacked=True)
        ax.set_xlabel(r"$\pi^0$ Multiplicity")
        ax.set_ylabel(r"Count [Scaled to 1.2e19 POT] / $\pi^0$")
        plt.legend()
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        ## PLOT: truth-level muon multiplicity
        fig, ax = plt.subplots(figsize=(6,4))
        hist_data = []
        hist_weights = []
        for i in range(len(list_of_sig_bkg_dfs)):
            data = np.array(list_of_sig_bkg_dfs[i]['muon_mult'])
            counts, bins = np.histogram(data, bins=np.linspace(0,10,11) )
            hist_data.append(bins[:-1])
            hist_weights.append(counts*scale_factor)
        hist_colors = cmr.take_cmap_colors('nipy_spectral', len(hist_labels))
        ax.hist(hist_data, bins=bins, weights = hist_weights, histtype='stepfilled', alpha=0.9, color=hist_colors, label=hist_labels, stacked=True)
        ax.set_xlabel(r"Muon Multiplicity")
        ax.set_ylabel(r"Count [Scaled to 1.2e19 POT] / Muon")
        plt.legend()
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        ## PLOT: truth-level primary electron multiplicity
        fig, ax = plt.subplots()
        hist_data = []
        hist_weights = []
        for i in range(len(list_of_sig_bkg_dfs)):
            data = np.array(list_of_sig_bkg_dfs[i]['primary_electron_mult'])
            counts, bins = np.histogram(data, bins=np.linspace(0,10,11) )
            hist_data.append(bins[:-1])
            hist_weights.append(counts*scale_factor)
        hist_colors = cmr.take_cmap_colors('nipy_spectral', len(hist_labels))
        ax.hist(hist_data, bins=bins, weights = hist_weights, histtype='stepfilled', alpha=0.9, color=hist_colors, label=hist_labels, stacked=True)
        ax.set_xlabel(r"Primary Electron Multiplicity")
        ax.set_ylabel(r"Count [Scaled to 1.2e19 POT] / Electron")
        plt.legend()
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        ## PLOT: truth-level primary photon multiplicity
        fig, ax = plt.subplots()
        hist_data = []
        hist_weights = []
        for i in range(len(list_of_sig_bkg_dfs)):
            data = np.array(list_of_sig_bkg_dfs[i]['primary_photon_mult'])
            counts, bins = np.histogram(data, bins=np.linspace(0,10,11) )
            hist_data.append(bins[:-1])
            hist_weights.append(counts*scale_factor)
        hist_colors = cmr.take_cmap_colors('nipy_spectral', len(hist_labels))
        ax.hist(hist_data, bins=bins, weights = hist_weights, histtype='stepfilled', alpha=0.9, color=hist_colors, label=hist_labels, stacked=True)
        ax.set_xlabel(r"Primary Photon Multiplicity")
        ax.set_ylabel(r"Count [Scaled to 1.2e19 POT] / Photon")
        plt.legend()
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        ## PLOT: truth-level primary proton multiplicity
        fig, ax = plt.subplots()
        hist_data = []
        hist_weights = []
        for i in range(len(list_of_sig_bkg_dfs)):
            data = np.array(list_of_sig_bkg_dfs[i]['primary_proton_mult'])
            counts, bins = np.histogram(data, bins=np.linspace(0,10,11) )
            hist_data.append(bins[:-1])
            hist_weights.append(counts*scale_factor)
        hist_colors = cmr.take_cmap_colors('nipy_spectral', len(hist_labels))
        ax.hist(hist_data, bins=bins, weights = hist_weights, histtype='stepfilled', alpha=0.9, color=hist_colors, label=hist_labels, stacked=True)
        ax.set_xlabel(r"Primary Proton Multiplicity")
        ax.set_ylabel(r"Count [Scaled to 1.2e19 POT] / Proton")
        plt.legend()
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        ## PLOT: truth-level primary charged pion multiplicity
        fig, ax = plt.subplots()
        hist_data = []
        hist_weights = []
        for i in range(len(list_of_sig_bkg_dfs)):
            data = np.array(list_of_sig_bkg_dfs[i]['primary_charged_pion_mult'])
            counts, bins = np.histogram(data, bins=np.linspace(0,10,11) )
            hist_data.append(bins[:-1])
            hist_weights.append(counts*scale_factor)
        hist_colors = cmr.take_cmap_colors('nipy_spectral', len(hist_labels))
        ax.hist(hist_data, bins=bins, weights = hist_weights, histtype='stepfilled', alpha=0.9, color=hist_colors, label=hist_labels, stacked=True)
        ax.set_xlabel(r"Primary Charged Pion Multiplicity")
        ax.set_ylabel(r"Count [Scaled to 1.2e19 POT] / Charged Pion")
        plt.legend()
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)


        ## PLOT: truth-level primary charged kaon multiplicity
        fig, ax = plt.subplots()
        hist_data = []
        hist_weights = []
        for i in range(len(list_of_sig_bkg_dfs)):
            data = np.array(list_of_sig_bkg_dfs[i]['primary_charged_kaon_mult'])
            counts, bins = np.histogram(data, bins=np.linspace(0,10,11) )
            hist_data.append(bins[:-1])
            hist_weights.append(counts*scale_factor)
        hist_colors = cmr.take_cmap_colors('nipy_spectral', len(hist_labels))
        ax.hist(hist_data, bins=bins, weights = hist_weights, histtype='stepfilled', alpha=0.9, color=hist_colors, label=hist_labels, stacked=True)
        ax.set_xlabel(r"Primary Charged Kaon Multiplicity")
        ax.set_ylabel(r"Count [Scaled to 1.2e19 POT] / Charged Kaon")
        plt.legend()
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        ## PLOT: truth-level muon angle
        fig, ax = plt.subplots()
        hist_data = []
        hist_weights = []
        for i in range(len(list_of_sig_bkg_dfs)):
            data = np.array(list_of_sig_bkg_dfs[i]['muon_angle'])
            counts, bins = np.histogram(data, bins=np.linspace(0,2,21) )
            hist_data.append(bins[:-1])
            hist_weights.append(counts*scale_factor)
        hist_colors = cmr.take_cmap_colors('nipy_spectral', len(hist_labels))
        ax.hist(hist_data, bins=bins, weights = hist_weights, histtype='stepfilled', alpha=0.9, color=hist_colors, label=hist_labels, stacked=True)
        ax.set_xlabel(r"Muon Angle w.r.t. Beam [rad]")
        ax.set_ylabel(r"Count [Scaled to 1.2e19 POT] / 0.1 rad")
        plt.legend()
        plt.yscale('log')
        plt.xscale('log')
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        ## PLOT: truth-level muon mometum
        fig, ax = plt.subplots()
        hist_data = []
        hist_weights = []
        for i in range(len(list_of_sig_bkg_dfs)):
            data = np.array(list_of_sig_bkg_dfs[i]['muon_momentum'])/1000.
            counts, bins = np.histogram(data, bins=np.linspace(0,20,81) )
            hist_data.append(bins[:-1])
            hist_weights.append(counts*scale_factor)
        hist_colors = cmr.take_cmap_colors('nipy_spectral', len(hist_labels))
        ax.hist(hist_data, bins=bins, weights = hist_weights, histtype='stepfilled', alpha=0.9, color=hist_colors, label=hist_labels, stacked=True)
        ax.set_xlabel(r"Muon Momentum [GeV]")
        ax.set_ylabel(r"Count [Scaled to 1.2e19 POT] / 0.25 GeV")
        plt.legend()
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

         ## PLOT: primary shower available energy
        fig, ax = plt.subplots()
        hist_data = []
        hist_weights = []
        for i in range(len(list_of_sig_bkg_dfs)):
            data = np.array(list_of_sig_bkg_dfs[i]['primary_shower_available_energy'])
            data = np.concatenate(data)
            counts, bins = np.histogram(data, bins=np.linspace(0,2000,81) )
            hist_data.append(bins[:-1])
            hist_weights.append(counts*21.504)
        hist_colors = cmr.take_cmap_colors('nipy_spectral', len(hist_labels))
        ax.hist(hist_data, bins=bins, weights = hist_weights, histtype='stepfilled', alpha=0.9, color=hist_colors, label=hist_labels, stacked=True)
        ax.set_xlabel(r"Primary Shower Available Energy [MeV]")
        ax.set_ylabel(r"Count [Scaled to 1.2e19 POT] / 25 MeV")
        plt.yscale('log')
        #plt.xscale('log')
        plt.legend()
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

         ## PLOT: pi0 child available energy
        fig, ax = plt.subplots()
        hist_data = []
        hist_weights = []
        for i in range(len(list_of_sig_bkg_dfs)):
            data = list(itertools.chain.from_iterable(list_of_sig_bkg_dfs[i]['pi0_child_available_energy']))
            if len(data) == 0:
                data = np.array([])
            else:
                data = np.concatenate(data)
            counts, bins = np.histogram(data, bins=np.linspace(0,2000,81) )
            hist_data.append(bins[:-1])
            hist_weights.append(counts*21.504)
        hist_colors = cmr.take_cmap_colors('nipy_spectral', len(hist_labels))
        ax.hist(hist_data, bins=bins, weights = hist_weights, histtype='stepfilled', alpha=0.9, color=hist_colors, label=hist_labels, stacked=True)
        ax.set_xlabel(r"$\pi^0$ Child Available Energy [MeV]")
        ax.set_ylabel(r"Count [Scaled to 1.2e19 POT] / 25 MeV")
        #plt.yscale('log')
        #plt.xscale('log')
        plt.legend()
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        ## PLOT: primary shower deposited energy
        fig, ax = plt.subplots()
        hist_data = []
        hist_weights = []
        for i in range(len(list_of_sig_bkg_dfs)):
            data = np.array(list_of_sig_bkg_dfs[i]['primary_shower_edep'])
            data = np.concatenate(data)
            counts, bins = np.histogram(data, bins=np.linspace(0,2000,81) )
            hist_data.append(bins[:-1])
            hist_weights.append(counts*21.504)
        hist_colors = cmr.take_cmap_colors('nipy_spectral', len(hist_labels))
        ax.hist(hist_data, bins=bins, weights = hist_weights, histtype='stepfilled', alpha=0.9, color=hist_colors, label=hist_labels, stacked=True)
        ax.set_xlabel(r"Primary Shower Deposited Energy [MeV]")
        ax.set_ylabel(r"Count [Scaled to 1.2e19 POT] / 25 MeV")
        plt.yscale('log')
        #plt.xscale('log')
        plt.legend()
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        ## PLOT: pi0 child deposited energy
        fig, ax = plt.subplots()
        hist_data = []
        hist_weights = []
        for i in range(len(list_of_sig_bkg_dfs)):
            data = list(itertools.chain.from_iterable(list_of_sig_bkg_dfs[i]['pi0_child_edep']))
            if len(data) == 0:
                data = np.array([])
            else:
                data = np.concatenate(data)
            counts, bins = np.histogram(data, bins=np.linspace(0,2000,81) )
            hist_data.append(bins[:-1])
            hist_weights.append(counts*21.504)
        hist_colors = cmr.take_cmap_colors('nipy_spectral', len(hist_labels))
        ax.hist(hist_data, bins=bins, weights = hist_weights, histtype='stepfilled', alpha=0.9, color=hist_colors, label=hist_labels, stacked=True)
        ax.set_xlabel(r"$\pi^0$ Child Deposited Energy [MeV]")
        ax.set_ylabel(r"Count [Scaled to 1.2e19 POT] / 25 MeV")
        #plt.yscale('log')
        #plt.xscale('log')
        plt.legend()
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)


    return

def main(scale_factor, sig_bkg_json_file):

    sig_bkg_file = open(sig_bkg_json_file)
    sig_bkg_dict=json.load(sig_bkg_file)
    sig_bkg_df = pd.DataFrame(sig_bkg_dict)
    sig_bkg_df = sig_bkg_df.transpose()
    #print(sig_bkg_df)

    plot_sig_and_bkg_CC_1pi0(sig_bkg_df, scale_factor)

    return


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-sf', '--scale_factor', default=1, required=True, type=float, \
                        help='''Scale factor related to input dictionary''')
    parser.add_argument('-d', '--sig_bkg_json_file', default=None, required=True, type=str, \
                        help='''string corresponding to the path of the signal and background info JSON file''')
    args = parser.parse_args()
    main(**vars(args))