################################################################################
##                                                                            ##
##    CONTAINS: Script to create plots of reconstructed ixn variables vs.     ##
##              truth match shower variables with root file from sig bkg dict.##
##                                                                            ##
################################################################################

# LOAD Packages 
import uproot
import pandas as pd
import numpy as np
import awkward as ak
from matplotlib.axes import Axes
import matplotlib as mpl
import awkward_pandas
from mpl_toolkits.mplot3d import Axes3D
#matplotlib.rcParams.update(matplotlib.rcParamsDefault)
#matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
from matplotlib import cm, colors
import matplotlib.patches as mpatches
import h5py
import argparse
from matplotlib.backends.backend_pdf import PdfPages

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
rasterize_plots()

def main(reco_sample_file):

    # Load reconstructed shower sample
    ixns = uproot.open(reco_sample_file)
    ak_array = ixns["RecoBenchmarkTree"].arrays(library="ak")
    df = ak.to_dataframe(ak_array)
    df = df.reset_index(drop=True)
    num_reco_ixns = len(df)
    num_truth_ixns = len(df['truth_dict_key'].unique())
    df_unique_match = df.groupby('truth_dict_key').filter(lambda x: len(x) == 1)
    df_multi_match = df.groupby('truth_dict_key').filter(lambda x: len(x) > 1)
    num_reco_ixns_unique_match = len(df_unique_match)
    num_reco_ixns_multi_match = len(df_multi_match)
    num_truth_ixns_multi_match = len(df_multi_match['truth_dict_key'].unique())

    print('\n----------------- File content -----------------')
    print('Reco Sample File:', reco_sample_file)
    print('Number of reco ixns (TOTAL):', num_reco_ixns)
    print('Number of truth ixns (TOTAL):', num_truth_ixns)
    print('Number of reco ixns w/ unique truth match:', num_reco_ixns_unique_match)
    print('Number of reco ixns 1:N Truth:Reco:', num_reco_ixns_multi_match)
    print('Number of truth ixns 1:N Truth:Reco:', num_truth_ixns_multi_match)
    print('------------------------------------------------\n')

    # Set up output file
    output_pdf_name = reco_sample_file.split('.')[0]+'_validations.pdf'

    # Make plots in output file
    with PdfPages(output_pdf_name) as output:

        # Plot 1:1 and 1:N ixn overlap distributions
        ovlp_bins = np.linspace(0.,1.,101)

        # Linear scale
        fig, ax = plt.subplots()
        plt.hist(df_unique_match['overlap'], histtype='stepfilled',  bins=ovlp_bins, label='1:1 Truth:Reco, Full Sample', color='royalblue', alpha=0.5)
        plt.hist(df_multi_match['overlap'], histtype='stepfilled', bins=ovlp_bins,  label='1:N Truth:Reco, Full Sample', color='darkorange', alpha=0.5)
        ax.set_xlabel(r"Interaction Overlap Value [0,1]")
        ax.set_ylabel(r"Count [Scaled to 1.2e19 POT] / Reco Interaction")
        plt.legend()
        plt.xlim(0.,1.)
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        # Log scale
        fig, ax = plt.subplots()
        plt.hist(df_unique_match['overlap'], histtype='stepfilled',  bins=ovlp_bins, label='1:1 Truth:Reco, Full Sample', color='royalblue', alpha=0.5)
        plt.hist(df_multi_match['overlap'], histtype='stepfilled',  bins=ovlp_bins, label='1:N Truth:Reco, Full Sample', color='darkorange', alpha=0.5)
        ax.set_xlabel(r"Interaction Overlap Value [0,1]")
        ax.set_ylabel(r"Count [Scaled to 1.2e19 POT] / Reco Interaction")
        plt.legend()
        ax.set_yscale('log')
        plt.xlim(0.,1.)
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        df_cc1pi0 = df[df['true_ixn_pi0_mult'] == 1]
        df_cc1pi0_unique_match = df_cc1pi0.groupby('truth_dict_key').filter(lambda x: len(x) == 1)
        df_cc1pi0_multi_match = df_cc1pi0.groupby('truth_dict_key').filter(lambda x: len(x) > 1)
        print("Number of 1:1 truth:reco interactions (1pi0):", len(df_cc1pi0_unique_match))
        print("Number of 1:N truth:reco interactions (entries) (1pi0):", len(df_cc1pi0_multi_match))
        print("Number of 1:N truth:reco interactions (unique truth_dict_keys) (1pi0):", len(df_cc1pi0_multi_match['truth_dict_key'].unique()))
        df_cc1pi0_no_chpi = df_cc1pi0[df_cc1pi0['true_ixn_chpi_mult'] == 0]
        df_cc1pi0_no_chmeson = df_cc1pi0_no_chpi[df_cc1pi0_no_chpi['true_ixn_chkaon_mult'] == 0]
        print("Number of CC 1pi0 0 charged meson ixns:", len(df_cc1pi0_no_chmeson['truth_dict_key'].unique()))
        df_cc1pi0_no_chmeson_unique_match = df_cc1pi0_no_chmeson.groupby('truth_dict_key').filter(lambda x: len(x) == 1)
        df_cc1pi0_no_chmeson_multi_match = df_cc1pi0_no_chmeson.groupby('truth_dict_key').filter(lambda x: len(x) > 1)
        print("Number of 1:1 truth:reco interactions (1pi0, 0chme):", len(df_cc1pi0_no_chmeson_unique_match))
        print("Number of 1:N truth:reco interactions (entries) (1pi0, 0chme):", len(df_cc1pi0_no_chmeson_multi_match))
        print("Number of 1:N truth:reco interactions (unique truth_dict_keys) (1pi0, 0chme):", len(df_cc1pi0_no_chmeson_multi_match['truth_dict_key'].unique()))

        # Histogram -- ixn overlap (1:1, 1:N) -- 1pi0
        ovlp_bins = np.linspace(0.,1.,101)
        fig, ax = plt.subplots()
        plt.hist(df_cc1pi0_unique_match['overlap'], histtype='stepfilled',  bins=ovlp_bins, label=r'1:1 Truth:Reco, CC1$\pi^0$', color='royalblue', alpha=0.5)
        plt.hist(df_cc1pi0_multi_match['overlap'], histtype='stepfilled', bins=ovlp_bins,  label=r'1:N Truth:Reco, CC1$\pi^0$', color='darkorange', alpha=0.5)
        ax.set_xlabel(r"Interaction Overlap Value [0,1]")
        ax.set_ylabel(r"Count [Scaled to 1.2e19 POT] / Reco Interaction")
        plt.legend()
        plt.xlim(0.,1.)
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        fig, ax = plt.subplots()
        plt.hist(df_cc1pi0_unique_match['overlap'], histtype='stepfilled',  bins=ovlp_bins, label=r'1:1 Truth:Reco, CC1$\pi^0$', color='royalblue', alpha=0.5)
        plt.hist(df_cc1pi0_multi_match['overlap'], histtype='stepfilled',  bins=ovlp_bins, label=r'1:N Truth:Reco, CC1$\pi^0$', color='darkorange', alpha=0.5)
        ax.set_xlabel(r"Interaction Overlap Value [0,1]")
        ax.set_ylabel(r"Count [Scaled to 1.2e19 POT] / Reco Interaction")
        plt.legend()
        ax.set_yscale('log')
        plt.xlim(0.,1.)
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        # Look at unique match distributions only -- CC 1pi0

        # Reco vs. true 
        # histogram -- muon mult cc1pi0
        muon_mult_bins = np.linspace(0,10,11)
        fig, ax = plt.subplots()
        plt.hist(df_cc1pi0_unique_match['true_ixn_muon_mult'], histtype='stepfilled',  bins=muon_mult_bins, label=r'Truth', color='blueviolet', alpha=0.5)
        plt.hist(df_cc1pi0_unique_match['reco_ixn_muon_mult'], histtype='stepfilled', bins=muon_mult_bins, label=r'Reco', color='orangered', alpha=0.5)
        ax.set_xlabel(r"Muon Multiplicity")
        ax.set_ylabel(r"Count [Scaled to 1.2e19 POT] / Interaction")
        plt.legend()
        plt.xlim(0,8)
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        # histogram -- muon mult difference
        muon_mult_diff_bins = np.linspace(-10,5,16)
        fig, ax = plt.subplots()
        plt.hist((df_cc1pi0_unique_match['true_ixn_muon_mult']-df_cc1pi0_unique_match['reco_ixn_muon_mult']), histtype='stepfilled', bins= muon_mult_diff_bins, label=r'Truth - Reco', color='forestgreen', alpha=0.5)
        ax.set_xlabel(r"Muon Multiplicity Difference (Truth-Reco)")
        ax.set_ylabel(r"Count [Scaled to 1.2e19 POT] / Interaction")
        #plt.legend()
        plt.xlim(-6,3)
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        # histogram -- charged pion mult cc1pi0
        mult_bins = np.linspace(0,20,21)
        fig, ax = plt.subplots()
        plt.hist(df_cc1pi0_unique_match['true_ixn_chpi_mult'], histtype='stepfilled',  bins=mult_bins, label=r'Truth', color='blueviolet', alpha=0.5)
        plt.hist(df_cc1pi0_unique_match['reco_ixn_chpi_mult'], histtype='stepfilled', bins=mult_bins, label=r'Reco', color='orangered', alpha=0.5)
        ax.set_xlabel(r"Charged Pion Multiplicity")
        ax.set_ylabel(r"Count [Scaled to 1.2e19 POT] / Interaction")
        plt.legend()
        plt.xlim(0,10)
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        # histogram -- charged pion difference cc1pi0
        mult_diff_bins = np.linspace(-20,20,41)
        fig, ax = plt.subplots()
        plt.hist((df_cc1pi0_unique_match['true_ixn_chpi_mult']-df_cc1pi0_unique_match['reco_ixn_chpi_mult']), histtype='stepfilled', bins= mult_diff_bins, label=r'Truth - Reco', color='forestgreen', alpha=0.5)
        ax.set_xlabel(r"Charged Pion  Multiplicity Difference (Truth-Reco)")
        ax.set_ylabel(r"Count [Scaled to 1.2e19 POT] / Interaction")
        #plt.legend()
        plt.xlim(-5,10)
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        # histogram -- charged kaon mult cc1pi0
        mult_bins = np.linspace(0,20,21)
        fig, ax = plt.subplots()
        plt.hist(df_cc1pi0_unique_match['true_ixn_chkaon_mult'], histtype='stepfilled',  bins=mult_bins, label=r'Truth', color='blueviolet', alpha=0.5)
        plt.hist(df_cc1pi0_unique_match['reco_ixn_chkaon_mult'], histtype='stepfilled', bins=mult_bins, label=r'Reco', color='orangered', alpha=0.5)
        ax.set_xlabel(r"Charged Kaon Multiplicity")
        ax.set_ylabel(r"Count [Scaled to 1.2e19 POT] / Interaction")
        plt.legend()
        plt.xlim(0,6)
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        # histogram -- charged kaon mult difference
        mult_diff_bins = np.linspace(-20,20,41)
        fig, ax = plt.subplots()
        plt.hist((df_cc1pi0_unique_match['true_ixn_chkaon_mult']-df_cc1pi0_unique_match['reco_ixn_chkaon_mult']), histtype='stepfilled', bins= mult_diff_bins, label=r'Truth - Reco', color='forestgreen', alpha=0.5)
        ax.set_xlabel(r"Charged Kaon  Multiplicity Difference (Truth-Reco)")
        ax.set_ylabel(r"Count [Scaled to 1.2e19 POT] / Interaction")
        #plt.legend()
        plt.xlim(-4,4)
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        # histogram -- proton mult cc1pi0
        mult_bins = np.linspace(0,20,21)
        fig, ax = plt.subplots()
        plt.hist(df_cc1pi0_unique_match['true_ixn_proton_mult'], histtype='stepfilled',  bins=mult_bins, label=r'Truth', color='blueviolet', alpha=0.5)
        plt.hist(df_cc1pi0_unique_match['reco_ixn_proton_mult'], histtype='stepfilled', bins=mult_bins, label=r'Reco', color='orangered', alpha=0.5)
        ax.set_xlabel(r"Proton Multiplicity")
        ax.set_ylabel(r"Count [Scaled to 1.2e19 POT] / Interaction")
        plt.legend()
        plt.xlim(0,20)
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        # histogram -- proton mult difference
        mult_diff_bins = np.linspace(-20,20,41)
        fig, ax = plt.subplots()
        plt.hist((df_cc1pi0_unique_match['true_ixn_proton_mult']-df_cc1pi0_unique_match['reco_ixn_proton_mult']), histtype='stepfilled', bins= mult_diff_bins, label=r'Truth - Reco', color='forestgreen', alpha=0.5)
        ax.set_xlabel(r"Proton Multiplicity Difference (Truth-Reco)")
        ax.set_ylabel(r"Count [Scaled to 1.2e19 POT] / Interaction")
        #plt.legend()
        plt.xlim(-5,20)
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        # histogram -- reco electron and gamma mult cc1pi0
        mult_bins = np.linspace(0,20,21)
        fig, ax = plt.subplots()
        plt.hist(df_cc1pi0_unique_match['reco_ixn_e_mult'], histtype='stepfilled',  bins=mult_bins, label=r'Reco Electrons', color='orange', alpha=0.5)
        plt.hist(df_cc1pi0_unique_match['reco_ixn_gamma_mult'], histtype='stepfilled', bins=mult_bins, label=r'Reco Photons', color='royalblue', alpha=0.5)
        ax.set_xlabel(r"Reco Interaction Shower Multiplicity")
        ax.set_ylabel(r"Count [Scaled to 1.2e19 POT] / Interaction")
        plt.legend()
        plt.xlim(0,7)
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        # histogram -- true primary electron and gamma mult difference
        #mult_diff_bins = np.linspace(-20,20,41)
        fig, ax = plt.subplots()
        plt.hist(df_cc1pi0_unique_match['true_ixn_e_mult'], histtype='stepfilled',  bins=mult_bins, label=r'Reco Electrons', color='orange', alpha=0.5)
        plt.hist(df_cc1pi0_unique_match['true_ixn_gamma_mult'], histtype='stepfilled', bins=mult_bins, label=r'Reco Photons', color='royalblue', alpha=0.5)
        ax.set_xlabel(r"Truth Interaction Primary Shower Multiplicity")
        ax.set_ylabel(r"Count [Scaled to 1.2e19 POT] / Interaction")
        #plt.legend()
        plt.xlim(0,5)
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        # Look at true vs. reco vertex location (and true - reco location) [LINEAR SCALE] for x,y,z
        fig, ax = plt.subplots(2,3, figsize=(16, 8))

        # Vertex X
        ax[0,0].hist(df_cc1pi0_unique_match['true_ixn_vtx_x_pos'], bins=np.linspace(-100, 100,101), alpha=0.5, label='True', color='b')
        ax[0,0].hist(df_cc1pi0_unique_match['reco_ixn_vtx_x_pos'], bins=np.linspace(-100, 100,101), alpha=0.5, label='Reco', color='r')
        ax[0,0].set_xlabel('Vertex Location - X [cm]')
        ax[0,0].set_ylabel('Interactions /  2 cm')
        ax[0,0].set_xlim(-75, 75)
        ax[0,0].set_ylim(0, 90)

        ax[1,0].hist((df_cc1pi0_unique_match['true_ixn_vtx_x_pos']-df_cc1pi0_unique_match['reco_ixn_vtx_x_pos']), bins=np.linspace(-150, 150,151), alpha=0.5, color='g')
        ax[1,0].set_xlabel('(True - Reco) Vertex Location - X [cm]')
        ax[1,0].set_ylabel('Interactions / 2 cm')
        ax[1,0].set_xlim(-75, 75)
        ax[1,0].set_ylim(0, 2000)

        # Vertex Y
        ax[0,1].hist(df_cc1pi0_unique_match['true_ixn_vtx_y_pos'], bins=np.linspace(-100, 100,101), alpha=0.5, label='True', color='b')
        ax[0,1].hist(df_cc1pi0_unique_match['reco_ixn_vtx_y_pos'], bins=np.linspace(-100, 100,101), alpha=0.5, label='Reco', color='r')
        ax[0,1].set_xlabel('Vertex Location - Y [cm]')
        ax[0,1].set_ylabel('Interactions /  2 cm')
        ax[0,1].legend()
        ax[0,1].set_xlim(-75, 75)
        ax[0,1].set_ylim(0, 90)
        ax[0,1].yaxis.set_visible(False) 

        ax[1,1].hist((df_cc1pi0_unique_match['true_ixn_vtx_y_pos']-df_cc1pi0_unique_match['reco_ixn_vtx_y_pos']), bins=np.linspace(-150, 150,151), alpha=0.5, color='g')
        ax[1,1].set_xlabel('(True - Reco) Vertex Location - Y [cm]')
        ax[1,1].set_ylabel('Interactions / 2 cm')
        ax[1,1].set_xlim(-75, 75)
        ax[1,1].set_ylim(0, 2000)
        ax[1,1].yaxis.set_visible(False) 

        # Vertex Z
        ax[0,2].hist(df_cc1pi0_unique_match['true_ixn_vtx_z_pos'],  bins=np.linspace(-100, 100,101), alpha=0.5, label='True', color='b')
        ax[0,2].hist(df_cc1pi0_unique_match['reco_ixn_vtx_z_pos'],  bins=np.linspace(-100, 100,101), alpha=0.5, label='Reco', color='r')
        ax[0,2].set_xlabel('Vertex Location - Z [cm]')
        ax[0,2].set_ylabel('Interactions /  2 cm')
        ax[0,2].set_xlim(-75, 75)
        ax[0,2].set_ylim(0, 90)
        ax[0,2].yaxis.set_visible(False) 

        ax[1,2].hist((df_cc1pi0_unique_match['true_ixn_vtx_z_pos']-df_cc1pi0_unique_match['reco_ixn_vtx_z_pos']),  bins=np.linspace(-150, 150,151), alpha=0.5, color='g')
        ax[1,2].set_xlabel('(True - Reco) Vertex Location - Z [cm]')
        ax[1,2].set_ylabel('Interactions / 2 cm')
        ax[1,2].set_xlim(-75, 75)
        ax[1,2].set_ylim(0, 2000)
        ax[1,2].yaxis.set_visible(False) 
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)


    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--reco_sample_file', default=None, type=str,\
                        help='''string corresponding to the path of the file containing a \
                                TTree with a sample of reconstructed shower characteristics''')
    #parser.add_argument('-c', '--overlap_cut', default=0.5, type=float,\
    #                    help='''float corresponding with truth/reco overlap cut value''')
    args = parser.parse_args()
main(**vars(args))