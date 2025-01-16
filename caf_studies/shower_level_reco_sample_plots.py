################################################################################
##                                                                            ##
##    CONTAINS: Script to create plots of reconstructed shower variables vs.  ##
##              truth match shower variables.                                 ##
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
    num_showers = len(df)

    print('\n----------------- File content -----------------')
    print('Reco Sample File:', reco_sample_file)
    print('Number of reconstructed showers:', num_showers)
    print('------------------------------------------------\n')

    # Set up output file
    output_pdf_name = reco_sample_file.split('.')[0]+'_validations.pdf'

    # Make plots in output file
    with PdfPages(output_pdf_name) as output:

        # Plot reco vs. true PDG ID
        fig, ax = plt.subplots(1,2, figsize=(10, 8))
        ax[0].pie(df['reco_pdg'].value_counts(), labels = df['reco_pdg'].value_counts().index, colors=['#ff7f0e', '#1f77b4'], autopct='%1.1f%%')
        ax[0].set_title('Reco PDG ID')

        ax[1].pie(df['true_pdg'].value_counts(), labels = df['true_pdg'].value_counts().index, autopct='%1.1f%%')
        ax[1].set_title('True PDG ID')
        fig.suptitle('PDG ID Distribution for Reconstructed and Matched True Showers')
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        # Plot reco PDG ID for showers with true PDG = 22 and true PDG ID for reco PDG = 22
        df_true_gamma = df[df['true_pdg'] == 22]
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.pie(df_true_gamma['reco_pdg'].value_counts(), labels = df_true_gamma['reco_pdg'].value_counts().index, autopct='%1.1f%%')
        ax.set_title('Reco PDG ID for Showers Matched to True Photons')
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        df_reco_gamma = df[df['reco_pdg'] == 22]
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.pie(df_reco_gamma['true_pdg'].value_counts(), labels = df_reco_gamma['true_pdg'].value_counts().index, autopct='%1.1f%%')
        ax.set_title('True PDG ID for Showers Reconstructed as Photons')
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)


        # Look at true vs. reco energy and true - reco energy [LOG SCALE]
        fig, ax = plt.subplots(1,2, figsize=(14, 6))
        ax[0].hist(df['true_energy']*1000, weights=np.full(num_showers, 1/num_showers), bins=np.linspace(0,2000,81), alpha=0.5, label='True Energy', color='b')
        ax[0].hist(df['reco_energy']*1000, weights=np.full(num_showers, 1/num_showers), bins=np.linspace(0,2000,81), alpha=0.5, label='Reco Energy', color='r')
        ax[0].set_xlabel('Energy [MeV]')
        ax[0].set_ylabel('Fraction of Showers / 25 MeV')
        ax[0].legend()
        ax[0].set_xlim(0, 2000)
        ax[0].set_yscale('log')

        ax[1].hist((df['true_energy']-df['reco_energy'])*1000, weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-2000,2000,161), alpha=0.5, color='g')
        ax[1].set_xlabel('(True - Reco) Energy [MeV]')
        ax[1].set_ylabel('Fraction of Showers / 25 MeV')
        ax[1].set_xlim(-2000, 2000)
        ax[1].set_yscale('log')
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        # Look at true vs. reco energy and true - reco energy [LINEAR SCALE]
        fig, ax = plt.subplots(1,2, figsize=(14, 6))
        ax[0].hist(df['true_energy']*1000, weights=np.full(num_showers, 1/num_showers), bins=np.linspace(0,2000,81), alpha=0.5, label='True Energy', color='b')
        ax[0].hist(df['reco_energy']*1000, weights=np.full(num_showers, 1/num_showers), bins=np.linspace(0,2000,81), alpha=0.5, label='Reco Energy', color='r')
        ax[0].set_xlabel('Energy [MeV]')
        ax[0].set_ylabel('Fraction of Showers / 25 MeV')
        ax[0].legend()
        ax[0].set_xlim(0, 2000)

        ax[1].hist((df['true_energy']-df['reco_energy'])*1000, weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-2000,2000,161), alpha=0.5, color='g')
        ax[1].set_xlabel('(True - Reco) Energy [MeV]')
        ax[1].set_ylabel('Fraction of Showers / 25 MeV')
        ax[1].set_xlim(-2000, 2000)
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        # Look at true vs. reco shower start to vertex distance (and true-reco distance) [LOG SCALE]
        fig, ax = plt.subplots(1,2, figsize=(14, 6))
        ax[0].hist(df['true_shower_to_vtx_dist'], weights=np.full(num_showers, 1/num_showers), bins=np.linspace(0,200,101), alpha=0.5, label='True', color='b')
        ax[0].hist(df['reco_shower_to_vtx_dist'], weights=np.full(num_showers, 1/num_showers), bins=np.linspace(0,200,101), alpha=0.5, label='Reco', color='r')
        ax[0].set_xlabel('Vertex to Shower Start Distance [cm]')
        ax[0].set_ylabel('Fraction of Showers / 2 cm')
        ax[0].legend()
        ax[0].set_xlim(0, 200)
        ax[0].set_yscale('log')

        ax[1].hist((df['true_shower_to_vtx_dist']-df['reco_shower_to_vtx_dist']), weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-200,200,201), alpha=0.5, color='g')
        ax[1].set_xlabel('(True - Reco) Vertex to Shower Start Distance [cm]')
        ax[1].set_ylabel('Fraction of Showers / 2 cm')
        ax[1].set_xlim(-200, 200)
        ax[1].set_yscale('log')
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        # Look at true vs. reco shower start to vertex distance (and true-reco distance) [LINEAR SCALE]
        fig, ax = plt.subplots(1,2, figsize=(14, 6))
        ax[0].hist(df['true_shower_to_vtx_dist'], weights=np.full(num_showers, 1/num_showers), bins=np.linspace(0,200,101), alpha=0.5, label='True', color='b')
        ax[0].hist(df['reco_shower_to_vtx_dist'], weights=np.full(num_showers, 1/num_showers), bins=np.linspace(0,200,101), alpha=0.5, label='Reco', color='r')
        ax[0].set_xlabel('Vertex to Shower Start Distance [cm]')
        ax[0].set_ylabel('Fraction of Showers / 2 cm')
        ax[0].legend()
        ax[0].set_xlim(0, 200)

        ax[1].hist((df['true_shower_to_vtx_dist']-df['reco_shower_to_vtx_dist']), weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-200,200,201), alpha=0.5, color='g')
        ax[1].set_xlabel('(True - Reco) Vertex to Shower Start Distance [cm]')
        ax[1].set_ylabel('Fraction of Showers / 2 cm')
        ax[1].set_xlim(-200, 200)
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        # Look at true vs. reco shower angle with respect to the beam (and true-reco angle) [LOG SCALE]
        fig, ax = plt.subplots(1,2, figsize=(14, 6))
        ax[0].hist(df['true_angle'], weights=np.full(num_showers, 1/num_showers), bins=np.linspace(0, 3.15,64), alpha=0.5, label='True', color='b')
        ax[0].hist(df['reco_angle'], weights=np.full(num_showers, 1/num_showers), bins=np.linspace(0, 3.15,64), alpha=0.5, label='Reco', color='r')
        ax[0].set_xlabel('Shower Angle with Respect to Beam [rad]')
        ax[0].set_ylabel('Fraction of Showers / 0.05 rad')
        ax[0].legend()
        ax[0].set_xlim(0, 3.15)
        ax[0].set_yscale('log')

        ax[1].hist((df['true_angle']-df['reco_angle']), weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-3.15, 3.15,127), alpha=0.5, color='g')
        ax[1].set_xlabel('(True - Reco) Shower Angle with Respect to Beam [rad]')
        ax[1].set_ylabel('Fraction of Showers / 0.05 rad')
        ax[1].set_xlim(-3.15, 3.15)
        ax[1].set_yscale('log')
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        # Look at true vs. reco shower angle with respect to the beam (and true-reco angle) [LINEAR SCALE]
        fig, ax = plt.subplots(1,2, figsize=(14, 6))
        ax[0].hist(df['true_angle'], weights=np.full(num_showers, 1/num_showers), bins=np.linspace(0, 3.15,64), alpha=0.5, label='True', color='b')
        ax[0].hist(df['reco_angle'], weights=np.full(num_showers, 1/num_showers), bins=np.linspace(0, 3.15,64), alpha=0.5, label='Reco', color='r')
        ax[0].set_xlabel('Shower Angle with Respect to Beam [rad] ')
        ax[0].set_ylabel('Fraction of Showers / 0.05 rad')
        ax[0].legend()
        ax[0].set_xlim(0, 3.15)

        ax[1].hist((df['true_angle']-df['reco_angle']), weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-3.15, 3.15,127), alpha=0.5, color='g')
        ax[1].set_xlabel('(True - Reco) Shower Angle with Respect to Beam [rad]')
        ax[1].set_ylabel('Fraction of Showers / 0.05 rad')
        ax[1].set_xlim(-3.15, 3.15)
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        # Look at true vs. reco vertex location (and true - reco location) [LINEAR SCALE] for x,y,z
        fig, ax = plt.subplots(2,3, figsize=(16, 8))

        # Vertex X
        ax[0,0].hist(df['true_ixn_vtx_x_pos'], weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-100, 100,101), alpha=0.5, label='True', color='b')
        ax[0,0].hist(df['reco_ixn_vtx_x_pos'], weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-100, 100,101), alpha=0.5, label='Reco', color='r')
        ax[0,0].set_xlabel('Vertex Location - X [cm]')
        ax[0,0].set_ylabel('Fraction of Showers /  2 cm')
        ax[0,0].set_xlim(-100, 100)
        ax[0,0].set_ylim(0, 0.04)

        ax[1,0].hist((df['true_ixn_vtx_x_pos']-df['reco_ixn_vtx_x_pos']), weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-150, 150,151), alpha=0.5, color='g')
        ax[1,0].set_xlabel('(True - Reco) Vertex Location - X [cm]')
        ax[1,0].set_ylabel('Fraction of Showers / 2 cm')
        ax[1,0].set_xlim(-150, 150)
        ax[1,0].set_ylim(0, 0.19)

        # Vertex Y
        ax[0,1].hist(df['true_ixn_vtx_y_pos'], weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-100, 100,101), alpha=0.5, label='True', color='b')
        ax[0,1].hist(df['reco_ixn_vtx_y_pos'], weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-100, 100,101), alpha=0.5, label='Reco', color='r')
        ax[0,1].set_xlabel('Vertex Location - Y [cm]')
        ax[0,1].set_ylabel('Fraction of Showers /  2 cm')
        ax[0,1].legend()
        ax[0,1].set_xlim(-100, 100)
        ax[0,1].set_ylim(0, 0.04)
        ax[0,1].yaxis.set_visible(False) 

        ax[1,1].hist((df['true_ixn_vtx_y_pos']-df['reco_ixn_vtx_y_pos']), weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-150, 150,151), alpha=0.5, color='g')
        ax[1,1].set_xlabel('(True - Reco) Vertex Location - Y [cm]')
        ax[1,1].set_ylabel('Fraction of Showers / 2 cm')
        ax[1,1].set_xlim(-150, 150)
        ax[1,1].set_ylim(0, 0.19)
        ax[1,1].yaxis.set_visible(False) 

        # Vertex Z
        ax[0,2].hist(df['true_ixn_vtx_z_pos'], weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-100, 100,101), alpha=0.5, label='True', color='b')
        ax[0,2].hist(df['reco_ixn_vtx_z_pos'], weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-100, 100,101), alpha=0.5, label='Reco', color='r')
        ax[0,2].set_xlabel('Vertex Location - Z [cm]')
        ax[0,2].set_ylabel('Fraction of Showers /  2 cm')
        ax[0,2].set_xlim(-100, 100)
        ax[0,2].set_ylim(0, 0.04)
        ax[0,2].yaxis.set_visible(False) 

        ax[1,2].hist((df['true_ixn_vtx_z_pos']-df['reco_ixn_vtx_z_pos']), weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-150, 150,151), alpha=0.5, color='g')
        ax[1,2].set_xlabel('(True - Reco) Vertex Location - Z [cm]')
        ax[1,2].set_ylabel('Fraction of Showers / 2 cm')
        ax[1,2].set_xlim(-150, 150)
        ax[1,2].set_ylim(0, 0.19)
        ax[1,2].yaxis.set_visible(False) 
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)


        # Look at true vs. reco Shower Start location (and true - reco location) [LINEAR SCALE] for x,y,z
        fig, ax = plt.subplots(2,3, figsize=(16, 8))

        # Shower Start X
        ax[0,0].hist(df['true_shower_start_x'], weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-100, 100,101), alpha=0.5, label='True', color='b')
        ax[0,0].hist(df['reco_shower_start_x'], weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-100, 100,101), alpha=0.5, label='Reco', color='r')
        ax[0,0].set_xlabel('Shower Start Location - X [cm]')
        ax[0,0].set_ylabel('Fraction of Showers /  2 cm')
        ax[0,0].set_xlim(-100, 100)
        ax[0,0].set_ylim(0, 0.05)

        ax[1,0].hist((df['true_shower_start_x']-df['reco_shower_start_x']), weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-150, 150,151), alpha=0.5, color='g')
        ax[1,0].set_xlabel('(True - Reco) Shower Start Location - X [cm]')
        ax[1,0].set_ylabel('Fraction of Showers / 2 cm')
        ax[1,0].set_xlim(-150, 150)
        ax[1,0].set_ylim(0, 0.175)

        # Shower Start Y
        ax[0,1].hist(df['true_shower_start_y'], weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-100, 100,101), alpha=0.5, label='True', color='b')
        ax[0,1].hist(df['reco_shower_start_y'], weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-100, 100,101), alpha=0.5, label='Reco', color='r')
        ax[0,1].set_xlabel('Shower Start Location - Y [cm]')
        ax[0,1].set_ylabel('Fraction of Showers /  2 cm')
        ax[0,1].legend()
        ax[0,1].set_xlim(-100, 100)
        ax[0,1].set_ylim(0, 0.05)
        ax[0,1].yaxis.set_visible(False) 

        ax[1,1].hist((df['true_shower_start_y']-df['reco_shower_start_y']), weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-150, 150,151), alpha=0.5, color='g')
        ax[1,1].set_xlabel('(True - Reco) Shower Start Location - Y [cm]')
        ax[1,1].set_ylabel('Fraction of Showers / 2 cm')
        ax[1,1].set_xlim(-150, 150)
        ax[1,1].set_ylim(0, 0.175)
        ax[1,1].yaxis.set_visible(False) 

        # Shower Start Z
        ax[0,2].hist(df['true_shower_start_z'], weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-100, 100,101), alpha=0.5, label='True', color='b')
        ax[0,2].hist(df['reco_shower_start_z'], weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-100, 100,101), alpha=0.5, label='Reco', color='r')
        ax[0,2].set_xlabel('Shower Start Location - Z [cm]')
        ax[0,2].set_ylabel('Fraction of Showers /  2 cm')
        ax[0,2].set_xlim(-100, 100)
        ax[0,2].set_ylim(0, 0.05)
        ax[0,2].yaxis.set_visible(False) 

        ax[1,2].hist((df['true_shower_start_z']-df['reco_shower_start_z']), weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-150, 150,151), alpha=0.5, color='g')
        ax[1,2].set_xlabel('(True - Reco) Shower Start Location - Z [cm]')
        ax[1,2].set_ylabel('Fraction of Showers / 2 cm')
        ax[1,2].set_xlim(-150, 150)
        ax[1,2].set_ylim(0, 0.175)
        ax[1,2].yaxis.set_visible(False) 
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        # Look at true vs. reco momentum components (and true - reco momentum) [LINEAR SCALE] for x,y,z
        fig, ax = plt.subplots(2,3, figsize=(16, 8))

        # Shower Momentum X
        ax[0,0].hist(df['true_p_x']*1000, weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-100, 100,101), alpha=0.5, label='True', color='b')
        ax[0,0].hist(df['reco_p_x']*1000, weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-100, 100,101), alpha=0.5, label='Reco', color='r')
        ax[0,0].set_xlabel('Shower Momentum - X [MeV/c]')
        ax[0,0].set_ylabel('Fraction of Showers /  2 MeV/c')
        ax[0,0].set_xlim(-100, 100)
        #ax[0,0].set_ylim(0, 0.04)
        #ax[0,0].set_yscale('log')

        # Look at true vs. reco energy
        ax[1,0].hist((df['true_p_x']-df['reco_p_x'])*1000, weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-150, 150,151), alpha=0.5, color='g')
        ax[1,0].set_xlabel('(True - Reco) Shower Momentum - X [MeV/c]')
        ax[1,0].set_ylabel('Fraction of Showers / 2 MeV/c')
        ax[1,0].set_xlim(-150, 150)
        #ax[1,0].set_ylim(0, 0.19)
        #ax[1,0].set_yscale('log')

        # Shower Momentum Y
        ax[0,1].hist(df['true_p_y'], weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-100, 100,101), alpha=0.5, label='True', color='b')
        ax[0,1].hist(df['reco_p_y'], weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-100, 100,101), alpha=0.5, label='Reco', color='r')
        ax[0,1].set_xlabel('Shower Momentum - Y [MeV/c]')
        ax[0,1].set_ylabel('Fraction of Showers /  2 MeV/c')
        ax[0,1].legend()
        ax[0,1].set_xlim(-100, 100)
        #ax[0,1].set_ylim(0, 0.04)
        #ax[0,1].yaxis.set_visible(False) 
        #ax[0,1].set_yscale('log')

        ax[1,1].hist((df['true_p_y']-df['reco_p_y'])*1000, weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-150, 150,151), alpha=0.5, color='g')
        ax[1,1].set_xlabel('(True - Reco) Shower Momentum - Y [MeV/c]')
        ax[1,1].set_ylabel('Fraction of Showers / 2 MeV/c')
        ax[1,1].set_xlim(-150, 150)
        #ax[1,1].set_ylim(0, 0.19)
        #ax[1,1].yaxis.set_visible(False) 
        #ax[1,1].set_yscale('log')

        # Shower Momentum Z
        ax[0,2].hist(df['true_p_z'], weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-100, 100,101), alpha=0.5, label='True', color='b')
        ax[0,2].hist(df['reco_p_z'], weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-100, 100,101), alpha=0.5, label='Reco', color='r')
        ax[0,2].set_xlabel('Shower Momentum - Z [MeV/c]')
        ax[0,2].set_ylabel('Fraction of Showers /  2 MeV/c')
        ax[0,2].set_xlim(-100, 100)
        #ax[0,2].set_ylim(0, 0.04)
        #ax[0,2].yaxis.set_visible(False) 
        #ax[0,2].set_yscale('log')

        ax[1,2].hist((df['true_p_z']-df['reco_p_z'])*1000, weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-150, 150,151), alpha=0.5, color='g')
        ax[1,2].set_xlabel('(True - Reco) Shower Momentum - Z [MeV/c]')
        ax[1,2].set_ylabel('Fraction of Showers / 2 MeV/c')
        ax[1,2].set_xlim(-150, 150)
        #ax[1,2].set_ylim(0, 0.19)
        #ax[1,2].yaxis.set_visible(False) 
        #ax[1,2].set_yscale('log')
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        # Look at true vs. reco momentum magnitude (and true - reco momentum) [LOG SCALE]
        fig, ax = plt.subplots(1,2, figsize=(14, 6))
        ax[0].hist(df['true_p_mag']*1000, weights=np.full(num_showers, 1/num_showers), bins=np.linspace(0,2000,81), alpha=0.5, label='True', color='b')
        ax[0].hist(df['reco_p_mag']*1000, weights=np.full(num_showers, 1/num_showers), bins=np.linspace(0,2000,81), alpha=0.5, label='Reco', color='r')
        ax[0].set_xlabel('Shower Momentum [MeV/c]')
        ax[0].set_ylabel('Fraction of Showers')
        ax[0].legend()
        ax[0].set_xlim(0, 2000)
        ax[0].set_yscale('log')

        ax[1].hist((df['true_p_mag']-df['reco_p_mag'])*1000, weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-2000,2000,161), alpha=0.5, color='g')
        ax[1].set_xlabel('(True - Reco) Shower Momentum [MeV/c]')
        ax[1].set_ylabel('Fraction of Showers')
        ax[1].set_xlim(-2000, 2000)
        ax[1].set_yscale('log')
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        # Look at true vs. reco momentum magnitude (and true - reco momentum) [LINEAR SCALE]
        fig, ax = plt.subplots(1,2, figsize=(14, 6))
        ax[0].hist(df['true_p_mag']*1000, weights=np.full(num_showers, 1/num_showers), bins=np.linspace(0,2000,81), alpha=0.5, label='True', color='b')
        ax[0].hist(df['reco_p_mag']*1000, weights=np.full(num_showers, 1/num_showers), bins=np.linspace(0,2000,81), alpha=0.5, label='Reco', color='r')
        ax[0].set_xlabel('Shower Momentum [MeV/c]')
        ax[0].set_ylabel('Fraction of Showers')
        ax[0].legend()
        ax[0].set_xlim(0, 2000)

        # Look at true vs. reco energy
        ax[1].hist((df['true_p_mag']-df['reco_p_mag'])*1000, weights=np.full(num_showers, 1/num_showers), bins=np.linspace(-2000,2000,161), alpha=0.5, color='g')
        ax[1].set_xlabel('(True - Reco) Shower Momentum [MeV/c]')
        ax[1].set_ylabel('Fraction of Showers')
        ax[1].set_xlim(-2000, 2000)   
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