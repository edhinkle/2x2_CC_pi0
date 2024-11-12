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
        sample_title = 'Single '+r'$\pi^0$ $\bar{\nu}_\mu$/$\nu_\mu$ CC'
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
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        # PLOT: Number of pi0 decay products with associated trajectories
        fig, ax = plt.subplots(figsize=(6,4))
        data = np.array([np.count_nonzero(np.array(d[key]['child_pdg'])) for key in d.keys()])
        counts, bins =np.histogram(data, bins=np.linspace(0,5,6))
        ax.hist(bins[:-1], bins=bins, weights = counts*scale_factor, histtype='stepfilled', alpha=0.5, color='blue')
        ax.set_xlabel(r"Multiplicity of $\pi^0$ Decay Products which Deposit Energy in 2x2")
        ax.set_ylabel("Events / Decay Product Multiplicity")
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        # PLOT: Number pi0 decay products/mode multiplicity
        colors = ['lightskyblue', 'cyan', 'teal', 'royalblue','navy']
        colors = ['navy', 'royalblue', 'cyan', 'teal', 'lightskyblue']
        fig, ax = plt.subplots(figsize=(6,4))
        pi0_decay_prod_pdg_list=[sorted(d[key]['child_pdg']) for key in d.keys()]
        pi0_decay_prod_set = set(tuple(pdg) for pdg in pi0_decay_prod_pdg_list)
        pi0_decay_prod_count = [(pdg_set, pi0_decay_prod_pdg_list.count(list(pdg_set))) for pdg_set in pi0_decay_prod_set]
        pi0_decay_prod_fraction = [100*(i[1]/len(pi0_decay_prod_pdg_list)) for i in pi0_decay_prod_count]
        pi0_decay_prod_labels = [''.join(str(pdg_defs.pi0_decay_prod_pdg_dict[j]) for j in i[0]) for i in pi0_decay_prod_count]
        for i in range(len(pi0_decay_prod_labels)):
            if pi0_decay_prod_labels[i] == '':
                pi0_decay_prod_labels[i] = 'None'
        ax.pie(pi0_decay_prod_fraction, labels=pi0_decay_prod_labels, colors=colors, autopct='%1.1f%%', textprops={'color': 'black'})
        ax.set_title(r"Visible $\pi^0$ Decay Products in "+sample_title+" Events")
        plt.tight_layout()
        output.savefig(fig)
        plt.close(fig)

        gg = {k:v for k,v in d.items() if len(v['child_pdg']) == 2 and v['child_pdg'][0] == 22 and v['child_pdg'][1] == 22}
        print("Number of diphoton events:", len(gg))
        ee = {k:v for k,v in d.items() if len(v['child_pdg']) == 2 and abs(v['child_pdg'][0]) == 11 and abs(v['child_pdg'][1]) == 11}
        print("Number of dielectron events:", len(ee))
        g = {k:v for k,v in d.items() if np.count_nonzero(np.array(v['child_pdg'])) == 1 and max(v['child_pdg']) == 22}
        print("Number of single photon events:", len(g))
        nodec = {k:v for k,v in d.items() if np.count_nonzero(np.array(v['child_pdg'])) == 0}
        print("Number of no decay product events:", len(nodec))
        eeg = {k:v for k,v in d.items() if len(v['child_pdg']) == 3 and set(v['child_pdg']) == {11,-11,22}}
        print("Number of single photon dielectron events:", len(eeg))

        # PLOT: pi0 start location by decay product (XY)
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        gg_pi0_start = np.array([np.array(gg[key]['pi0_start_loc']) for key in gg.keys()])
        ax.scatter(gg_pi0_start[:,2], gg_pi0_start[:,0], gg_pi0_start[:,1], c=colors[4], marker='o', label='Diphoton')
        ee_pi0_start = np.array([np.array(ee[key]['pi0_start_loc']) for key in ee.keys()])
        ax.scatter(ee_pi0_start[:,2], ee_pi0_start[:,0], ee_pi0_start[:,1], c=colors[1], marker='x', label=r'$e^+ e^-$')
        g_pi0_start = np.array([np.array(g[key]['pi0_start_loc']) for key in g.keys()])
        ax.scatter(g_pi0_start[:,2], g_pi0_start[:,0], g_pi0_start[:,1], c=colors[2], marker='*', label='Single Photon')
        nodec_pi0_start = np.array([np.array(nodec[key]['pi0_start_loc']) for key in nodec.keys()])
        ax.scatter(nodec_pi0_start[:,2], nodec_pi0_start[:,0], nodec_pi0_start[:,1], c=colors[3], marker='^', label='No Decay Products')
        eeg_pi0_start = np.array([np.array(eeg[key]['pi0_start_loc']) for key in eeg.keys()])
        ax.scatter(eeg_pi0_start[:,2], eeg_pi0_start[:,0], eeg_pi0_start[:,1], c=colors[0], marker='v', label=r'Single Photon $e^+ e^-$')
        ax.set_xlabel('Beam Axis [cm]')
        ax.set_ylabel('Drift Axis [cm]')
        ax.set_zlabel('Vertical Axis [cm]')
        ax.set_title(r"$\pi^0$ Start Location by Decay Product")
        ax.legend()
        output.savefig(fig)
        plt.close(fig)
        #print("TPC Bounds 0: ", geo_defs.tpc_bounds(0))
        #counts, bins =np.histogram(data, bins=np.linspace(0,5,6))
        #ax.hist(bins[:-1], bins=bins, weights = counts*scale_factor, histtype='stepfilled', alpha=0.5, color='blue')

        # PLOT: pi0 start location by decay product (ZX, ZY, XY 2D projections)
        fig = plt.figure(figsize=(21, 6))
        # First subplot -- ZX
        ax1 = fig.add_subplot(131)
        ax1.scatter(gg_pi0_start[:,2], gg_pi0_start[:,0], c=colors[4], marker='o', label='Diphoton')
        ax1.scatter(ee_pi0_start[:,2], ee_pi0_start[:,0], c=colors[1], marker='x', label=r'$e^+ e^-$')
        ax1.scatter(g_pi0_start[:,2], g_pi0_start[:,0], c=colors[2], marker='*', label='Single Photon')
        ax1.scatter(nodec_pi0_start[:,2], nodec_pi0_start[:,0], c=colors[3], marker='^', label='No Decay Products')
        ax1.scatter(eeg_pi0_start[:,2], eeg_pi0_start[:,0], c=colors[0], marker='v', label=r'Single Photon $e^+ e^-$')
        ax1.set_xlabel('Beam Axis [cm]')
        ax1.set_ylabel('Drift Axis [cm]')
        ax1.set_title(r"$\pi^0$ Start Location by Decay Product")
        #ax1.legend()

        # Second subplot -- ZY
        ax2 = fig.add_subplot(132)
        ax2.scatter(gg_pi0_start[:,2], gg_pi0_start[:,1], c=colors[4], marker='o', label='Diphoton')
        ax2.scatter(ee_pi0_start[:,2], ee_pi0_start[:,1], c=colors[1], marker='x', label=r'$e^+ e^-$')
        ax2.scatter(g_pi0_start[:,2], g_pi0_start[:,1], c=colors[2], marker='*', label='Single Photon')
        ax2.scatter(nodec_pi0_start[:,2], nodec_pi0_start[:,1], c=colors[3], marker='^', label='No Decay Products')
        ax2.scatter(eeg_pi0_start[:,2], eeg_pi0_start[:,1], c=colors[0], marker='v', label=r'Single Photon $e^+ e^-$')
        ax2.set_xlabel('Beam Axis [cm]')
        ax2.set_ylabel('Vertical Axis [cm]')
        ax2.set_title(r"$\pi^0$ Start Location by Decay Product")
        #ax2.legend()

        # Third subplot -- XY
        ax3 = fig.add_subplot(133)
        ax3.scatter(gg_pi0_start[:,0], gg_pi0_start[:,1], c=colors[4], marker='o', label='Diphoton')
        ax3.scatter(ee_pi0_start[:,0], ee_pi0_start[:,1], c=colors[1], marker='x', label=r'$e^+ e^-$')
        ax3.scatter(g_pi0_start[:,0], g_pi0_start[:,1], c=colors[2], marker='*', label='Single Photon')
        ax3.scatter(nodec_pi0_start[:,0], nodec_pi0_start[:,1], c=colors[3], marker='^', label='No Decay Products')
        ax3.scatter(eeg_pi0_start[:,0], eeg_pi0_start[:,1], c=colors[0], marker='v', label=r'Single Photon $e^+ e^-$')
        ax3.set_xlabel('Drift Axis [cm]')
        ax3.set_ylabel('Vertical Axis [cm]')
        ax3.set_title(r"$\pi^0$ Start Location by Decay Product")
        #ax3.legend()

        output.savefig(fig)
        plt.close(fig)

        # PLOT: pi0 end location by decay product (XYZ)
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        gg_pi0_end = np.array([np.array(gg[key]['pi0_end_loc']) for key in gg.keys()])
        ax.scatter(gg_pi0_end[:,2], gg_pi0_end[:,0], gg_pi0_end[:,1], c=colors[4], marker='o', label='Diphoton')
        ee_pi0_end = np.array([np.array(ee[key]['pi0_end_loc']) for key in ee.keys()])
        ax.scatter(ee_pi0_end[:,2], ee_pi0_end[:,0], ee_pi0_end[:,1], c=colors[1], marker='x', label=r'$e^+ e^-$')
        g_pi0_end = np.array([np.array(g[key]['pi0_end_loc']) for key in g.keys()])
        ax.scatter(g_pi0_end[:,2], g_pi0_end[:,0], g_pi0_end[:,1], c=colors[2], marker='*', label='Single Photon')
        nodec_pi0_end = np.array([np.array(nodec[key]['pi0_end_loc']) for key in nodec.keys()])
        ax.scatter(nodec_pi0_end[:,2], nodec_pi0_end[:,0], nodec_pi0_end[:,1], c=colors[3], marker='^', label='No Decay Products')
        eeg_pi0_end = np.array([np.array(eeg[key]['pi0_end_loc']) for key in eeg.keys()])
        ax.scatter(eeg_pi0_end[:,2], eeg_pi0_end[:,0], eeg_pi0_end[:,1], c=colors[0], marker='v', label=r'Single Photon $e^+ e^-$')
        ax.set_xlabel('Beam Axis [cm]')
        ax.set_ylabel('Drift Axis [cm]')
        ax.set_zlabel('Vertical Axis [cm]')
        ax.set_title(r"$\pi^0$ End Location by Decay Product")
        ax.legend()
        output.savefig(fig)
        plt.close(fig)

        # PLOT: pi0 end location by decay product (ZX, ZY, XY 2D projections)
        fig = plt.figure(figsize=(21, 6))
        # First subplot -- ZX
        ax1 = fig.add_subplot(131)
        ax1.scatter(gg_pi0_end[:,2], gg_pi0_end[:,0], c=colors[4], marker='o', label='Diphoton')
        ax1.scatter(ee_pi0_end[:,2], ee_pi0_end[:,0], c=colors[1], marker='x', label=r'$e^+ e^-$')
        ax1.scatter(g_pi0_end[:,2], g_pi0_end[:,0], c=colors[2], marker='*', label='Single Photon')
        ax1.scatter(nodec_pi0_end[:,2], nodec_pi0_end[:,0], c=colors[3], marker='^', label='No Decay Products')
        ax1.scatter(eeg_pi0_end[:,2], eeg_pi0_end[:,0], c=colors[0], marker='v', label=r'Single Photon $e^+ e^-$')
        ax1.set_xlabel('Beam Axis [cm]')
        ax1.set_ylabel('Drift Axis [cm]')
        ax1.set_title(r"$\pi^0$ End Location by Decay Product")
        #ax1.legend()

        # Second subplot -- ZY
        ax2 = fig.add_subplot(132)
        ax2.scatter(gg_pi0_end[:,2], gg_pi0_end[:,1], c=colors[4], marker='o', label='Diphoton')
        ax2.scatter(ee_pi0_end[:,2], ee_pi0_end[:,1], c=colors[1], marker='x', label=r'$e^+ e^-$')
        ax2.scatter(g_pi0_end[:,2], g_pi0_end[:,1], c=colors[2], marker='*', label='Single Photon')
        ax2.scatter(nodec_pi0_end[:,2], nodec_pi0_end[:,1], c=colors[3], marker='^', label='No Decay Products')
        ax2.scatter(eeg_pi0_end[:,2], eeg_pi0_end[:,1], c=colors[0], marker='v', label=r'Single Photon $e^+ e^-$')
        ax2.set_xlabel('Beam Axis [cm]')
        ax2.set_ylabel('Vertical Axis [cm]')
        ax2.set_title(r"$\pi^0$ End Location by Decay Product")
        #ax2.legend()

        # Third subplot -- XY
        ax3 = fig.add_subplot(133)
        ax3.scatter(gg_pi0_end[:,0], gg_pi0_end[:,1], c=colors[4], marker='o', label='Diphoton')
        ax3.scatter(ee_pi0_end[:,0], ee_pi0_end[:,1], c=colors[1], marker='x', label=r'$e^+ e^-$')
        ax3.scatter(g_pi0_end[:,0], g_pi0_end[:,1], c=colors[2], marker='*', label='Single Photon')
        ax3.scatter(nodec_pi0_end[:,0], nodec_pi0_end[:,1], c=colors[3], marker='^', label='No Decay Products')
        ax3.scatter(eeg_pi0_end[:,0], eeg_pi0_end[:,1], c=colors[0], marker='v', label=r'Single Photon $e^+ e^-$')
        ax3.set_xlabel('Drift Axis [cm]')
        ax3.set_ylabel('Vertical Axis [cm]')
        ax3.set_title(r"$\pi^0$ End Location by Decay Product")
        #ax3.legend()

        output.savefig(fig)
        plt.close(fig)

        # PLOT: pi0 start direction by decay product (wrt to beam axis)
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111)
        print("gg_pi0_start_mom_unit_vec:", [np.array(gg[key]['pi0_start_mom']/np.sqrt(np.sum(np.array(gg[key]['pi0_start_mom'])**2))) for key in gg.keys()][:5])
        print("gg_pi0_start_dir_a_dot_dot:", [np.sum(np.array(gg[key]['pi0_start_mom']/np.sqrt(np.sum(np.array(gg[key]['pi0_start_mom'])**2)))*np.array([0.,0.,1.])) for key in gg.keys()][:5])
        gg_pi0_start_dir = np.array([np.sum(np.array(gg[key]['pi0_start_mom']/np.sqrt(np.sum(np.array(gg[key]['pi0_start_mom'])**2)))*np.array([0.,0.,1.])) for key in gg.keys()])
        print("gg_pi0_start_dir: ", gg_pi0_start_dir[:5])
        ggpi0_start_dir_counts, ggpi0_start_dir_bins = np.histogram(gg_pi0_start_dir, bins=np.linspace(-1, 1,41))
        ax.hist(ggpi0_start_dir_bins[:-1], bins=ggpi0_start_dir_bins, weights = ggpi0_start_dir_counts*scale_factor, histtype='stepfilled', alpha=0.45, color=colors[4], label='Diphoton', linestyle='-', linewidth=1.5, edgecolor=colors[4])
        ee_pi0_start_dir = np.array([np.sum(np.array(ee[key]['pi0_start_mom']/np.sqrt(np.sum(np.array(ee[key]['pi0_start_mom'])**2)))*np.array([0.,0.,1.])) for key in ee.keys()])
        eepi0_start_dir_counts, eepi0_start_dir_bins = np.histogram(ee_pi0_start_dir, bins=np.linspace(-1, 1,41))
        ax.hist(eepi0_start_dir_bins[:-1], bins=eepi0_start_dir_bins, weights = eepi0_start_dir_counts*scale_factor, histtype='stepfilled', alpha=0.45, color=colors[1], label=r'$e^+ e^-$', linestyle='-.', linewidth=1.5, edgecolor=colors[1])
        g_pi0_start_dir = np.array([np.sum(np.array(g[key]['pi0_start_mom']/np.sqrt(np.sum(np.array(g[key]['pi0_start_mom'])**2)))*np.array([0.,0.,1.])) for key in g.keys()])
        gpi0_start_dir_counts, gpi0_start_dir_bins = np.histogram(g_pi0_start_dir, bins=np.linspace(-1, 1,41))
        ax.hist(gpi0_start_dir_bins[:-1], bins=gpi0_start_dir_bins, weights = gpi0_start_dir_counts*scale_factor, histtype='stepfilled', alpha=0.45, color=colors[2], label='Single Photon', linestyle='--', linewidth=1.5, edgecolor=colors[2])
        nodec_pi0_start_dir = np.array([np.sum(np.array(nodec[key]['pi0_start_mom']/np.sqrt(np.sum(np.array(nodec[key]['pi0_start_mom'])**2)))*np.array([0.,0.,1.])) for key in nodec.keys()])
        nodecpi0_start_dir_counts, nodecpi0_start_dir_bins = np.histogram(nodec_pi0_start_dir, bins=np.linspace(-1, 1,41))
        ax.hist(nodecpi0_start_dir_bins[:-1], bins=nodecpi0_start_dir_bins, weights = nodecpi0_start_dir_counts*scale_factor, histtype='stepfilled', alpha=0.45, color=colors[3], label='No Decay Products', linestyle=':', linewidth=1.5, edgecolor=colors[3])
        eeg_pi0_start_dir = np.array([np.sum(np.array(eeg[key]['pi0_start_mom']/np.sqrt(np.sum(np.array(eeg[key]['pi0_start_mom'])**2)))*np.array([0.,0.,1.])) for key in eeg.keys()])
        eegpi0_start_dir_counts, eegpi0_start_dir_bins = np.histogram(eeg_pi0_start_dir, bins=np.linspace(-1, 1,41))
        ax.hist(eegpi0_start_dir_bins[:-1], bins=eegpi0_start_dir_bins, weights = eegpi0_start_dir_counts*scale_factor, histtype='stepfilled', alpha=0.45, color=colors[0], label=r'Single Photon $e^+ e^-$', linestyle='-.', linewidth=1.5, edgecolor=colors[0])
        ax.set_xlabel(r"Cosine of $\pi^0$ Start Direction w.r.t. Beam Axis")
        ax.set_ylabel("Count / 0.05")
        ax.set_title(r"$\pi^0$ Start Direction by Decay Product")
        ax.legend()
        output.savefig(fig)
        plt.close(fig)






        # PLOT: pi0 

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