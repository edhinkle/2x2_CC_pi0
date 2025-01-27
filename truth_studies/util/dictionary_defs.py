################################################################################
##                                                                            ##
##    CONTAINS: Methods to create python dictionaries with information to     ##
##              describe general signal events (get_signal_dict, muons in     ##
##              signal events (muon_characterization), and hadrons in signal  ##
##              events (hadron_characterization). Modified from               ##
##              https://github.com/edhinkle/mesonless_numubarCC/blob/main/    ##
##              truth_kinematics/file_parsing/signal_characterization.py      ##
##                                                                            ##
################################################################################

import sys
sys.path.append('../../')
import truth_studies.util.location_class as loc_class
import truth_studies.util.particlePDG_defs as pdg_defs
import truth_studies.util.truth_ixn_methods as truth
import truth_studies.util.singleParticleAssociation_methods as particle_assoc
import truth_studies.util.kinematicVariable_methods as kinematics
import numpy as np

''' TO DO: Add other hadron mult over threshold? '''


'''
NOTE: While script is labelled as 'signal' characterization, methods can be used for
characterizing certain backgrounds as well.

INCLUDED METHODS:
 - get_truth_dict(spill_id, vert_id, ghdr, gstack, traj, seg, signal_dict)
 - muon_characterization(spill_id, vert_id, ghdr, gstack, traj, seg, muon_dict)
 - hadron_characterization(spill_id, vert_id, ghdr, gstack, traj, seg, hadron_dict)
 '''


''' Purpose: Fill Python dictionary with basic GENIE truth-level vertex information
             and 
    Inputs : Spill ID (INT), Vertex ID (INT), genie_hdr dataset (HDF5 DATASET), 
             genie_stack dataset (HDF5 DATASET), edep-sim trajectories dataset (HDF5 DATASET), 
             vertex dataset (HDF5 DATASET), edep-sim segements dataset (HDF5 DATASET), 
             empty Python dictionary (DICT)
    Outputs: Nothing returned, but signal_dict (DICT) is full after
             method runs'''
def get_truth_dict(spill_id, vert_id, ghdr, gstack, traj, vert, seg, signal_dict):

    ghdr_vert_mask = ghdr['vertex_id']==vert_id
    truth_level_summ = ghdr[ghdr_vert_mask]

    mom = truth_level_summ['lep_mom'] # Truth-level outgoing muon momentum
    ang = truth_level_summ['lep_ang'] *np.pi / 180. # Truth-level muon angle with beam
    nu_energy = truth_level_summ['Enu'] # Truth-level neutrino energy
    q2 = truth_level_summ['Q2'] # Truth-level interaction 4-momentum squared
    #print(truth_level_summ.dtype.names)
    try:
        vtx = truth_level_summ['vertex'] # Truth-level vertex information
        vtx_x = vtx[0][0]
        vtx_y = vtx[0][1]
        vtx_z = vtx[0][2]
    except:
        vtx_x = truth_level_summ['x_vert']
        vtx_y = truth_level_summ['y_vert']
        vtx_z = truth_level_summ['z_vert']
                                                                                                                                          
    signal_dict[(spill_id,vert_id)]=dict(
        nu_energy=float(nu_energy),
        q2 = float(q2),
        mom=float(mom), 
        ang=float(ang),
        vtx_x = float(vtx_x), 
        vtx_y = float(vtx_y), 
        vtx_z = float(vtx_z))
    return

''' Purpose: Fill Python dictionary with overall truth-level information
             and 
    Inputs : Spill ID (INT), Vertex ID (INT), genie_hdr dataset (HDF5 DATASET), 
             genie_stack dataset (HDF5 DATASET), edep-sim trajectories dataset (HDF5 DATASET), 
             vertex dataset (HDF5 DATASET), edep-sim segements dataset (HDF5 DATASET), 
             empty Python dictionary (DICT)
    Outputs: Nothing returned, but signal_dict (DICT) is full after
             method runs'''
def get_truth_sig_bkg_dict(spill_id, vert_id, ghdr, gstack, traj, vert, seg, sig_bkg_dict, sim_file, flow_event_id, file_number, sig_or_bkg, is_cc):

    # Truth level summary information
    ghdr_vert_mask = ghdr['vertex_id']==vert_id
    truth_level_summ = ghdr[ghdr_vert_mask]

    nu_energy = truth_level_summ['Enu'] # Truth-level neutrino energy
    q2 = truth_level_summ['Q2'] # Truth-level interaction 4-momentum squared
    nu_pdg = truth_level_summ['nu_pdg']
    nu_int_type = truth.nu_int_type(ghdr, vert_id) # Neutrino interaction mechanism

    try:
        vtx = truth_level_summ['vertex'] # Truth-level vertex information
        vtx_x = vtx[0][0]
        vtx_y = vtx[0][1]
        vtx_z = vtx[0][2]
    except:
        vtx_x = truth_level_summ['x_vert']
        vtx_y = truth_level_summ['y_vert']
        vtx_z = truth_level_summ['z_vert']

    if (abs(nu_pdg) == 14):
        muon_momentum = truth_level_summ['lep_mom']
        muon_angle = truth_level_summ['lep_ang'] *np.pi / 180. # Truth-level muon angle with beam
    else:
        muon_momentum = -99999999. # Initialization; will change once we go through final states
        muon_angle = -99999999. # Initialization; will change once we go through final states

    # Get primary particles associated with vertex
    gstack_vert_mask = gstack['vertex_id']==vert_id
    gstack_vert = gstack[gstack_vert_mask] # Particle ID information associated with vertex

    # Get final state particles associated with vertex and their multiplicities
    gstack_vert_fs_mask = gstack_vert['part_status']==1 # Excludes initial state particles
    gstack_vert_fs = gstack_vert[gstack_vert_fs_mask]['part_pdg'] # Final state particle PDG IDs
    gstack_vert_fs_traj = gstack_vert[gstack_vert_fs_mask]['traj_id'] # Final state particle trajectory IDs
    gstack_vert_fs_abs = [abs(fsp) for fsp in gstack_vert_fs] # Absolute value of final state particle PDG IDs
    muon_mult = gstack_vert_fs_abs.count(13) # Muon multiplicity
    pi0_mult = gstack_vert_fs_abs.count(111) # Pi0 multiplicity
    primary_photon_mult = gstack_vert_fs_abs.count(22) # Photon multiplicity
    primary_electron_mult = gstack_vert_fs_abs.count(11) # Electron multiplicity
    primary_proton_mult = gstack_vert_fs_abs.count(2212) # Proton multiplicity
    primary_charged_pion_mult = gstack_vert_fs_abs.count(211) # Charged pion multiplicity
    primary_neutron_mult = gstack_vert_fs_abs.count(2112) # Neutron multiplicity
    primary_charged_kaon_mult = gstack_vert_fs_abs.count(321) # Charged kaon multiplicity
    primary_eta_mult = gstack_vert_fs_abs.count(221) # Eta multiplicity

    # Get trajectories associated with vertex
    traj_vert_mask = traj['vertex_id']==vert_id 
    final_states = traj[traj_vert_mask] # Get trajectories associated with vertex

    # Get primary pi0 and eta information
    gstack_pi0_mask = gstack_vert['part_pdg']==111
    gstack_pi0_traj = gstack_vert[gstack_pi0_mask]['traj_id']
    if len(gstack_pi0_traj) > 0:
        pi0_traj_mask = np.where(np.isin(final_states['traj_id'], gstack_pi0_traj))
        fs_pi0s = final_states[pi0_traj_mask] # Get list of pi0s associated with vertex
        pi0_start_momentum = fs_pi0s['pxyz_start']
        pi0_start_location = fs_pi0s['xyz_start']
        pi0_ke = fs_pi0s['E_start'] - pdg_defs.rest_mass_dict[111]
    else:
        pi0_start_momentum = []
        pi0_start_location = []
        pi0_ke = []

    gstack_eta_mask = gstack_vert['part_pdg']==221
    gstack_eta_traj = gstack_vert[gstack_eta_mask]['traj_id']
    if len(gstack_eta_traj) > 0:
        eta_traj_mask = final_states['traj_id']==gstack_eta_traj
        fs_etas = final_states[eta_traj_mask] # Get list of etas associated with vertex
        eta_start_momentum = fs_etas['pxyz_start']
        eta_start_location = fs_etas['xyz_start']
        eta_ke = fs_etas['E_start'] - pdg_defs.rest_mass_dict[221]
    else:
        eta_start_momentum = []
        eta_start_location = []
        eta_ke = []

    # Initialize lists for pi0 and eta children, primary showers
    pi0_child_pdg = [[] for i in range(pi0_mult)]
    pi0_child_edep = [[] for i in range(pi0_mult)]
    pi0_child_available_energy = [[] for i in range(pi0_mult)]
    pi0_child_start_location = [[] for i in range(pi0_mult)]
    pi0_child_end_location = [[] for i in range(pi0_mult)]
    pi0_child_start_pxyz = [[] for i in range(pi0_mult)]
    pi0_child_conv_dist = [[] for i in range(pi0_mult)]

    eta_child_pdg = [[] for i in range(primary_eta_mult)]
    eta_child_edep = [[] for i in range(primary_eta_mult)]
    eta_child_available_energy = [[] for i in range(primary_eta_mult)]
    eta_child_start_location = [[] for i in range(primary_eta_mult)]
    eta_child_end_location = [[] for i in range(primary_eta_mult)]
    eta_child_start_pxyz = [[] for i in range(primary_eta_mult)]
    eta_child_conv_dist = [[] for i in range(primary_eta_mult)]

    primary_shower_pdg = []
    primary_shower_edep = []
    primary_shower_available_energy = []
    primary_shower_start_location = []
    primary_shower_end_location = []
    primary_shower_start_pxyz = []
    primary_shower_conv_dist = []


    exclude_track_ids = set(gstack_pi0_traj) | set(gstack_eta_traj) # Create set of track IDs to exclude to eliminate redundancies

    # Loop through final states to get relevant information
    for fs in final_states:
        
        # Get track id/figure out if it is excluded
        track_id = fs['traj_id']
        if track_id in exclude_track_ids: continue
        else: 
            start_tid = fs['traj_id']

        # Figure out track PDG and figure out what to do next accordingly
        pdg = fs['pdg_id'] # *** pdg ***   

        # If the track is a muon, get its momentum, and endpoint containment
        if abs(pdg) == 13:

            # Get track ID set associated with muon
            track_id_set = particle_assoc.same_pdg_connected_trajectories(pdg, start_tid, final_states, traj, ghdr)
            exclude_track_ids.update(track_id_set)

            # Make sure the muon is a primary particle
            is_primary = truth.is_primary_particle(track_id_set, final_states, traj, ghdr)
            if is_primary == False: continue

            if abs(nu_pdg) == 12:

                # Get muon momentum and angle for electron neutrino interactions
                tid_at_vertex = particle_assoc.find_trajectory_at_vertex(track_id_set, final_states, traj, ghdr)
                tid_at_vertex_mask = final_states['traj_id'] == tid_at_vertex
                fs_at_vertex = final_states[tid_at_vertex_mask]
                muon_momentum_three_vector = fs_at_vertex['pxyz_start']
                muon_momentum = np.linalg.norm(muon_momentum_three_vector)
                muon_angle = np.arccos(np.dot(kinematics.beam_direction, muon_momentum_three_vector) / muon_momentum)

            # Characterize Muon Endpoint/Containment
            if len(track_id_set)>1:
                start_tid = particle_assoc.find_trajectory_at_vertex(track_id_set, final_states,traj, ghdr)
                start_tid_mask = final_states['traj_id']==start_tid
                muon_start_traj = final_states[start_tid_mask]
                start_pt = muon_start_traj['xyz_start']

                end_tid = particle_assoc.find_forward_primary_particle_end_trajectory(track_id_set, final_states,traj, ghdr)
                end_tid_mask = final_states['traj_id']==end_tid
                muon_end_traj = final_states[end_tid_mask]
                end_pt = muon_end_traj['xyz_end']
            else:
                end_pt = fs['xyz_end']
                start_pt = fs['xyz_start']
            muon_end_pt_loc = loc_class.particle_end_loc(start_pt, end_pt)

        # If the track is a shower particle (photon or electron/positron)
        elif abs(pdg) == 11 or pdg == 22:

            # Primary shower particles
            if start_tid in gstack_vert_fs_traj:

                primary_shower_pdg.append(pdg)

                # Get track ID set associated with primary shower particle
                initial_track_id_set = particle_assoc.same_pdg_connected_trajectories(pdg, start_tid, final_states, traj, ghdr)
                track_id_set = particle_assoc.connect_shower_traj(pdg, track_id, final_states, traj, ghdr)
                exclude_track_ids.update(track_id_set) # Exclude track IDs associated with same particle from future counting

                if pdg == 22:
                    total_edep = fs['E_start'] # total available energy at start of shower
                elif abs(pdg) == 11:
                    total_edep = fs['E_start'] - pdg_defs.rest_mass_dict[abs(pdg)] # total available energy at start of shower
                primary_shower_available_energy.append(total_edep)
        
                contained_edep = 0.
                for tid in track_id_set:
                
                    contained_edep+= kinematics.fv_edep_e(tid, vert_id, seg)
                
                primary_shower_edep.append(contained_edep)
        
                # Characterize shower endpoints
                if len(initial_track_id_set)>1:
                    start_tid_mask = final_states['traj_id']==start_tid
                    start_traj = final_states[start_tid_mask]
                    start_pt = start_traj['xyz_start']
                    start_pxyz = start_traj['pxyz_start']
                    if len(start_pt) > 1:
                        print("Multiple start points found for trajectory ID:", start_tid)
                    else:
                        start_pt = start_pt[0]
                        start_pxyz = start_pxyz[0]
                    
                    end_tid = particle_assoc.find_secondary_particle_end_trajectory(initial_track_id_set, final_states,start_tid)
                    end_tid_mask = final_states['traj_id']==end_tid
                    end_traj = final_states[end_tid_mask]
                    end_pt = end_traj['xyz_end']
                    if len(end_pt) > 1:
                        print("Multiple end points found for trajectory ID:", end_tid)
                    else:
                        end_pt = end_pt[0]
                else:
                    end_tid = start_tid
                    end_pt = fs['xyz_end']
                    start_pt = fs['xyz_start']
                    start_pxyz = fs['pxyz_start']

                primary_shower_conv_dist.append(particle_assoc.distance_between_points(start_pt, end_pt))
                primary_shower_start_location.append(start_pt)
                primary_shower_end_location.append(end_pt)
                primary_shower_start_pxyz.append(start_pxyz)

            # Particle is potentially a pi0 child
            else:

                # Choose trajectories associated with pi0s as parents
                if (fs['parent_id'] in gstack_pi0_traj):

                    fs_parent_id = fs['parent_id']
                    pi0_idx = np.argwhere(gstack_pi0_traj == fs_parent_id)[0][0]

                    pi0_child_pdg[pi0_idx].append(pdg)   

                    initial_child_track_id_set = particle_assoc.same_pdg_connected_trajectories(pdg, start_tid, final_states, traj, ghdr)

                    track_id_set = particle_assoc.connect_shower_traj(pdg, track_id, final_states, traj, ghdr)
                    exclude_track_ids.update(track_id_set) # Exclude track IDs associated with same particle from future counting

                    if pdg == 22:
                        total_edep = fs['E_start'] # total available energy at start of shower
                    elif abs(pdg) == 11:
                        total_edep = fs['E_start'] - pdg_defs.rest_mass_dict[abs(pdg)] # total available energy at start of shower

                    contained_edep = 0.
                    for tid in track_id_set:
                    
                        contained_edep+= kinematics.fv_edep_e(tid, vert_id, seg)

                    pi0_child_edep[pi0_idx].append(contained_edep)
                    pi0_child_available_energy[pi0_idx].append(total_edep)

                    # Characterize pi0 decay product Endpoints/Containments
                    if len(initial_child_track_id_set)>1:
                        start_tid_mask = final_states['traj_id']==start_tid
                        child_start_traj = final_states[start_tid_mask]
                        start_pt = child_start_traj['xyz_start']
                        start_pxyz = child_start_traj['pxyz_start']
                        if len(start_pt) > 1:
                            print("Multiple start points found for trajectory ID:", start_tid)
                        else:
                            start_pt = start_pt[0]
                            start_pxyz = start_pxyz[0]

                        end_tid = particle_assoc.find_secondary_particle_end_trajectory(track_id_set, final_states,start_tid)
                        end_tid_mask = final_states['traj_id']==end_tid
                        child_end_traj = final_states[end_tid_mask]
                        end_pt = child_end_traj['xyz_end']
                        if len(end_pt) > 1:
                            print("Multiple end points found for trajectory ID:", end_tid)
                        else:
                            end_pt = end_pt[0]
                    else:
                        end_tid = start_tid
                        end_pt = fs['xyz_end']
                        start_pt = fs['xyz_start']
                        start_pxyz = fs['pxyz_start']

                    pi0_child_conv_dist[pi0_idx].append(particle_assoc.distance_between_points(start_pt, end_pt))
                    pi0_child_start_location[pi0_idx].append(start_pt)
                    pi0_child_end_location[pi0_idx].append(end_pt)
                    pi0_child_start_pxyz[pi0_idx].append(start_pxyz)

        # Choose trajectories associated with etas as parents
        if (fs['parent_id'] in gstack_eta_traj):

            fs_parent_id = fs['parent_id']
            eta_idx = np.argwhere(gstack_eta_traj == fs_parent_id)[0][0]

            eta_child_pdg[eta_idx].append(pdg)
            print("Eta Child PDG:", pdg)
        

    if sig_or_bkg == 1:
        sig_bkg_label = r"$\nu_{\mu}$ CC 1$\pi^0$, 0$\pi^{\pm}$"
    else:
        if is_cc == 1 and abs(nu_pdg) == 14:
            if pi0_mult == 0:
                if primary_photon_mult == 1:
                    sig_bkg_label = r"$\nu_{\mu}$ CC 0$\pi^0$, 1$\gamma$"
                elif primary_photon_mult == 2:
                    sig_bkg_label = r"$\nu_{\mu}$ CC 0$\pi^0$, 2$\gamma$"
                elif primary_photon_mult > 2:
                    sig_bkg_label = r"$\nu_{\mu}$ CC 0$\pi^0$, >2$\gamma$"
                else:
                    sig_bkg_label = r"$\nu_{\mu}$ CC 0$\pi^0$, 0$\gamma$"
            elif pi0_mult == 1:
                if primary_charged_pion_mult > 0:
                    sig_bkg_label = r"$\nu_{\mu}$ CC 1$\pi^0$, N$\pi^{\pm}$"
                else:
                    sig_bkg_label = r"$\nu_{\mu}$ CC 1$\pi^0$, 0$\pi^{\pm}$, OTHER"
            elif pi0_mult == 2:
                if primary_charged_pion_mult == 0:
                    sig_bkg_label = r"$\nu_{\mu}$ CC 2$\pi^0$, 0$\pi^{\pm}$"
                elif primary_charged_pion_mult >= 1:
                    sig_bkg_label = r"$\nu_{\mu}$ CC 2$\pi^0$, N$\pi^{\pm}$"
            elif pi0_mult > 2:
                if primary_charged_pion_mult == 0:
                    sig_bkg_label = r"$\nu_{\mu}$ CC >2$\pi^0$, 0$\pi^{\pm}$"
                elif primary_charged_pion_mult >= 1:
                    sig_bkg_label = r"$\nu_{\mu}$ CC >2$\pi^0$, N$\pi^{\pm}$"
            else:
                if primary_charged_pion_mult == 0:
                    sig_bkg_label = r"$\nu_{\mu}$ CC >2$\pi^0$, 0$\pi^{\pm}$"
                elif primary_charged_pion_mult >= 1:
                    sig_bkg_label = r"$\nu_{\mu}$ CC >2$\pi^0$, N$\pi^{\pm}$"
        elif is_cc == 1 and abs(nu_pdg) == 12:
            sig_bkg_label = r"$\nu_{e}$ CC"      
        else:
            sig_bkg_label = "NC Background"
                                                                                            
    sig_bkg_dict[(spill_id,vert_id)]=dict(
        sig = bool(sig_or_bkg),
        sig_bkg_label = str(sig_bkg_label),
        is_cc = bool(is_cc),
        nu_pdg = int(nu_pdg),
        muon_angle = float(muon_angle),
        muon_momentum = float(muon_momentum),
        muon_end_pt_loc = str(muon_end_pt_loc),
        nu_int_type = str(nu_int_type),
        muon_mult = int(muon_mult),
        pi0_mult = int(pi0_mult),
        primary_photon_mult = int(primary_photon_mult),
        primary_electron_mult = int(primary_electron_mult),
        primary_proton_mult = int(primary_proton_mult),
        primary_charged_pion_mult = int(primary_charged_pion_mult),
        primary_neutron_mult = int(primary_neutron_mult),
        primary_charged_kaon_mult = int(primary_charged_kaon_mult),
        primary_eta_mult = int(primary_eta_mult),
        pi0_start_momentum = [[float(j) for j in pi0_start_momentum[i]] for i in range(len(pi0_start_momentum))],
        pi0_start_location = [[float(j) for j in pi0_start_location[i]] for i in range(len(pi0_start_location))],
        pi0_ke = [float(i) for i in pi0_ke],
        eta_start_momentum = [[float(j) for j in eta_start_momentum[i]] for i in range(len(eta_start_momentum))],
        eta_start_location = [[float(j) for j in eta_start_location[i]] for i in range(len(eta_start_location))],
        eta_ke = [float(i) for i in eta_ke],
        pi0_child_pdg = [[int(j) for j in pi0_child_pdg[i]] for i in range(len(pi0_child_pdg))],
        pi0_child_edep = [[float(j) for j in pi0_child_edep[i]] for i in range(len(pi0_child_edep))],
        pi0_child_available_energy = [[float(j) for j in pi0_child_available_energy[i]] for i in range(len(pi0_child_available_energy))],
        pi0_child_start_location = [[[float(k) for k in j] for j in pi0_child_start_location[i]] for i in range(len(pi0_child_start_location))],
        pi0_child_end_location = [[[float(k) for k in j] for j in pi0_child_end_location[i]] for i in range(len(pi0_child_end_location))],
        pi0_child_start_pxyz = [[[float(k) for k in j] for j in pi0_child_start_pxyz[i]] for i in range(len(pi0_child_start_pxyz))],
        pi0_child_conv_dist = [[float(j) for j in pi0_child_conv_dist[i]] for i in range(len(pi0_child_conv_dist))],
        eta_child_pdg = [[int(j) for j in eta_child_pdg[i]] for i in range(len(eta_child_pdg))],
        #eta_child_edep = [[float(j) for j in eta_child_edep[i]] for i in range(len(eta_child_edep))],
        #eta_child_available_energy = [[float(j) for j in eta_child_available_energy[i]] for i in range(len(eta_child_available_energy))],
        #eta_child_start_location = [[[float(k) for k in j] for j in eta_child_start_location[i]] for i in range(len(eta_child_start_location))],
        #eta_child_end_location = [[[float(k) for k in j] for j in eta_child_end_location[i]] for i in range(len(eta_child_end_location))],
        #eta_child_start_pxyz = [[[float(k) for k in j] for j in eta_child_start_pxyz[i]] for i in range(len(eta_child_start_pxyz))],
        #eta_child_conv_dist = [[[float(k) for k in j] for j in eta_child_conv_dist[i]] for i in range(len(eta_child_conv_dist))],
        primary_shower_start_location = [[float(j) for j in primary_shower_start_location[i]] for i in range(len(primary_shower_start_location))],
        primary_shower_end_location = [[float(j) for j in primary_shower_end_location[i]] for i in range(len(primary_shower_end_location))],
        primary_shower_pdg = [int(i) for i in primary_shower_pdg],
        primary_shower_edep = [float(i) for i in primary_shower_edep],
        primary_shower_available_energy = [float(i) for i in primary_shower_available_energy],
        primary_shower_start_pxyz = [[float(j) for j in primary_shower_start_pxyz[i]] for i in range(len(primary_shower_start_pxyz))],
        primary_shower_conv_dist = [float(i) for i in primary_shower_conv_dist],
        nu_energy=float(nu_energy),
        q2 = float(q2),
        vtx_x = float(vtx_x), 
        vtx_y = float(vtx_y), 
        vtx_z = float(vtx_z), 
        filepath = str(sim_file),
        flow_event_id = int(flow_event_id),
        file_number = int(file_number)
        )
    return


''' Purpose: Fill Python dictionary with basic GENIE and edep-sim truth-level vertex information
             related to FS muon
    Inputs : Spill ID (INT), Vertex ID (INT), genie_hdr dataset (HDF5 DATASET), 
             genie_stack dataset (HDF5 DATASET), edep-sim trajectories dataset (HDF5 DATASET), 
             vertex dataset (HDF5 DATASET), edep-sim segements dataset (HDF5 DATASET), 
             empty Python dictionary (DICT)
    Outputs: Nothing returned, but muon_dict (DICT) is full after
             method runs'''
def muon_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, muon_dict):

    traj_vert_mask = traj['vertex_id']==vert_id 
    final_states = traj[traj_vert_mask] # Get trajectories associated with vertex

    ghdr_vert_mask = ghdr['vertex_id']==vert_id
    truth_level_summ = ghdr[ghdr_vert_mask] # Get GENIE truth info associated with vertex

    mom = truth_level_summ['lep_mom'] # Truth-level outgoing muon momentum
    ang = truth_level_summ['lep_ang'] *np.pi / 180. # Truth-level muon angle with beam
    nu_energy = truth_level_summ['Enu'] # Truth-level neutrino energy
    q2 = truth_level_summ['Q2'] # Truth-level interaction 4-momentum squared
    nu_int_type = truth.nu_int_type(ghdr, vert_id) # Neutrino interaction mechanism
    nu_pdg = truth_level_summ['nu_pdg']

    total_edep=0.; contained_edep=0.; total_length=0.; contained_length=0. # Set contained and total track energies and lengths to 0

    gstack_vert_mask = gstack['vertex_id']==vert_id
    gstack_pdg_set = set(gstack[gstack_vert_mask]['part_pdg']) # Get set of PDG IDs for particles associated with vertex

    exclude_track_ids = set() # Create set of track IDs to exclude to eliminate redundancies
                              # (i.e. if they've been identified as being from the same particle as an earlier track)
    for fs in final_states:

        # Choose nu_mu_bar or nu_mu vertices
        if (abs(fs['pdg_id']) != 13): continue
        
        track_id = fs['traj_id']
        if track_id in exclude_track_ids: continue

        pdg = fs['pdg_id'] # *** pdg ***     

        track_id_set = particle_assoc.same_pdg_connected_trajectories(pdg, track_id, final_states, traj, ghdr)
        exclude_track_ids.update(track_id_set) # Exclude track IDs associated with same particle from future counting

        is_primary = truth.is_primary_particle(track_id_set, final_states, traj, ghdr) 

        if is_primary == False: continue # Only look at final state particles

        track_id_at_vertex = particle_assoc.find_trajectory_at_vertex(track_id_set, final_states,traj, ghdr)
        final_state_vertex_tid_mask = final_states['traj_id'] == track_id_at_vertex
        fs_at_vertex = final_states[final_state_vertex_tid_mask]

        parent_pdg = truth.find_parent_pdg(fs_at_vertex['parent_id'],vert_id, traj, ghdr)# *** parent pdg ***

        if len(track_id_set)>1:
            print("Length of Track ID Set:", len(track_id_set))

        for tid in track_id_set:

            total_edep += kinematics.total_edep_charged_e(tid,vert_id,seg) # *** total visible energy ***
            contained_edep+= kinematics.fv_edep_charged_e(tid, vert_id, seg)
            
            contained_length+=kinematics.fv_edep_charged_length(tid, vert_id, seg) # *** total visible track length ***
            total_length+=kinematics.total_edep_charged_length(tid, vert_id, seg)
        
        # Characterize Muon Endpoint/Containment
        if len(track_id_set)>1:
            start_tid = particle_assoc.find_trajectory_at_vertex(track_id_set, final_states,traj, ghdr)
            start_tid_mask = final_states['traj_id']==start_tid
            muon_start_traj = final_states[start_tid_mask]
            start_pt = muon_start_traj['xyz_start']
            
            end_tid = particle_assoc.find_forward_primary_particle_end_trajectory(track_id_set, final_states,traj, ghdr)
            end_tid_mask = final_states['traj_id']==end_tid
            muon_end_traj = final_states[end_tid_mask]
            end_pt = muon_end_traj['xyz_end']
        else:
            end_pt = fs['xyz_end']
            start_pt = fs['xyz_start']
        end_pt_loc = loc_class.particle_end_loc(start_pt, end_pt)

    # Save collected info in input muon_dict                                                                                                                         
    muon_dict[(spill_id,vert_id)]=dict(
        pdg=int(pdg),
        nu_pdg = int(nu_pdg), 
        parent_pdg=int(parent_pdg),
        total_edep=float(total_edep),
        contained_edep=float(contained_edep),
        total_length=float(total_length),
        contained_length=float(contained_length),
        mom=float(mom), 
        ang=float(ang),
        nu_energy=float(nu_energy),
        q2 = float(q2),
        end_pt_loc = str(end_pt_loc),
        muon_start = [float(i) for i in start_pt],
        muon_end = [float(i) for i in end_pt],
        nu_int_type=str(nu_int_type))
    return


''' Purpose: Fill Python dictionary with basic GENIE and edep-sim truth-level vertex information
             related to hadrons 
    Inputs : Spill ID (INT), Vertex ID (INT), genie_hdr dataset (HDF5 DATASET), 
             genie_stack dataset (HDF5 DATASET), edep-sim trajectories dataset (HDF5 DATASET), 
             vertex dataset (HDF5 DATASET), edep-sim segements dataset (HDF5 DATASET), threshold length in cm (FLOAT)
             empty Python dictionary (DICT)
    Outputs: Nothing returned, but hadron_dict (DICT) is full after
             method runs'''
def hadron_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, threshold, hadron_dict):
        
    #print("\nHADRONS:")
    traj_vert_mask = traj['vertex_id']==vert_id
    final_states = traj[traj_vert_mask] # Trajectories associated with vertex

    leptons_abs_pdg = [11, 12, 13, 14, 15, 16] # List of lepton PDG IDs (abs value)

    ghdr_vert_mask = ghdr['vertex_id']==vert_id
    truth_level_summ = ghdr[ghdr_vert_mask] # Get GENIE truth info associated with vertex
    nu_pdg = truth_level_summ['nu_pdg']
    
    gstack_vert_mask = gstack['vertex_id']==vert_id
    gstack_vert = gstack[gstack_vert_mask] # Particle ID information associated with vertex

    gstack_vert_fs_mask = gstack_vert['part_status']==1 # Excludes initial state particles
    gstack_vert_fs = gstack_vert[gstack_vert_fs_mask]['part_pdg'] # Final state particle PDG IDs

    gstack_vert_fs_hadrons = [fsp for fsp in gstack_vert_fs if abs(fsp) not in leptons_abs_pdg and fsp != 22 and fsp != 111] # LIST of f.s. hadron PDG IDs 
    gstack_vert_fs_pdg_set = set(gstack_vert_fs_hadrons) # SET of f.s. hadron PDG IDs

    #print("Vertex PDG Stack:", gstack_vert['part_pdg'])
    #print("Final State Hadrons:", gstack_vert_fs_hadrons)
    #print("Final State PDG Stack:", gstack_vert_fs)
    #print("Vertex PDG Stack Set:", gstack_pdg_set)
    nu_int_type = truth.nu_int_type(ghdr, vert_id) # Neutrino interaction mechanism (Truth)

    hadron_mult = len(gstack_vert_fs_hadrons) # F.S. hadron multiplicity
    n_mult = 0
    p_mult = 0
    other_had_mult = 0

    # Get f.s. hadron sub-category multiplicities
    for had in range(hadron_mult):

        if gstack_vert_fs_hadrons[had] == 2112:
            n_mult+=1 # Neutron multiplicity
        elif gstack_vert_fs_hadrons[had] == 2212:
            p_mult+=1 # Proton multiplicity
        else:
            other_had_mult+=1

    # Set contained and total track energies and lengths to 0; Set hadron + proton multiplicities over threshold to 0.
    total_edep=0.; contained_edep=0.; total_length=0.; contained_length=0.
    total_edep_over_thresh=0.; contained_edep_over_thresh=0.
    max_proton_contained_length=0.; max_proton_total_length=0. 
    lead_proton_ang_wrt_beam=0.; lead_proton_momentum=0.; sub_lead_proton_ang_wrt_beam=0.; sub_lead_proton_momentum=0.
    lead_proton_traj_at_vertex = 0.; sub_lead_proton_traj_at_vertex = 0.
    hadron_mult_over_thresh = 0.; p_mult_over_thresh = 0.
    p_ke = 0.

    exclude_track_ids = set() # Create set of track IDs to exclude to eliminate redundancies
                              # (i.e. if they've been identified as being from the same particle as an earlier track)
    for fs in final_states:

        if abs(fs['pdg_id']) in leptons_abs_pdg: continue # No leptons
        if fs['pdg_id'] > 1000000000: continue # No nuclei
        if fs['pdg_id'] == 22: continue # No photons
        if fs['pdg_id'] == 111: continue # No pi0s
        
        track_id = fs['traj_id']
        if track_id in exclude_track_ids: 
            #print("Excluding a track because it has already been studied.")
            continue 

        track_id_set = particle_assoc.same_pdg_connected_trajectories(fs['pdg_id'], track_id, final_states, traj, ghdr)
        exclude_track_ids.update(track_id_set) # Exclude track IDs associated with same particle from future counting
        #if fs['pdg_id'] == 2112: print("\nTrack ID Set:", track_id_set)

        proton_contained_length = 0.; proton_total_length=0. # Reset proton track lengths
        fs_total_edep =0.; fs_contained_edep=0. # Reset Edep for individual final states
        #print("Final state PDG: ", fs['pdg_id'])
        for tid in track_id_set:

            total_edep_temp = kinematics.total_edep_charged_e(tid,vert_id,seg)
            contained_edep_temp = kinematics.fv_edep_charged_e(tid, vert_id, seg)
            fs_total_edep+=total_edep_temp
            total_edep += total_edep_temp# *** total visible energy ***
            fs_contained_edep+= contained_edep_temp
            contained_edep+= contained_edep_temp

            contained_length+=kinematics.fv_edep_charged_length(tid, vert_id, seg) # *** total contained length for protons ***
            total_length+=kinematics.total_edep_charged_length(tid, vert_id, seg)
            
            if fs['pdg_id'] == 2212:
                proton_contained_length+=kinematics.fv_edep_charged_length(tid, vert_id, seg) # *** total contained length for protons ***
                proton_total_length+=kinematics.total_edep_charged_length(tid, vert_id, seg)

        if truth.is_primary_particle(track_id_set, final_states, traj, ghdr) and contained_length > threshold \
            and fs['pdg_id'] not in pdg_defs.neutral_hadron_pdg_dict.keys():
            hadron_mult_over_thresh +=1
            if fs['pdg_id'] == 2212: 
                p_mult_over_thresh += 1
                p_traj_id_at_vertex = particle_assoc.find_trajectory_at_vertex(track_id_set, final_states, traj, ghdr)
                p_mom = kinematics.truth_primary_particle_momentum(track_id_set, final_states, traj, ghdr)
                p_ke += kinematics.truth_primary_particle_kinetic_energy(fs['pdg_id'],track_id_set, final_states, traj, ghdr)
                total_edep_over_thresh += fs_total_edep # *** total visible energy ***
                contained_edep_over_thresh+= fs_contained_edep
                if p_mom > lead_proton_momentum:
                    sub_lead_proton_traj_at_vertex = lead_proton_traj_at_vertex
                    sub_lead_proton_momentum = lead_proton_momentum
                    sub_lead_proton_ang_wrt_beam = lead_proton_ang_wrt_beam

                    lead_proton_traj_at_vertex = p_traj_id_at_vertex
                    lead_proton_momentum = p_mom
                    lead_proton_ang_wrt_beam = kinematics.angle_wrt_beam_direction(track_id_set, final_states, traj, ghdr)
                elif p_mom <= lead_proton_momentum and p_mom > sub_lead_proton_momentum:
                    sub_lead_proton_traj_at_vertex = p_traj_id_at_vertex
                    sub_lead_proton_momentum = p_mom
                    sub_lead_proton_ang_wrt_beam = kinematics.angle_wrt_beam_direction(track_id_set, final_states, traj, ghdr)
                
                if proton_contained_length > max_proton_contained_length:
                    max_proton_contained_length = proton_contained_length # Update max contained proton length in vertex
                    max_proton_total_length = proton_total_length # Update max total proton length in vertex

    angle_between_lead_and_sublead_protons = 0. # Initialize angle between leading and subleading protons
    if p_mult_over_thresh >= 2:
        angle_between_lead_and_sublead_protons = kinematics.angle_between_two_trajectories(lead_proton_traj_at_vertex, \
                                                                                          sub_lead_proton_traj_at_vertex, final_states)

    # Save collected info in input hadron_dict
    hadron_dict[(spill_id,vert_id)]=dict(
        nu_pdg = int(nu_pdg),
        hadron_mult = int(hadron_mult),
        neutron_mult = int(n_mult),
        proton_mult = int(p_mult),
        other_had_mult = int(other_had_mult),
        hadron_mult_over_thresh = int(hadron_mult_over_thresh),
        proton_mult_over_thresh = int(p_mult_over_thresh),
        hadron_pdg = [int(i) for i in gstack_vert_fs_hadrons],
        hadron_pdg_set = [int(i) for i in list(gstack_vert_fs_pdg_set)],
        total_edep=float(total_edep),
        contained_edep=float(contained_edep),
        max_p_total_length=float(max_proton_total_length),
        max_p_contained_length=float(max_proton_contained_length),
        lead_proton_momentum = float(lead_proton_momentum), 
        sub_lead_proton_momentum = float(sub_lead_proton_momentum), 
        lead_proton_ang_wrt_beam = float(lead_proton_ang_wrt_beam), 
        sub_lead_proton_ang_wrt_beam = float(sub_lead_proton_ang_wrt_beam),
        sub_lead_proton_angle_with_lead_proton = float(angle_between_lead_and_sublead_protons),
        primary_protons_total_ke = float(p_ke),
        total_edep_over_thresh = float(total_edep_over_thresh),
        contained_edep_over_thresh = float(contained_edep_over_thresh),
        nu_int_type=str(nu_int_type))
    return



''' Purpose: Fill Python dictionary with basic GENIE and edep-sim truth-level vertex information
             related to FS pi0(s)
    Inputs : Spill ID (INT), Vertex ID (INT), genie_hdr dataset (HDF5 DATASET), 
             genie_stack dataset (HDF5 DATASET), edep-sim trajectories dataset (HDF5 DATASET), 
             vertex dataset (HDF5 DATASET), edep-sim segements dataset (HDF5 DATASET), 
             empty Python dictionary (DICT)
    Outputs: Nothing returned, but pi0_dict (DICT) is full after
             method runs'''
def pi0_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, pi0_dict, sim_file, flow_event_id):

    traj_vert_mask = traj['vertex_id']==vert_id 
    final_states = traj[traj_vert_mask] # Get trajectories associated with vertex
    #print("Final States:", final_states['pdg_id'])

    ghdr_vert_mask = ghdr['vertex_id']==vert_id
    truth_level_summ = ghdr[ghdr_vert_mask] # Get GENIE truth info associated with vertex

    nu_energy = truth_level_summ['Enu'] # Truth-level neutrino energy
    q2 = truth_level_summ['Q2'] # Truth-level interaction 4-momentum squared
    nu_int_type = truth.nu_int_type(ghdr, vert_id) # Neutrino interaction mechanism
    nu_pdg = truth_level_summ['nu_pdg']
    #print("Reaction Type:", truth_level_summ['reaction'])


    gstack_vert_mask = gstack['vertex_id']==vert_id
    gstack_for_vertex = gstack[gstack_vert_mask] # Get list of PDG IDs for particles associated with vertex
    #print("Vertex PDG Stack:", gstack_for_vertex)
    gstack_pi0_mask = gstack_for_vertex['part_pdg']==111
    gstack_pi0_traj = gstack_for_vertex[gstack_pi0_mask]['traj_id'] #
    #print("Pi0 Trajectories:", gstack_pi0_traj)
    pi0_traj_mask = final_states['traj_id']==gstack_pi0_traj
    #print("Final States:", final_states['traj_id'])
    #print("Primary PDG:", final_states[final_states['parent_id']==-1]['pdg_id'])
    fs_pi0s = final_states[pi0_traj_mask] # Get list of pi0s associated with vertex
    #print("Pi0 Trajectory:", fs_pi0s)
    pi0_start_mom = fs_pi0s['pxyz_start']
    pi0_end_process = fs_pi0s['end_process']
    pi0_end_subprocess = fs_pi0s['end_subprocess']
    pi0_ke = fs_pi0s['E_start'] - pdg_defs.rest_mass_dict[111]
    if len(pi0_start_mom) > 1:
        print("Multiple start momenta found for pi0 trajectory ID:", gstack_pi0_traj)
    else:
        #print("Pi0 Start Momentum:", pi0_start_mom)
        pi0_start_mom = pi0_start_mom[0]
    pi0_end_loc = fs_pi0s['xyz_end']
    if len(pi0_end_loc) > 1:
        print("Multiple end points found for pi0 trajectory ID:", gstack_pi0_traj)
    else:
        pi0_end_loc = pi0_end_loc[0]
    pi0_start_loc = fs_pi0s['xyz_start']
    if len(pi0_start_loc) > 1:
        print("Multiple start points found for pi0 trajectory ID:", gstack_pi0_traj)
    else:
        pi0_start_loc = pi0_start_loc[0]

    exclude_track_ids = set(gstack_pi0_traj) # Create set of track IDs to exclude to eliminate redundancies
                              # (i.e. if they've been identified as being from the same particle as an earlier track)
    child_pdg = []
    child_total_edep = []
    child_contained_edep = []
    child_total_length = []
    child_contained_length = []
    child_start_pt = []
    child_end_pt = []
    child_end_pt_loc = []
    child_start_pxyz = []
    child_start_dir = []
    child_conv_dist = []
    child_fraction_traj_id_contained = []
    child_shower_electron_traj_id_fraction = []
                    
    for fs in final_states:

        total_edep=0.; contained_edep=0.; total_length=0.; contained_length=0. # Set contained and total track energies and lengths to 0

        # Choose trajectories associated with pi0s as parents
        if (fs['parent_id'] not in gstack_pi0_traj): continue
        #print("Found a pi0 child!")
        
        track_id = fs['traj_id']
        if track_id in exclude_track_ids: continue
        else: 
            start_tid = fs['traj_id']

        pdg = fs['pdg_id'] # *** pdg ***
        child_pdg.append(pdg)
        #print("Child PDG:", pdg)     

        initial_child_track_id_set = particle_assoc.same_pdg_connected_trajectories(pdg, start_tid, final_states, traj, ghdr)

        track_id_set = particle_assoc.connect_shower_traj(pdg, track_id, final_states, traj, ghdr)
        exclude_track_ids.update(track_id_set) # Exclude track IDs associated with same particle from future counting
        #print("Track ID Set:", track_id_set)
        #for track in track_id_set:
        #    print("Track PDG:", final_states[final_states['traj_id']==track]['pdg_id'])
        #    #print("Track Parent ID:", final_states[final_states['traj_id']==track]['parent_id'])
        #    print("Track End Process:", final_states[final_states['traj_id']==track]['end_process'])
        #    print("Track End Subprocess:", final_states[final_states['traj_id']==track]['end_subprocess'])
        #    print("Track Children:", final_states[final_states['parent_id']==track]['pdg_id'])
        #if len(track_id_set)>1:
        #    print("Length of Track ID Set:", len(track_id_set))

        if pdg == 22:
            total_edep = fs['E_start'] # total available energy at start of shower
        elif abs(pdg) == 11:
            total_edep = fs['E_start'] - pdg_defs.rest_mass_dict[abs(pdg)] # total available energy at start of shower

        for tid in track_id_set:

            #total_edep += kinematics.total_edep_charged_e(tid,vert_id,seg) # *** total visible energy ***
            contained_edep+= kinematics.fv_edep_e(tid, vert_id, seg)
            
            contained_length+=kinematics.fv_edep_length(tid, vert_id, seg) # *** total visible track length ***
            total_length+=kinematics.total_edep_charged_length(tid, vert_id, seg)
        
        child_total_edep.append(total_edep)
        child_contained_edep.append(contained_edep)
        child_total_length.append(total_length)
        child_contained_length.append(contained_length)
        child_fraction_traj_id_contained.append(loc_class.particle_containment_traj_set_fraction(final_states, track_id_set))
        child_shower_electron_traj_id_fraction.append(particle_assoc.shower_pdg_fraction(track_id_set, final_states, traj, ghdr))

        # Characterize pi0 decay product Endpoints/Containments
        if len(initial_child_track_id_set)>1:
            start_tid_mask = final_states['traj_id']==start_tid
            child_start_traj = final_states[start_tid_mask]
            start_pt = child_start_traj['xyz_start']
            start_pxyz = child_start_traj['pxyz_start']
            if len(start_pt) > 1:
                print("Multiple start points found for trajectory ID:", start_tid)
            else:
                start_pt = start_pt[0]
                start_pxyz = start_pxyz[0]
            
            end_tid = particle_assoc.find_secondary_particle_end_trajectory(track_id_set, final_states,start_tid)
            end_tid_mask = final_states['traj_id']==end_tid
            child_end_traj = final_states[end_tid_mask]
            end_pt = child_end_traj['xyz_end']
            if len(end_pt) > 1:
                print("Multiple end points found for trajectory ID:", end_tid)
            else:
                end_pt = end_pt[0]
        else:
            end_tid = start_tid
            end_pt = fs['xyz_end']
            start_pt = fs['xyz_start']
            start_pxyz = fs['pxyz_start']

        #end_pt_loc = loc_class.particle_end_loc(start_pt, end_pt)
        child_conv_dist.append(particle_assoc.distance_between_points(start_pt, end_pt))
        child_start_pt.append(start_pt)
        child_end_pt.append(end_pt)
        #child_end_pt_loc.append(end_pt_loc)
        child_start_pxyz.append(start_pxyz)

    if len(child_pdg) < 2: 
        #print("Spill ID:", spill_id)
        too_few_by = 2 - len(child_pdg)
        for i in range(too_few_by):
            child_pdg.append(0)
            child_total_edep.append(0.)
            child_contained_edep.append(0.)
            child_total_length.append(0.)
            child_contained_length.append(0.)
            child_start_pt.append([0.,0.,0.])
            child_end_pt.append([0.,0.,0.])
            child_end_pt_loc.append('0')
            child_start_pxyz.append([0.,0.,0.])

    child_total_p = [np.sqrt(pxyz[0]**2 + pxyz[1]**2 + pxyz[2]**2) for pxyz in child_start_pxyz]
    highest_p_index = child_total_p.index(max(child_total_p))
    if child_pdg[highest_p_index] == 0:
        highest_p_child_dir = np.array([0.,0.,0.])
    else:
        highest_p_child_dir = np.array([child_start_pxyz[highest_p_index][i]/child_total_p[highest_p_index] for i in range(3)])

    for c in range(len(child_pdg)):
        if child_pdg[c]==0:
            child_start_dir.append(-1.)
            continue
        elif c == highest_p_index: 
            child_start_dir.append(0.)
            continue
        else:
            child_dir = np.array([child_start_pxyz[c][i]/child_total_p[c] for i in range(3)])
            child_vs_highest_p_child_angle = np.arccos(np.sum(abs(child_dir * highest_p_child_dir)) / \
                                                   (np.sqrt(np.sum(abs(child_dir)**2)) * np.sqrt(np.sum(abs(highest_p_child_dir)**2)))) # angle is returned in radians using a \dot b = |a||b|cos\theta
            child_start_dir.append(child_vs_highest_p_child_angle)

    #print("Excluded Track IDs:", exclude_track_ids)
    #print("Non-Excluded Track IDs:", set(final_states['traj_id']) - exclude_track_ids)
    non_excluded_track_ids = set(final_states['traj_id']) - exclude_track_ids
    for tid in non_excluded_track_ids:
        #print("Non-Excluded PDG:", final_states[final_states['traj_id']==tid]['pdg_id'])
        non_excluded_particle_parent_id = final_states[final_states['traj_id']==tid]['parent_id']
        if np.isin(non_excluded_particle_parent_id, exclude_track_ids):
            print("------------------ SOMETHING IS WRONG ------------------")
        #print("Non-Excluded Parent PDG:", final_states[final_states['traj_id'] == final_states[final_states['traj_id']==tid]['parent_id']]['pdg_id'])
    # Save collected info in input pi0_dict                                                                                                                         
    pi0_dict[(spill_id,vert_id)]=dict(
        filepath = str(sim_file),
        flow_event_id = str(flow_event_id),
        pi0_start_mom = [float(i) for i in pi0_start_mom],
        pi0_end_loc = [float(i) for i in pi0_end_loc],
        pi0_start_loc = [float(i) for i in pi0_start_loc],
        pi0_end_process = [int(i) for i in pi0_end_process],
        pi0_end_subprocess = [int(i) for i in pi0_end_subprocess],
        pi0_ke = float(pi0_ke),
        child_pdg=[int(pdg) for pdg in child_pdg],
        nu_pdg = int(nu_pdg), 
        child_available_energy=[float(total_edep) for total_edep in child_total_edep],
        child_contained_edep=[float(contained_edep) for contained_edep in child_contained_edep],
        child_total_length=[float(total_length) for total_length in child_total_length],
        child_contained_length=[float(contained_length) for contained_length in child_contained_length],
        nu_energy=float(nu_energy),
        q2 = float(q2),
        #end_pt_loc = [str(end_pt_loc) for end_pt_loc in child_end_pt_loc],
        child_start = [[float(i) for i in start_pt] for start_pt in child_start_pt],
        child_end = [[float(i) for i in end_pt] for end_pt in child_end_pt],
        child_start_dir = [float(dir) for dir in child_start_dir],
        child_conv_dist = [float(dist) for dist in child_conv_dist],
        child_fraction_traj_id_contained = [float(frac) for frac in child_fraction_traj_id_contained],
        child_shower_electron_traj_id_fraction = [float(frac) for frac in child_shower_electron_traj_id_fraction],
        nu_int_type=str(nu_int_type))
    return