################################################################################
##                                                                            ##
##    CONTAINS: Script to create JSON dictionaries describing muons and       ##
##              hadrons in signal events, plots describing muons and hadron   ##
##              in signal events, and a text file listing signal event count  ##
##              information. ONLY TRUTH INFORMATION IS USED.                  ##
##              Modified from https://github.com/edhinkle/mesonless_numubarCC/##
##              blob/main/truth_kinematics/file_parsing/signal_kinematics.py  ##
##                                                                            ##
################################################################################

import h5py, glob, argparse
import numpy as np
import sys
sys.path.append('../../')
import truth_studies.util.file_parsing as file_parsing
import truth_studies.util.location_class as loc_class 
import truth_studies.util.truth_ixn_methods as truth
import truth_studies.util.kinematicVariable_methods as kinematics
import truth_studies.util.dictionary_defs as dict_defs
 
def main(sim_dir, input_type, n_files_processed, flow_dir='/global/cfs/cdirs/dune/www/data/2x2/simulation/productions/MiniRun6.1_1E19_RHC/MiniRun6.1_1E19_RHC.flow/FLOW/'):

    test_count = 0
    files_per_section = 100
    total_sections = str(n_files_processed//files_per_section +1)

    ### NOTE: Current POT scaling is based on MiniRun4 larnd file situation
    if int(n_files_processed) < 1024.: 
        scale_factor = (1./(int(n_files_processed)/1024.))*1.05
    else:
        scale_factor = (1./(1023./1024.))*1.05

    pi0_dict = dict() # Initialize muon dictionary
    muon_dict = dict() # Initialize muon dictionary
    hadron_dict = dict() # Initialize hadron dictionary

    # Dictionaries for combining with other background explorations
    signal_dict = dict() # Initialize dictionary for signal muons for full comparison
    sig_bkg_dict = dict() # Initialize dictionary for signal and w.s. bkg muons for full comparison
    
    file_ext = '' ## Changes based on input type
    event_spill_id = '' ## Changes based on input type

    if input_type == 'larnd': 
        file_ext = '.LARNDSIM.hdf5'
        event_spill_id = 'event_id' #new for MiniRun 4
    elif input_type == 'edep':
        file_ext = '.EDEPSIM.hdf5'
        event_spill_id = 'event_id' #new for MiniRun 4

    signal_count = 0 #(144 in first 100 files originally )
    nc_one_pi0_count = 0
    nc_npi0_count = 0
    for sim_file in glob.glob(sim_dir+'/000*000/*'+file_ext): # Loop over simulation files

        if test_count ==int(n_files_processed) : break
        test_count+=1
        print("Looking at file: ", sim_file)
        #if test_count <149: continue
        if sim_file.find('0000912') != -1: 
            print("---------------SKIPPING PROBLEM FILE---------------")
            continue # Skip MiniRun 5 file 0000912 due to bug (for now)
        #if sim_file.find('0000216') != -1: 
        #    print("---------------SKIPPING PROBLEM FILE---------------")
        #    continue # Skip MiniRun 5 file 0000912 due to bug (for now)

        if (test_count % 5 == 0):
            print("Processing file: ", str(test_count), "/", str(n_files_processed))

        try:
            sim_h5 = h5py.File(sim_file,'r')
        except:
            print("Error opening file: ", sim_file)
            continue
        #print(sim_h5.keys(),'\n')

        # Get flow event information
        file_number = sim_file.split("/")[-1].split(".")[-3]
        file_number_int = int(file_number)
        flow_file = glob.glob(flow_dir+'000*000/MiniRun6.1_1E19_RHC.flow.'+file_number+'.FLOW.hdf5')[0]
        f = h5py.File(flow_file, 'r')

        events = f['charge/events/data']
        mc_truth_ixns = f['mc_truth/interactions/data']
        raw_events = f['charge/raw_events/data']
        raw_events_ref = f['charge/raw_events/ref/charge/events/ref']
        raw_events_region = f['charge/raw_events/ref/charge/events/ref_region']
        mc_truth_ixn_ref = f['charge/raw_events/ref/mc_truth/interactions/ref']
        detector_bounds = f['geometry_info'].attrs['module_RO_bounds']

        ### partition file by spill
        unique_spill = np.unique(sim_h5['trajectories'][event_spill_id])
        #print("Unique spill IDs in file:", unique_spill)
        #print("Number of unique spills in file:", len(unique_spill))
        for spill_id in unique_spill:

            ghdr, gstack, traj, vert, seg = file_parsing.get_spill_data(sim_h5, spill_id, input_type)

            ### partition by vertex ID within beam spill
            #print("Number of unique vertices in spill:", len(vert['vertex_id']))
            for v_i in range(len(vert['vertex_id'])):

                vert_pos= [vert['x_vert'][v_i], vert['y_vert'][v_i], vert['z_vert'][v_i]] 
                vert_in_active_LAr = loc_class.fiducialized_vertex( vert_pos ) # Check vertex location relative to FV

                ##### REQUIRE: neutrino vertex in LAr active volume #####
                if vert_in_active_LAr==False: 
                    continue

                vert_id = vert['vertex_id'][v_i]
                nu_mu = truth.signal_nu_pdg(ghdr, vert_id) # nu_mu OR nu_mu_bar
                is_cc = truth.signal_cc(ghdr, vert_id)
                no_charged_mesons = truth.non_pi0_meson_status(gstack, vert_id)
                one_pi0 = truth.single_pi0_status(gstack, vert_id)
                multi_pi0 = truth.multi_pi0_status(gstack, vert_id)
                fv_particle_origin=loc_class.fiducialized_particle_origin(traj, vert_id)

                if is_cc==True: continue

                ### REQUIRE: (A) nu_mu(_bar), (B) CC interaction, (C) NO final state mesons, (D) final state particle start point in FV
                #if nu_mu==True and is_cc==True and no_charged_mesons==True and one_pi0==True and fv_particle_origin==True:
                if one_pi0==True and fv_particle_origin==True:
                    sig_or_bkg = True
                    signal_count+=1
                    nc_one_pi0_count+=1
                elif multi_pi0==True and fv_particle_origin==True:
                    sig_or_bkg = True
                    signal_count+=1
                    nc_npi0_count+=1
                else:
                    continue
                    #sig_or_bkg = False

                #print("Signal or Background: ", sig_or_bkg)

                #print("Sim file: ", sim_file)
                #print("Spill ID: ", spill_id)
                #print("Spill ID index:", np.where(unique_spill==spill_id))
                #print("Length of unique spill array:", np.shape(unique_spill))
                
                # Get flow event id information
                #print("Raw Events:", raw_events)
                #print("Vertex Position:", vert_pos)
                #print("Detector Bounds: ",detector_bounds)
                '''gstack_vert_mask = gstack['vertex_id']==vert_id
                gstack_vert = gstack[gstack_vert_mask] # Particle ID information associated with vertex
                # Get final state particles associated with vertex and their multiplicities
                gstack_vert_fs_mask = gstack_vert['part_status']==1 # Excludes initial state particles
                gstack_vert_fs = gstack_vert[gstack_vert_fs_mask]['part_pdg']
                #print("Genie stack primaries:",gstack_vert_fs)
                truth_ev_id = spill_id
                #print("Truth Event ID: ",truth_ev_id)
                mc_truth_ev_id_idx = np.where(mc_truth_ixns['event_id'] == truth_ev_id)[0]
                #print("MC Truth Event ID Index: ",mc_truth_ev_id_idx)
                mc_truth_ixn_ref_for_ixn = []#np.zeros(len(mc_truth_ev_id_idx), dtype=int)
                raw_event_ref_from_mc_truth = []#np.zeros(len(mc_truth_ev_id_idx), dtype=int)
                for i in range(len(mc_truth_ev_id_idx)):
    
                    ev_id_idx = mc_truth_ev_id_idx[i]
                    mc_truth_ixn_ref_values = mc_truth_ixn_ref[mc_truth_ixn_ref[:,1] == ev_id_idx,  1]
                    if len(mc_truth_ixn_ref_values) > 0:
                        mc_truth_ixn_ref_for_ixn.append(mc_truth_ixn_ref_values[0])
                        raw_event_ref_from_mc_truth.append(mc_truth_ixn_ref[mc_truth_ixn_ref[:,1] == ev_id_idx, 0])
                    else: continue
                #print("Raw event ref from mc_truth before np unique: ",raw_event_ref_from_mc_truth)
                raw_event_ref_from_mc_truth = np.unique((np.concatenate(raw_event_ref_from_mc_truth)))
                #print("Raw event ref from MC Truth: ",raw_event_ref_from_mc_truth)
                mc_truth_ixn_id = mc_truth_ixns[mc_truth_ixn_ref_for_ixn]['event_id']
                #print("MC Truth Interaction ID: ",mc_truth_ixn_id)
                raw_event_id = raw_events[raw_event_ref_from_mc_truth]['id'][0]
                #print("Raw Event ID: ",raw_event_id)

                raw_event_ref = raw_events_ref[raw_events_region[raw_event_id,'start']:raw_events_region[raw_event_id,'stop']]
                raw_event_ref_for_raw_event = np.sort(raw_event_ref[raw_event_ref[:,0] == raw_event_id, 1])
                raw_event_ref_for_event = np.sort(raw_event_ref[raw_event_ref[:,0] == raw_event_id, 0])
                event_id = events[raw_event_ref_for_raw_event]['id'][0]
                #print("Event ID: ",event_id)

                dict_defs.get_truth_sig_bkg_dict(spill_id, vert_id, ghdr, gstack, traj, vert, seg, sig_bkg_dict, sim_file, event_id, file_number_int, sig_or_bkg, is_cc)
                #dict_defs.pi0_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, pi0_dict, sim_file, event_id)
                #dict_defs.muon_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, muon_dict)
                #dict_defs.hadron_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, kinematics.threshold, hadron_dict)
                #dict_defs.get_truth_dict(spill_id, vert_id, ghdr, gstack, traj, vert, seg, signal_dict)'''

        if test_count % files_per_section == 0:
            print("Signal Event Count: ", signal_count)
            print("NC 1pi0 Event Count: ", nc_one_pi0_count)
            print("NC npi0 Event Count: ", nc_npi0_count)
            print("Number of files processed: ", test_count)
            section = str(test_count//files_per_section)
            #file_parsing.save_dict_to_json(sig_bkg_dict, "NC_sig_bkg_dict_"+section+"_of_"+total_sections, True)
            #sig_bkg_dict = dict()

    # Save all Python dictionaries to JSON files
    #file_parsing.save_dict_to_json(signal_dict, "signal_dict", True)
    #file_parsing.save_dict_to_json(pi0_dict, "pi0_dict", True)
    #file_parsing.save_dict_to_json(muon_dict, "muon_dict", True)
    #file_parsing.save_dict_to_json(hadron_dict, "hadron_dict", True)
    #file_parsing.save_dict_to_json(sig_bkg_dict, "NC_sig_bkg_dict_"+total_sections+"_of_"+total_sections, True)

    # Save full signal and w.s. bkg counts to TXT file
    signal_count_final = signal_count*scale_factor
    print("NC 1pi0 Event Count Final [UNSCALED]: ", nc_one_pi0_count)
    print("NC Npi0 Event Count Final [UNSCALED]: ", nc_npi0_count)
    outfile = open('signal_event_counts.txt', "w")
    outfile.writelines(["Signal Events (scaled to 1.05e19 POT): "+str(signal_count_final)+"\n", \
                        "NC 1pi0 Events (scaled to 1.05e19 POT): "+str(nc_one_pi0_count*scale_factor)+"\n", \
                        "NC npi0 Events (scaled to 1.05e19 POT): "+str(nc_npi0_count*scale_factor)+"\n", \
                        "Number of files used to get count: "+str(n_files_processed)+"\n", \
                        "Scale factor:"+str(scale_factor)+"\n"])
    outfile.close()



if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--sim_dir', default=None, required=True, type=str, \
                        help='''string corresponding to the path of the directory containing edep-sim or larnd ouput simulation file to be considered''')
    parser.add_argument('-t', '--input_type', default='larnd', choices=['edep', 'larnd'], type=str, \
                        help='''string corresponding to the output file type: edep or larnd''')
    parser.add_argument('-n', '--n_files_processed', default=1, required=True, type=int, \
                        help='''File count of number of files processed in production sample''')
    parser.add_argument('-f', '--flow_dir', default='/global/cfs/cdirs/dune/www/data/2x2/simulation/productions/MiniRun6.1_1E19_RHC/MiniRun6.1_1E19_RHC.flow/FLOW/', type=str, \
                        help='''string corresponding to the path of the directory containing flow files''')
    args = parser.parse_args()
    main(**vars(args))