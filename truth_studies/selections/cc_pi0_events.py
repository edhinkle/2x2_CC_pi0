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

def main(sim_dir, input_type, n_files_processed):

    test_count = 0

    ### NOTE: Current POT scaling is based on MiniRun4 larnd file situation
    if int(n_files_processed) < 1024.: 
        scale_factor = (1./(int(n_files_processed)/1024.))*1.2
    else:
        scale_factor = 1.2

    pi0_dict = dict() # Initialize muon dictionary
    muon_dict = dict() # Initialize muon dictionary
    hadron_dict = dict() # Initialize hadron dictionary

    # Dictionaries for combining with other background explorations
    signal_dict = dict() # Initialize dictionary for signal muons for full comparison
    
    file_ext = '' ## Changes based on input type
    event_spill_id = '' ## Changes based on input type

    if input_type == 'larnd': 
        file_ext = '.LARNDSIM.hdf5'
        event_spill_id = 'event_id' #new for MiniRun 4
    elif input_type == 'edep':
        file_ext = '.EDEPSIM.hdf5'
        event_spill_id = 'event_id' #new for MiniRun 4

    for sim_file in glob.glob(sim_dir+'/000*000/*'+file_ext): # Loop over simulation files

        if test_count ==int(n_files_processed) : break
        test_count+=1
        #if test_count <369: continue

        if (test_count % 5 == 0):
            print("Processing file: ", str(test_count), "/", str(n_files_processed))

        try:
            sim_h5 = h5py.File(sim_file,'r')
        except:
            print("Error opening file: ", sim_file)
            continue
        #print(sim_h5.keys(),'\n')

        ### partition file by spill
        unique_spill = np.unique(sim_h5['trajectories'][event_spill_id])
        print("Number of unique spills in file:", len(unique_spill))
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
                fv_particle_origin=loc_class.fiducialized_particle_origin(traj, vert_id)

                ### REQUIRE: (A) nu_mu(_bar), (B) CC interaction, (C) NO final state mesons, (D) final state particle start point in FV
                if nu_mu==True and is_cc==True and no_charged_mesons==True and one_pi0==True and fv_particle_origin==True:
                    print("Sim file: ", sim_file)
                    print("Spill ID: ", spill_id)
                    print("Spill ID index:", np.where(unique_spill==spill_id))
                    dict_defs.pi0_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, pi0_dict)
                    dict_defs.muon_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, muon_dict)
                    dict_defs.hadron_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, kinematics.threshold, hadron_dict)
                    dict_defs.get_truth_dict(spill_id, vert_id, ghdr, gstack, traj, vert, seg, signal_dict)

    # Save all Python dictionaries to JSON files
    file_parsing.save_dict_to_json(signal_dict, "signal_dict", True)
    file_parsing.save_dict_to_json(pi0_dict, "pi0_dict", True)
    file_parsing.save_dict_to_json(muon_dict, "muon_dict", True)
    file_parsing.save_dict_to_json(hadron_dict, "hadron_dict", True)

    # Save full signal and w.s. bkg counts to TXT file
    signal_count = len(signal_dict)*scale_factor
    outfile = open('signal_event_counts.txt', "w")
    outfile.writelines(["Signal Events (scaled to 1.2e19 POT): "+str(signal_count)+"\n", \
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
    args = parser.parse_args()
    main(**vars(args))