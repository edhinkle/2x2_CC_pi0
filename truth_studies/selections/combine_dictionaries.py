import h5py, glob, argparse
import numpy as np
import sys
import pandas as pd
import json 
sys.path.append('../../')
import truth_studies.util.file_parsing as file_parsing

def main(directory):


    sig_bkg_dict_FINAL = {}

    for i, dict_file in enumerate(glob.glob(directory+'/sig_bkg_dict_*_of_11.json')):
        with open(dict_file) as sig_bkg_file:
            sig_bkg_dict = json.load(sig_bkg_file)
            sig_bkg_dict_FINAL.update(sig_bkg_dict)

    # Save the combined dictionary to a new JSON file
    with open(directory + '/combined_sig_bkg_dict.json', 'w') as outfile:
        json.dump(sig_bkg_dict_FINAL, outfile, indent=4)



    #file_parsing.save_dict_to_json(sig_bkg_dict, directory+'/sig_bkg_dict_FINAL', True)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--directory', default=None, required=True, type=str, \
                        help='''string corresponding to the path of the directory containing edep-sim or larnd ouput simulation file to be considered''')
    args = parser.parse_args()
    main(**vars(args))