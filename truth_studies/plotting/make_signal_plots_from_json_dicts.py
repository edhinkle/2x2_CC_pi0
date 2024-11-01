################################################################################
##                                                                            ##
##    CONTAINS: Script to create plots describing muons, hadrons, pi0s in     ##
##              signal events, using  dictionaries created using the          ##
##              methods in /truth_studies/selections/dictionary_defs.py       ##
##              and a scale factor for scaling event counts to those expected ##
##              with 1.2e19 POT. Adapted from https://github.com/edhinkle/    ##
##              mesonless_numubarCC/blob/main/truth_kinematics/plotting/      ##
##              signal_plots_from_json_dicts.py                               ##
##                                                                            ##
################################################################################

import json
import argparse
import sys
sys.path.append('../../')
from truth_studies.plotting.plot_pi0_in_single_pi0_events import plot_pi0s
from truth_studies.plotting.plot_muons_in_single_pi0_events import plot_muons
from truth_studies.plotting.plot_hadrons_in_single_pi0_events import plot_hadrons

def main(scale_factor, pi0_json_file, muon_json_file, hadron_json_file):

    pi0_file = open(pi0_json_file)
    pi0_dict=json.load(pi0_file)

    muon_file = open(muon_json_file)
    muon_dict=json.load(muon_file)

    hadron_file = open(hadron_json_file)
    hadron_dict=json.load(hadron_file)

    plot_pi0s(pi0_dict, scale_factor, sig_bkg = 0)
    plot_muons(muon_dict, scale_factor, sig_bkg = 0)
    plot_hadrons(hadron_dict, scale_factor, sig_bkg = 0)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-sf', '--scale_factor', default=1, required=True, type=float, \
                        help='''Scale factor related to input dictionary''')
    parser.add_argument('-npi', '--pi0_json_file', default=None, required=True, type=str, \
                        help='''string corresponding to the path of the pi0 info JSON file''')
    parser.add_argument('-mu', '--muon_json_file', default=None, required=True, type=str, \
                        help='''string corresponding to the path of the muon info JSON file''')
    parser.add_argument('-had', '--hadron_json_file', default=None, required=True, type=str, \
                        help='''string corresponding to the path of the hadron info JSON file''')
    args = parser.parse_args()
    main(**vars(args))