/* #include "/cvmfs/dune.opensciencegrid.org/products/dune/srproxy/v00.43/include/SRProxy/BasicTypesProxy.h" */
#include "duneanaobj/StandardRecord/StandardRecord.h"
#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "TEfficiency.h"
#include <iostream>
#include <fstream>
#include <string>
#include <nlohmann/json.hpp>
using json = nlohmann::json;
//Requires a file containing a list of input CAF files and returns an int code for success/error
//The list should be one file/path per line, and files/lines can be commented out using #
//Boolean flag to switch behavior for reading structured/flat CAFs (defaults to structured CAFs)
int general_CC1pi0_selection(const std::string& file_list, const std::string& json_file_path, bool is_mc = true, bool is_flat = true)
{
    // Reads list of ROOT files in external file and writes to new list
    std::vector<std::string> root_list;
    std::ifstream fin(file_list, std::ios::in);
    if(!fin.is_open())
    {
        std::cerr << "Failed to open " << file_list << std::endl;
        std::cerr << "Exiting" << std::endl;
        return 111;
    }
    else
    {
        std::cout << "Reading " << file_list << " for input ROOT files." << std::endl;
        std::string name;
        while(std::getline(fin, name))
        {
            if(name.front() == '#')
                continue;
            root_list.push_back(name);
        }
    }

    if(root_list.empty())
    {
        std::cerr << "No input ROOT files. Exiting." << std::endl;
        return 121;
    }

    // Adds files in ROOT list to TChain
    TChain* caf_chain = new TChain("cafTree");
    for(const auto& file : root_list)
    {
        std::cout << "Adding " << file << " to TChain." << std::endl;
        caf_chain->Add(file.c_str());
    }
    std::cout << "Finished adding files..." << std::endl;

    // Load JSON dictionary of true events
    std::ifstream json_file(json_file_path); // replace with your JSON file path
    if (!json_file.is_open()) {
        std::cerr << "Failed to open JSON file." << std::endl;
        return 131;
    }

    json json_data;
    json_file >> json_data;

    // DEFINE: Vectors to hold information to keep in output TTree file (Reco Tree)
    std::vector< double >  reco_lead_shower_energy;
    std::vector< double >  reco_lead_shower_p_x; 
    std::vector< double >  reco_lead_shower_p_y; 
    std::vector< double >  reco_lead_shower_p_z;
    std::vector< double >  reco_sublead_shower_energy;
    std::vector< double >  reco_sublead_shower_p_x; 
    std::vector< double >  reco_sublead_shower_p_y; 
    std::vector< double >  reco_sublead_shower_p_z;
    std::vector< double >  reco_lead_shower_start_x;
    std::vector< double >  reco_lead_shower_start_y;
    std::vector< double >  reco_lead_shower_start_z;
    std::vector< double >  reco_sublead_shower_start_x;
    std::vector< double >  reco_sublead_shower_start_y;
    std::vector< double >  reco_sublead_shower_start_z;
    std::vector< double >  reco_muon_angle;
    std::vector< double >  reco_muon_angle_Mx2_match;
    std::vector< int >     reco_lead_shower_pdg;
    std::vector< int >     reco_sublead_shower_pdg;
    std::vector< int >     reco_ixn_gamma_mult;
    std::vector< int >     reco_ixn_e_mult;
    std::vector< int >     reco_ixn_cont_gamma_mult;
    std::vector< int >     reco_ixn_e_cont_mult;
    std::vector< int >     reco_ixn_muon_mult;
    std::vector< int >     reco_ixn_chpi_mult;
    std::vector< int >     reco_ixn_proton_mult;
    std::vector< int >     reco_ixn_chkaon_mult;
    std::vector< double >     reco_ixn_vtx_x_pos;
    std::vector< double >     reco_ixn_vtx_y_pos;
    std::vector< double >     reco_ixn_vtx_z_pos;

    std::vector< double >  true_match_lead_shower_energy;
    std::vector< double >  true_match_lead_shower_p_x; 
    std::vector< double >  true_match_lead_shower_p_y; 
    std::vector< double >  true_match_lead_shower_p_z;
    std::vector< double >  true_match_sublead_shower_energy;
    std::vector< double >  true_match_sublead_shower_p_x; 
    std::vector< double >  true_match_sublead_shower_p_y; 
    std::vector< double >  true_match_sublead_shower_p_z;
    std::vector< double >  true_match_lead_shower_start_x;
    std::vector< double >  true_match_lead_shower_start_y;
    std::vector< double >  true_match_lead_shower_start_z;
    std::vector< double >  true_match_sublead_shower_start_x;
    std::vector< double >  true_match_sublead_shower_start_y;
    std::vector< double >  true_match_sublead_shower_start_z;
    std::vector< double >  true_match_muon_angle;
    std::vector< int >     true_match_lead_shower_pdg;
    std::vector< int >     true_match_sublead_shower_pdg;
    std::vector< double >  true_match_lead_ovlp;
    std::vector< double >  true_match_sublead_ovlp;
    std::vector< int >     true_ixn_pi0_mult;
    std::vector< int >     true_ixn_e_mult;
    std::vector< int >     true_ixn_gamma_mult;
    std::vector< int >     true_ixn_cont_pi0_mult;
    std::vector< int >     true_ixn_muon_mult;
    std::vector< int >     true_ixn_chpi_mult;
    std::vector< int >     true_ixn_proton_mult;
    std::vector< int >     true_ixn_chkaon_mult;
    std::vector< double >     true_ixn_vtx_x_pos;
    std::vector< double >     true_ixn_vtx_y_pos;
    std::vector< double >     true_ixn_vtx_z_pos;
    std::vector< bool > true_ixn_is_cc;
    std::vector< bool > true_ixn_is_in_fv;
    std::vector< double > true_ixn_target;
    std::vector< std::string > true_ixn_mode; 
    std::vector< int > true_ixn_nu_pdg;
    std::vector< bool >  true_ixn_is_signal;

    std::vector< double >  overlap;
    std::vector< double >  true_ixn_index;
    std::vector< double >  reco_ixn_index;
    std::vector< int >     spill_index;
    std::vector< int >     file_index;
    std::vector< int >     event;
    std::vector< int >     run;
    std::vector< int >     subrun;
    std::vector< bool >  in_truth_dict;
    std::vector< std::string > caf_file_name;
    std::vector< std::string > truth_dict_key;
    
    // DEFINE: Vectors to hold information to keep in output TTree file (Truth Tree)
    std::vector< double > truth_vtx_x;
    std::vector< double > truth_vtx_y;
    std::vector< double > truth_vtx_z;
    std::vector< double > truth_muon_energy;
    std::vector< double > truth_muon_angle;
    std::vector< double > truth_nPionsPrimary;
    std::vector< double > truth_nProtonsPrimary;
    std::vector< double > truth_nKaonsPrimary;
    std::vector< double > truth_nMuonsPrimary;
    std::vector< double > truth_nPi0sPrimary;
    std::vector< double > truth_nElectronsPrimary;
    std::vector< double > truth_nPhotonsPrimary;
    std::vector< double > truth_found_file;
    std::vector< double > truth_file_index;

    // DEFINE: Vectors to hold information to keep in output TTree file (Purity/Efficiency Tree)
    std::vector< int > total_ixns;
    std::vector< int > total_signal;
    std::vector< int > signal_reco_no_cuts;
    std::vector< int > all_reco_no_cuts;
    std::vector< int > signal_post_vtx_cut;
    std::vector< int > all_post_vtx_cut;
    std::vector< int > signal_post_muon_cut;
    std::vector< int > all_post_muon_cut;
    std::vector< int > signal_post_mx2_cut;
    std::vector< int > all_post_mx2_cut;
    std::vector< int > signal_post_shower_cut;
    std::vector< int > all_post_shower_cut;
    std::vector< int > signal_post_pion_cut;
    std::vector< int > all_post_pion_cut;
    std::vector< int > signal_post_kaon_cut;
    std::vector< int > all_post_kaon_cut;

    // DEFINE: TTree and TBranches to go in output ROOT file

    // Define branches for Purity/Efficiency Info Tree
    TTree *fPurEffTree=new TTree("PurEffTree", "Purity and Efficiency Counting");
    fPurEffTree->Branch("total_ixns", &total_ixns);
    fPurEffTree->Branch("total_signal", &total_signal);
    fPurEffTree->Branch("signal_reco_no_cuts", &signal_reco_no_cuts);
    fPurEffTree->Branch("all_reco_no_cuts", &all_reco_no_cuts);
    fPurEffTree->Branch("signal_post_vtx_cut", &signal_post_vtx_cut);
    fPurEffTree->Branch("all_post_vtx_cut", &all_post_vtx_cut);
    fPurEffTree->Branch("signal_post_muon_cut", &signal_post_muon_cut);
    fPurEffTree->Branch("all_post_muon_cut", &all_post_muon_cut);
    fPurEffTree->Branch("signal_post_mx2_cut", &signal_post_mx2_cut);
    fPurEffTree->Branch("all_post_mx2_cut", &all_post_mx2_cut);
    fPurEffTree->Branch("signal_post_shower_cut", &signal_post_shower_cut);
    fPurEffTree->Branch("all_post_shower_cut", &all_post_shower_cut);
    fPurEffTree->Branch("signal_post_pion_cut", &signal_post_pion_cut);
    fPurEffTree->Branch("all_post_pion_cut", &all_post_pion_cut);
    fPurEffTree->Branch("signal_post_kaon_cut", &signal_post_kaon_cut);
    fPurEffTree->Branch("all_post_kaon_cut", &all_post_kaon_cut);

    // Define branches for Truth Info Tree
    TTree *fTruthTree=new TTree("TruthTree", "Truth Variables");
    fTruthTree->Branch("truth_vtx_x", &truth_vtx_x);
    fTruthTree->Branch("truth_vtx_y", &truth_vtx_y);
    fTruthTree->Branch("truth_vtx_z", &truth_vtx_z);
    fTruthTree->Branch("truth_muon_energy", &truth_muon_energy);
    fTruthTree->Branch("truth_muon_angle", &truth_muon_angle);
    fTruthTree->Branch("truth_nPionsPrimary", &truth_nPionsPrimary);
    fTruthTree->Branch("truth_nProtonsPrimary", &truth_nProtonsPrimary);
    fTruthTree->Branch("truth_nKaonsPrimary", &truth_nKaonsPrimary);
    fTruthTree->Branch("truth_nMuonsPrimary", &truth_nMuonsPrimary);
    fTruthTree->Branch("truth_nPi0sPrimary", &truth_nPi0sPrimary);
    fTruthTree->Branch("truth_nElectronsPrimary", &truth_nElectronsPrimary);
    fTruthTree->Branch("truth_nPhotonsPrimary", &truth_nPhotonsPrimary);
    fTruthTree->Branch("truth_found_file", &truth_found_file);
    fTruthTree->Branch("truth_file_index", &truth_file_index);

    // Define branches for Reco Info Tree
    TTree *fRecoTree=new TTree("RecoTree", "Reco Variables");
    fRecoTree->Branch("reco_lead_shower_energy", &reco_lead_shower_energy);
    fRecoTree->Branch("reco_lead_shower_p_x", &reco_lead_shower_p_x);
    fRecoTree->Branch("reco_lead_shower_p_y", &reco_lead_shower_p_y);
    fRecoTree->Branch("reco_lead_shower_p_z", &reco_lead_shower_p_z);
    fRecoTree->Branch("reco_sublead_shower_energy", &reco_sublead_shower_energy);
    fRecoTree->Branch("reco_sublead_shower_p_x", &reco_sublead_shower_p_x);
    fRecoTree->Branch("reco_sublead_shower_p_y", &reco_sublead_shower_p_y);
    fRecoTree->Branch("reco_sublead_shower_p_z", &reco_sublead_shower_p_z);
    fRecoTree->Branch("reco_lead_shower_start_x", &reco_lead_shower_start_x);
    fRecoTree->Branch("reco_lead_shower_start_y", &reco_lead_shower_start_y);
    fRecoTree->Branch("reco_lead_shower_start_z", &reco_lead_shower_start_z);
    fRecoTree->Branch("reco_sublead_shower_start_x", &reco_sublead_shower_start_x);
    fRecoTree->Branch("reco_sublead_shower_start_y", &reco_sublead_shower_start_y);
    fRecoTree->Branch("reco_sublead_shower_start_z", &reco_sublead_shower_start_z);
    fRecoTree->Branch("reco_muon_angle", &reco_muon_angle);
    fRecoTree->Branch("reco_muon_angle_Mx2_match", &reco_muon_angle_Mx2_match);
    fRecoTree->Branch("reco_lead_shower_pdg", &reco_lead_shower_pdg);
    fRecoTree->Branch("reco_sublead_shower_pdg", &reco_sublead_shower_pdg);
    fRecoTree->Branch("reco_ixn_gamma_mult", &reco_ixn_gamma_mult);
    fRecoTree->Branch("reco_ixn_e_mult", &reco_ixn_e_mult);
    fRecoTree->Branch("reco_ixn_index", &reco_ixn_index);
    fRecoTree->Branch("reco_ixn_cont_gamma_mult", &reco_ixn_cont_gamma_mult);
    fRecoTree->Branch("reco_ixn_e_cont_mult", &reco_ixn_e_cont_mult);
    fRecoTree->Branch("reco_ixn_muon_mult", &reco_ixn_muon_mult);
    fRecoTree->Branch("reco_ixn_chpi_mult", &reco_ixn_chpi_mult);
    fRecoTree->Branch("reco_ixn_proton_mult", &reco_ixn_proton_mult);
    fRecoTree->Branch("reco_ixn_chkaon_mult", &reco_ixn_chkaon_mult);
    fRecoTree->Branch("reco_ixn_vtx_x_pos", &reco_ixn_vtx_x_pos);
    fRecoTree->Branch("reco_ixn_vtx_y_pos", &reco_ixn_vtx_y_pos);
    fRecoTree->Branch("reco_ixn_vtx_z_pos", &reco_ixn_vtx_z_pos);


    fRecoTree->Branch("true_match_lead_shower_energy", &true_match_lead_shower_energy);
    fRecoTree->Branch("true_match_lead_shower_p_x", &true_match_lead_shower_p_x);
    fRecoTree->Branch("true_match_lead_shower_p_y", &true_match_lead_shower_p_y);
    fRecoTree->Branch("true_match_lead_shower_p_z", &true_match_lead_shower_p_z);
    fRecoTree->Branch("true_match_sublead_shower_energy", &true_match_sublead_shower_energy);
    fRecoTree->Branch("true_match_sublead_shower_p_x", &true_match_sublead_shower_p_x);
    fRecoTree->Branch("true_match_sublead_shower_p_y", &true_match_sublead_shower_p_y);
    fRecoTree->Branch("true_match_sublead_shower_p_z", &true_match_sublead_shower_p_z);
    fRecoTree->Branch("true_match_lead_shower_start_x", &true_match_lead_shower_start_x);
    fRecoTree->Branch("true_match_lead_shower_start_y", &true_match_lead_shower_start_y);
    fRecoTree->Branch("true_match_lead_shower_start_z", &true_match_lead_shower_start_z);
    fRecoTree->Branch("true_match_sublead_shower_start_x", &true_match_sublead_shower_start_x);
    fRecoTree->Branch("true_match_sublead_shower_start_y", &true_match_sublead_shower_start_y);
    fRecoTree->Branch("true_match_sublead_shower_start_z", &true_match_sublead_shower_start_z);
    fRecoTree->Branch("true_match_muon_angle", &true_match_muon_angle);
    fRecoTree->Branch("true_match_lead_shower_pdg", &true_match_lead_shower_pdg);
    fRecoTree->Branch("true_match_sublead_shower_pdg", &true_match_sublead_shower_pdg);
    fRecoTree->Branch("true_match_lead_ovlp", &true_match_lead_ovlp);
    fRecoTree->Branch("true_match_sublead_ovlp", &true_match_sublead_ovlp);
    fRecoTree->Branch("true_ixn_pi0_mult", &true_ixn_pi0_mult);
    fRecoTree->Branch("true_ixn_cont_pi0_mult", &true_ixn_cont_pi0_mult);
    fRecoTree->Branch("true_ixn_e_mult", &true_ixn_e_mult);
    fRecoTree->Branch("true_ixn_gamma_mult", &true_ixn_gamma_mult);
    fRecoTree->Branch("true_ixn_index", &true_ixn_index);
    fRecoTree->Branch("true_ixn_muon_mult", &true_ixn_muon_mult);
    fRecoTree->Branch("true_ixn_chpi_mult", &true_ixn_chpi_mult);
    fRecoTree->Branch("true_ixn_proton_mult", &true_ixn_proton_mult);
    fRecoTree->Branch("true_ixn_chkaon_mult", &true_ixn_chkaon_mult);
    fRecoTree->Branch("true_ixn_vtx_x_pos", &true_ixn_vtx_x_pos);
    fRecoTree->Branch("true_ixn_vtx_y_pos", &true_ixn_vtx_y_pos);
    fRecoTree->Branch("true_ixn_vtx_z_pos", &true_ixn_vtx_z_pos);
    fRecoTree->Branch("true_ixn_is_cc", &true_ixn_is_cc);
    fRecoTree->Branch("true_ixn_is_in_fv", &true_ixn_is_in_fv);
    fRecoTree->Branch("true_ixn_target", &true_ixn_target);
    fRecoTree->Branch("true_ixn_mode", &true_ixn_mode);
    fRecoTree->Branch("true_ixn_nu_pdg", &true_ixn_nu_pdg);
    fRecoTree->Branch("true_ixn_is_signal", &true_ixn_is_signal);

    fRecoTree->Branch("spill_index", &spill_index);
    fRecoTree->Branch("file_index", &file_index);
    fRecoTree->Branch("event", &event);
    fRecoTree->Branch("run", &run);
    fRecoTree->Branch("subrun", &subrun);
    fRecoTree->Branch("caf_file_name", &caf_file_name);
    fRecoTree->Branch("in_truth_dict", &in_truth_dict);
    fRecoTree->Branch("truth_dict_key", &truth_dict_key);

    fRecoTree->Branch("overlap", &overlap);

    // DEFINE: Constants for the detector geometry (2x2)
    //Beam direction -3.343 degrees in y
    const auto beam_dir = TVector3(0, -0.05836, 1.0);

    //z-direction (roughly beam dir)
    const auto z_plus_dir = TVector3(0, 0, 1.0);
    const auto y_plus_dir = TVector3(0, 1.0, 0.0);
    const auto x_plus_dir = TVector3(1.0, 0, 0.0);

    //negative y-direction 
    const auto y_minus_dir = TVector3(0, -1.0, 0.0);

    //Center of the 2x2 LAr
    //For MR4.5 CAFs the 2x2 coordinates are actually centered somewhere else
    const float tpc_x_center = 0.0;
    const float tpc_y_center = 0.0; // -268.0;
    const float tpc_z_center = 0.0; //(1333.5 + 1266.5) / 2.0;

    // 2x2 LAr Boundaries for data, MR6 (TO DO: update so as not to hard code)
    // TO-DO: Add FV cuts? 
    // from: https://github.com/DUNE/ndlar_flow/blob/dce6e26a4be8837907e8a378e7becb2573cc416a/event_display/interactive_event_display/display_utils.py#L370C4-L372C73
    const float det_x_min = -63.931;
    const float mod23_x_max = -3.069;
    const float mod12_x_min = 3.069;
    const float det_x_max = 63.931;

    const float det_y_min = -61.8543;
    const float det_y_max = 61.8543;

    const float det_z_min = -64.3163;
    const float ups_z_max = -2.6837;
    const float downs_z_min = 2.6837;
    const float det_z_max = 64.3163;

    const float max_drift_distance = 30.27225; // cm
    const float cathode_thickness = ((det_x_max - mod12_x_min) - (2.0 * max_drift_distance));
    const float cathode_abs_x = det_x_max - max_drift_distance - (cathode_thickness / 2.0); // cm

    // Mx2 Offsets
    double mnvOffsetX=-10; 
    double mnvOffsetY=5;
    double mnvOffsetZDS=0; //cm
    double mnvDistanceAllowance=15; //cm 
    double mnvAngleAllowance=0.06; //rad
    if (is_mc){ mnvOffsetX=0; mnvOffsetY=0; mnvOffsetZDS=0;}
    double mnvMatchDotProdCut=0.99; //cosine of angle between two vectors

    // Muon energy and angle cuts
    double muon_energy_cut = 0.0; 
    double muon_cos_angle_cut = 0.0;

    // Track length cut
    double minTrkLength = 0.0; // cm

    // Vtx allowance true vs. reco
    double vtx_allowance = 5.0; // cm

    // Fiducial volume cut
    double fv_cut_x = 2.0; // cm
    double fv_cut_y = 2.0; // cm
    double fv_cut_z = 3.0; // cm
    double fv_cut_cathode_x = 2.0; // cm

    // Module geometry limits after FV cuts
    const float fv_abs_x_min = abs(mod12_x_min) + fv_cut_x;
    const float fv_abs_x_max = abs(det_x_max) - fv_cut_x;

    const float fv_abs_y_max = abs(det_y_max) - fv_cut_y;

    const float fv_abs_z_min = abs(downs_z_min) + fv_cut_z;
    const float fv_abs_z_max = abs(det_z_max) - fv_cut_z;

    const float fv_abs_cathode_x_high = cathode_abs_x + (cathode_thickness / 2.0) + fv_cut_cathode_x;
    const float fv_abs_cathode_x_low = cathode_abs_x - (cathode_thickness / 2.0) - fv_cut_cathode_x;



    // Initialize variables for purity/efficiency counting
    int trueInteractions = 0.; int ixnsNoCuts = 0.; int ixnsVtxCut = 0.; int ixnsMuonCut = 0.; int ixnsMx2Cut = 0.; int ixnsShowerCut = 0.; int ixnsPionCut = 0.; int ixnsKaonCut = 0.;
    int trueSignal = 0.; int signalNoCuts = 0.; int signalVtxCut = 0.; int signalMuonCut = 0.; int signalMx2Cut = 0.; int signalShowerCut = 0.; int signalPionCut = 0.; int signalKaonCut = 0.;
   

    
    //Attach SR object to the tree, not using SRProxy (SR Proxy used below)
    //auto sr = new caf::StandardRecord;
    //caf_chain->SetBranchAddress("rec", &sr);

    //const unsigned long nspills = caf_chain->GetEntries();
    //const unsigned int incr = nspills / 10;
    //std::cout << "Looping over " << nspills << " entries/spills..." << std::endl;

    //Most outer loop over files for SRProxy
    const auto t_start{std::chrono::steady_clock::now()};
    auto file_num = 0;
    for(const auto& f : root_list)
    {
    std::string current_file = f;
    std::cout << "Processing " << current_file << std::endl;
    file_num++;
    
    // Get information about file to check with JSON dictionary
    // Get file number 
    bool found_file = false;
    std::vector <int> event_ids;
    std::vector <double> true_dict_vtx_x; 
    std::vector <double> true_dict_vtx_y;
    std::vector <double> true_dict_vtx_z;
    std::vector <string> truth_dict_key_strings;
    if(is_mc){
        std::string prefix = "RHC.caf_v2.";
        std::string suffix = ".CAF.flat.root";
        size_t start_pos = current_file.find(prefix);
        if (start_pos != std::string::npos) {
            start_pos += prefix.length();
            size_t end_pos = current_file.find(suffix, start_pos);
            if (end_pos != std::string::npos) {
                std::string file_num_str = current_file.substr(start_pos, end_pos - start_pos);
                file_num = std::stoi(file_num_str);
                //std::cout << "File number: " << file_num << std::endl;
            } else {
                std::cerr << "ERROR: could not find " << suffix << " in file name " << current_file << std::endl;
                return 1;
            }
        } else {
            std::cerr << "ERROR: could not find " << prefix << " in file name " << current_file << std::endl;
            return 1;
        }  
    
        // Check if this file shows up in the JSON file
        // Filter JSON data by "file_number"
        event_ids.clear();
        true_dict_vtx_x.clear();
        true_dict_vtx_y.clear();
        true_dict_vtx_z.clear();
        for (auto& it : json_data.items()) {
            auto item = it.value();
            auto it_key = it.key();
            auto spill_id = 0; // Initialize spill_id
            if (item.contains("file_number") && item["file_number"].is_number_integer()) {
                int json_file_number = item["file_number"];
                if (json_file_number == file_num) {
                    std::cout << "Found matching file_number: " << json_file_number << std::endl;
                    found_file = true;
                    std::string s = it_key; // Convert to string
                    size_t dash_pos = s.find('-');
                    if (dash_pos != std::string::npos) {
                        spill_id = std::stoi(s.substr(0, dash_pos));
                        // spill_id now holds 1118
                    }
                    event_ids.push_back(spill_id);
                    std::cout << "Spill ID: " << spill_id << std::endl;
                    true_dict_vtx_x.push_back(item["vtx_x"]);
                    true_dict_vtx_y.push_back(item["vtx_y"]);
                    true_dict_vtx_z.push_back(item["vtx_z"]);
                    truth_dict_key_strings.push_back(it_key);
                }
            }
        }
    }
    //if (!found_file) {
    //    //std::cout << "No events from File " << file_num << " in JSON dictionary." << std::endl;
    //    continue;
    //}
    
    //Open file and attach SRProxy object
    //Different tree name for structured and flat CAFs
    TFile* caf_file = TFile::Open(f.c_str(), "READ");
    TTree* caf_tree = (TTree*)caf_file->Get("cafTree");
    std::string tree_name = is_flat ? "rec" : "";
    auto sr = new caf::SRProxy(caf_tree, tree_name);

    // Get number of spills for looping
    const unsigned long nspills = caf_tree->GetEntries();
    const unsigned int incr = nspills / 10;

    // If MC, loop over each spill to get truth level information independent of reco 
    if(is_mc) {
        bool skipEvent=false;
        for(unsigned long n = 0; n < nspills; ++n){ 
            double trueVtxX=-9999; double trueVtxY=-999; double trueVtxZ=-9999;
            skipEvent=false;
            //if(n%10 == 0) std::cout << Form("Truth Level Event %ld of %ld", n, nspills) << std::endl;
            caf_tree->GetEntry(n); //Get spill from tree
            bool hasANeutrino=false;
            const auto nnu = sr->mc.nnu; 
            const auto nnu_value = static_cast<unsigned long>(nnu);
            //std::cout << "DEBUG: Number of neutrinos: " << nnu_value << std::endl;
            for(unsigned long ntrue=0; ntrue<nnu_value; ++ntrue){ 

                //std::cout << Form("DEBUG: In true vertex loop %ld of %ld", ntrue, nnu_value) << std::endl;
                const auto* truth_ixn = &sr->mc.nu[ntrue];
                const auto* vertex=&truth_ixn->vtx;
                //std::cout << "DEBUG: Got vertex" << std::endl;
                trueVtxX=vertex->x; trueVtxY=vertex->y; trueVtxZ=vertex->z;
                const auto* truePrimary=&truth_ixn->prim;
                //std::cout << "DEBUG: Got primaries" << std::endl;
                //std::cout << "DEBUG: Vertex X: " << trueVtxX << std::endl;
                //std::cout << "DEBUG: Target PDG: " << truth_ixn->targetPDG << std::endl;
                int truePart=0; int truePartTrkOnly=0; int truePartNoG4=0;
                if (truth_ixn->targetPDG!=1000180400) continue; // keep if target is argon
                //std::cout << "DEBUG: Past target cut" << std::endl;
                trueInteractions++;
                if (truth_ixn->iscc==false) continue; // keep if is CC event
                //std::cout << "DEBUG: Past CC cut" << std::endl;
                if (abs(truth_ixn->pdg)!=14) continue; // keep if neutrino is muon neutrino
                //std::cout << "DEBUG: Past PDG cut" << std::endl;
                // Make sure vertex is in the detector
                if (abs(trueVtxX)>fv_abs_x_max || abs(trueVtxX)<fv_abs_x_min) continue;
                //std::cout << "DEBUG: Past X cut" << std::endl;
                if (abs(trueVtxZ)>fv_abs_z_max || abs(trueVtxZ)<fv_abs_z_min) continue;
                //std::cout << "DEBUG: Past Y cut" << std::endl;
                if (abs(trueVtxY)>fv_abs_y_max) continue;    
                if (abs(trueVtxX) > fv_abs_cathode_x_low && abs(trueVtxX) < fv_abs_cathode_x_high) continue;
                //std::cout << "DEBUG: Past Z cut" << std::endl;

                // Particle multiplicity checks
                int npip=truth_ixn->npip;
                int npim=truth_ixn->npim;
                int np=truth_ixn->nproton;
                int nneutron=truth_ixn->nneutron;
                int npi0=truth_ixn->npi0;
                int npipPrimaries=0;
                int npimPrimaries=0;
                int npPrimaries=0;
                int nneutronPrimaries=0;
                int npi0Primaries=0;
                //std::cout << "DEBUG: About to enter true primaries loop #1" << std::endl;
                for (int k=0; k<sr->mc.nu[ntrue].prim.size(); k++){ 
                    //std::cout << "DEBUG: In true primaries loop #1" << std::endl;
                    if(sr->mc.nu[ntrue].prim[k].pdg==211) npipPrimaries++; 
                    if(sr->mc.nu[ntrue].prim[k].pdg==2212) npPrimaries++;  
                    if(sr->mc.nu[ntrue].prim[k].pdg==-211) npimPrimaries++;  
                    if(sr->mc.nu[ntrue].prim[k].pdg==2112) nneutronPrimaries++;
                    if(sr->mc.nu[ntrue].prim[k].pdg==111) npi0Primaries++;
                }
                if (npi0!=1) continue; // keep if there is one pi0
                //if (npipPrimaries>0) continue; // keep if there are no primary pi+ in the event
                //if (npimPrimaries>0) continue; // keep if there are no primary pi- in the event
            
                int npipSecondaries = npipPrimaries-npip;
                int npimSecondaries = npimPrimaries-npim;
                int npSecondaries = npPrimaries-np;
                int nneutronSecondaries = nneutronPrimaries-nneutron;
                int npi0Secondaries = npi0Primaries-npi0;

                int nProton=0; int nPion=0; int nPi0=0; int escapingPions=0; int nKaon=0;
                int nPhoton=0; int nElectron=0; int nMuon=0;
                double  Elep=-999; double cosL=-999; double muon_angle=-999;
                for(long unsigned primaries=0; primaries<truePrimary->size(); primaries++){
                    //std::cout << "DEBUG: In true primaries loop #2" << std::endl;
                    int pdg=(*truePrimary)[primaries].pdg;
                    double E=(*truePrimary)[primaries].p.E;
                    double px=(*truePrimary)[primaries].p.px; double py=(*truePrimary)[primaries].p.py; double pz=(*truePrimary)[primaries].p.pz;
                    double totP=TMath::Sqrt(px*px+py*py+pz*pz);

                    if (pdg==22) nPhoton++;
                    if (abs(pdg)==11) nElectron++;
                    if (pdg==111) nPi0++;
                
                    if ((abs(pdg)==13 || abs(pdg)==2212 || abs(pdg)==211 || abs(pdg)==321)){
                        //std::cout << "DEBUG: Look at track-like true particles" << std::endl;
                        truePartNoG4++;
                        if (abs(pdg)==211)  nPion++; 
                        if (abs(pdg)==2212)  nProton++; 
                        if (abs(pdg)==321)  nKaon++;
                        if (abs(pdg)==13)  nMuon++;
                        auto& start_pos=truth_ixn->prim[primaries].start_pos;
                        auto& end_pos=truth_ixn->prim[primaries].end_pos;
                        //if (std::isnan(start_pos.z)){  continue;}
                        auto& p=truth_ixn->prim[primaries].p; 
                        double dX=abs(end_pos.x-start_pos.x);
                        double dY=abs(end_pos.y-start_pos.y);
                        double dZ=abs(end_pos.z-start_pos.z);
                        double dZ_signed = end_pos.z-start_pos.z;
                        double length=TMath::Sqrt(dX*dX+dY*dY+dZ*dZ);
                        double lengthP=TMath::Sqrt(p.px*p.px+p.py*p.py+p.pz*p.pz);
                        if (truth_ixn->iscc && abs(truth_ixn->pdg)==14 && abs(pdg)==13){
                        
                            //std::cout << "DEBUG: Look at muon" << std::endl;
                            Elep=sr->mc.nu[ntrue].prim[primaries].p.E;
                            cosL=dZ/length;
                            muon_angle = TMath::ACos(dZ_signed/length);
                            //std::cout<<n<<","<<ntrue<<","<<primaries<<","<<dX<<","<<dY<<","<<dZ<<","<<p.pz/lengthP<<std::endl;

                            //trueDiffPosvsPDirZ->Fill(cosL-p.pz/lengthP);
                        
                        }
                        
                        // TO-DO: Look at contained vs. escaping pions? 
                        //if ((abs(end_pos.z)>60 || abs(end_pos.x)>60) && abs(pdg)==211) escapingPions++;
                        //else if (abs(pdg)==211)  containPiLen->Fill(length);
                        //if (abs(end_pos.z)<60 && abs(end_pos.x)<60 && abs(pdg)==2212) containPLen->Fill(length);
                        //if (cosL<0.9) continue;
                        if (length>minTrkLength){ truePartTrkOnly++;
                        }
                    }
                    truePart++;
                
                } //END OF PRIMARIES LOOP
                //std::cout << "DEBUG: End of primaries loop; start of cuts" << std::endl;
                //if (nKaon>0) continue; // keep if there are no kaons in the event
                //if (nPion>0) continue; // keep if there are no pions in the event
                if (cosL<muon_cos_angle_cut && Elep<muon_energy_cut) continue;
                //escapePi->Fill(escapingPions); 
                //if (truePartTrkOnly>0) true_multTrkOnly->Fill(truePartTrkOnly);    true_multGENIE->Fill(truePartNoG4); responseGenieToG4->Fill(truePart,truePartTrkOnly);
                hasANeutrino=true;  
                trueSignal++;

                //std::cout << "DEBUG: Fill branches in truth tree" << std::endl;
                truth_vtx_x.push_back(trueVtxX);
                truth_vtx_y.push_back(trueVtxY);
                truth_vtx_z.push_back(trueVtxZ);
                truth_muon_energy.push_back(Elep);
                truth_muon_angle.push_back(muon_angle);
                truth_nPionsPrimary.push_back(nPion);
                truth_nProtonsPrimary.push_back(nProton);
                truth_nKaonsPrimary.push_back(nKaon);
                truth_nMuonsPrimary.push_back(nMuon);
                truth_nPi0sPrimary.push_back(nPi0);
                truth_nElectronsPrimary.push_back(nElectron);
                truth_nPhotonsPrimary.push_back(nPhoton);
                truth_found_file.push_back(found_file);
                truth_file_index.push_back(file_num);
                std::cout << Form("Index in cafTree (cafTree->GetEntries()) %ld of %ld", n, nspills) << std::endl;
                std::cout << "DEBUG: Index in sr->mc.nu: " << ntrue << std::endl;
                std::cout << "DEBUG: Spill ID: " << sr->meta.lar2x2.event << std::endl;
                std::cout << "DEBUG: Vertex X: " << trueVtxX << std::endl;
                std::cout << "DEBUG: Vertex Y: " << trueVtxY << std::endl;
                std::cout << "DEBUG: Vertex Z: " << trueVtxZ << std::endl;
                std::cout << "DEBUG: True muon energy: " << Elep << std::endl;
                std::cout << "DEBUG: True muon angle: " << muon_angle << std::endl;
                //std::cout << "DEBUG: nPion: " << nPion << std::endl;
                //std::cout << "DEBUG: nProton: " << nProton << std::endl;
                //std::cout << "DEBUG: nKaon: " << nKaon << std::endl;
                //std::cout << "DEBUG: nMuon: " << nMuon << std::endl;
                //std::cout << "DEBUG: nPi0: " << nPi0 << std::endl;
                //std::cout << "DEBUG: nElectron: " << nElectron << std::endl;
                //std::cout << "DEBUG: nPhoton: " << nPhoton << std::endl;
                //std::cout << "DEBUG: Finished filling branches in truth tree" << std::endl;
            }
        } 
    }
    // RECO first loop over each spill, then each reco interaction, then each reco particle
    std::cout << "Tree name " << tree_name << std::endl;
    std::cout << "Looping over " << nspills << " entries/spills..." << std::endl;
    for(unsigned long i = 0; i < nspills; ++i)
    {
        caf_tree->GetEntry(i);

        auto spill_num = i;

        if(i % incr == 0) {
            std::cout << "Spill #: " << i << std::endl;
        }


        bool save_spill = false;
        std::vector<size_t> event_indices;
        if(is_mc) {// Check if spill_num is in event_ids (for check with JSON dictionary)
            event_indices.clear();
            if (found_file==true and std::find(event_ids.begin(), event_ids.end(), sr->meta.lar2x2.event) != event_ids.end()) {
                std::cout << "Spill is in event " << sr->meta.lar2x2.event << " in event_ids." << std::endl;
                save_spill = true;
                auto it = event_ids.begin();
                while ((it = std::find(it, event_ids.end(), sr->meta.lar2x2.event)) != event_ids.end()) {
                    event_indices.push_back(std::distance(event_ids.begin(), it));
                    ++it; // Move iterator to the next element
                }
            }
        }
        // For truth sample, only save if spill_num is in event_ids
        //if (!save_spill) {
        //    continue;
        //}

        // Loop over each reco interaction
        const auto num_ixn = sr->common.ixn.ndlp;
        bool is_signal = false;
        for(unsigned long ixn = 0; ixn < num_ixn; ++ixn)
        {
            ixnsNoCuts++;
            const auto& vtx = sr->common.ixn.dlp[ixn].vtx;


            // Require RECO vertex to be within the TPCs
            if (!(std::isfinite(vtx.x) && std::isfinite(vtx.y) && std::isfinite(vtx.z)))
                continue;
            if (abs(vtx.x) > fv_abs_x_max)
                continue;
            if (abs(vtx.x) < fv_abs_x_min)
                continue;
            if (abs(vtx.y) > fv_abs_y_max)
                continue;
            if (abs(vtx.z) > fv_abs_z_max)
                continue;
            if (abs(vtx.z) < fv_abs_z_min)
                continue;
            if (abs(vtx.x) > fv_abs_cathode_x_low && abs(vtx.x) < fv_abs_cathode_x_high)
                continue;
            //std::cout << "DEBUG: Passed vertex cuts" << std::endl;
            ixnsVtxCut++;

            // Initialize truth values
            auto truth_idx = -1;
            auto true_ixn_pi0s = -1;
            auto true_ixn_pi0s_contained = -1;
            auto true_ixn_electrons = -1;
            auto true_ixn_gammas = -1;
            auto true_ixn_muons = -1;
            auto true_ixn_chpi = -1;
            auto true_ixn_proton = -1;
            auto true_ixn_chkaon = -1;
            auto true_ixn_vtx_x = -1;
            auto true_ixn_vtx_y = -1;
            auto true_ixn_vtx_z = -1;
            auto true_ixn_isCC = false;
            auto true_ixn_isInFV = false;
            auto true_ixn_targetPDGValue = -1;
            std::string true_ixn_ixnMode = "Other";
            auto true_ixn_nuPDG = -1;
            double current_max = -1; // no overlap bc no truth match

            std::string truth_ixn_dict_key_str = "None";

            // IF MC Get the truth interaction(s) corresponding to this reco interaction
            bool save_ixn = false;
            if(is_mc){
                const auto& vec_truth_ixn = sr->common.ixn.dlp[ixn].truth;
                const auto& vec_overlap_ixn = sr->common.ixn.dlp[ixn].truthOverlap;

                // For truth study, only want to keep if there is a truth match
                // For reco study, I want to know if there is no truth match
                //if(vec_overlap_ixn.empty())
                //    continue;

                //Find the truth interaction with the largest overlap
                //auto result = std::max_element(vec_overlap_ixn.begin(), vec_overlap_ixn.end());
                //auto max_overlap = std::distance(vec_overlap_ixn.begin(), result);

                //Do this manually since SRProxy wraps all types and loses some functionality in the process
                //Maybe there is a way to get the object inside the Proxy, but this works for now
                // Get truth interaction information
                if (!(vec_overlap_ixn.empty()))
                {

                    current_max = 0;
                    unsigned int max_overlap = 0;
                    for(unsigned int i = 0; i < vec_overlap_ixn.size(); ++i)
                    {
                        auto val = vec_overlap_ixn.at(i);
                        if(val > current_max)
                        {
                            current_max = val;
                            max_overlap = i;
                        }
                    }
                    truth_idx = vec_truth_ixn.at(max_overlap);
                    const auto& truth_ixn = sr->mc.nu[truth_idx];

                    // Get truth vertex
                    true_ixn_vtx_x = truth_ixn.vtx.x;
                    true_ixn_vtx_y = truth_ixn.vtx.y;
                    true_ixn_vtx_z = truth_ixn.vtx.z;

                    // Check if vtx is associated with vtx from dictionary truth event
                    if (save_spill == true) {
                        for (size_t event_index : event_indices) {
                            if (std::abs(true_ixn_vtx_x - true_dict_vtx_x[event_index]) < 1.0 && std::abs(true_ixn_vtx_y - true_dict_vtx_y[event_index]) < 1.0 && std::abs(true_ixn_vtx_z - true_dict_vtx_z[event_index]) < 1.0) {
                                save_ixn = true;
                                truth_ixn_dict_key_str = truth_dict_key_strings[event_index];
                            }
                        }
                    }
                
                    //if (!save_ixn) {
                    //    continue;
                    //}
                    //std::cout << "Saving interaction with truth match from dictionary key: " << truth_ixn_dict_key_str << std::endl;'''

                    // Get number of pi0s in truth interaction
                    true_ixn_pi0s = truth_ixn.npi0;
                    true_ixn_pi0s_contained = 0;

                    // Get number of muons, electrons, gammas in truth interaction
                    true_ixn_electrons = 0;
                    true_ixn_gammas = 0;
                    true_ixn_muons = 0;
                    true_ixn_chpi = 0;
                    true_ixn_proton = 0;
                    true_ixn_chkaon = 0;
                    for(unsigned long ipart = 0; ipart < truth_ixn.prim.size(); ++ipart)
                    {
                        const auto& true_part = truth_ixn.prim[ipart];
                        //Put true particle counters
                        if((abs(true_part.pdg) == 11)) {
                            true_ixn_electrons++;
                        }

                        if((abs(true_part.pdg) == 22)) {
                            true_ixn_gammas++;
                        }

                        if((abs(true_part.pdg) == 13))
                            true_ixn_muons++;

                        if((abs(true_part.pdg) == 211))
                            true_ixn_chpi++;

                        if((abs(true_part.pdg) == 2212))
                            true_ixn_proton++;

                        if((abs(true_part.pdg) == 321))
                            true_ixn_chkaon++;

                        // check pi0 containment (?)
                        const auto& true_part_end = true_part.end_pos;
                        if((abs(true_part.pdg) == 111))
                            if (abs(true_part_end.x) > fv_abs_x_min && abs(true_part_end.x) < fv_abs_x_max && abs(true_part_end.y) < fv_abs_y_max && abs(true_part_end.z) > fv_abs_z_min && abs(true_part_end.z) < fv_abs_z_max) {
                                true_ixn_pi0s_contained++;
                            }
                    }
                    if (abs(truth_ixn.pdg)==14 && truth_ixn.iscc==true && truth_ixn.targetPDG==1000180400 && true_ixn_pi0s==1 && abs(vtx.x-true_ixn_vtx_x)<vtx_allowance && abs(vtx.y-true_ixn_vtx_y)<vtx_allowance && abs(vtx.z-true_ixn_vtx_z)<vtx_allowance &&
                        (abs(true_ixn_vtx_x)>abs(mod12_x_min) && abs(true_ixn_vtx_x)<det_x_max && abs(true_ixn_vtx_z)>abs(downs_z_min) && abs(true_ixn_vtx_z)<det_z_max && abs(true_ixn_vtx_y)<det_y_max)) { // && true_ixn_chpi==0 && true_ixn_chkaon==0
                        if (abs(true_ixn_vtx_x)>fv_abs_x_min && abs(true_ixn_vtx_x)<fv_abs_x_max && abs(true_ixn_vtx_z)>fv_abs_z_min && abs(true_ixn_vtx_z)<fv_abs_z_max && abs(true_ixn_vtx_y)<fv_abs_y_max && (!(abs(true_ixn_vtx_x) > fv_abs_cathode_x_low && abs(true_ixn_vtx_x) < fv_abs_cathode_x_high))) {
                                signalNoCuts++;
                                signalVtxCut++;
                                is_signal = true;
                            }
                    }
                    if (abs(true_ixn_vtx_x)>fv_abs_x_min && abs(true_ixn_vtx_x)<fv_abs_x_max && abs(true_ixn_vtx_z)>fv_abs_z_min && abs(true_ixn_vtx_z)<fv_abs_z_max && abs(true_ixn_vtx_y)<fv_abs_y_max &&(!(abs(true_ixn_vtx_x) > fv_abs_cathode_x_low && abs(true_ixn_vtx_x) < fv_abs_cathode_x_high))) {
                        true_ixn_isInFV = true;
                    }
                    true_ixn_isCC = truth_ixn.iscc;
                    true_ixn_targetPDGValue = truth_ixn.targetPDG;
                    true_ixn_nuPDG = truth_ixn.pdg;
                    if (true_ixn_isCC) {
                        if (truth_ixn.mode == 1) true_ixn_ixnMode = "QE";
                        else if (truth_ixn.mode == 10) true_ixn_ixnMode = "MEC";
                        else if (truth_ixn.mode == 4) true_ixn_ixnMode = "RES";
                        else if (truth_ixn.mode == 3) true_ixn_ixnMode = "DIS";
                        else if (truth_ixn.mode == 5) true_ixn_ixnMode = "COH";
                    }
                }
            } // end of MC only section getting truth info

            // Get number of pi0s in reco interaction
            auto reco_ixn_gammas = 0;
            auto reco_ixn_electrons = 0;
            auto reco_ixn_gammas_contained = 0;
            auto reco_ixn_electrons_contained = 0;
            auto reco_ixn_muons = 0;
            auto reco_ixn_chpi = 0;
            auto reco_ixn_proton = 0;
            auto reco_ixn_chkaon = 0;

            int minervaTracks=0; int minervaThrough=0;
            double dirZExiting=-999;   double startZMuonCand=-999;
            double maxDotProductDS=-999; double maxDotProductUS=-999;
            int maxEventPar=-999; int maxEventTyp=-9999; int maxEventIxn=-999;
            double maxShowerEnergy=-999; 
            // Loop over reco particles in the interaction **"The big for loop"**
            for(unsigned long ipart = 0; ipart < sr->common.ixn.dlp[ixn].part.dlp.size(); ++ipart)
            {
                //Store current reco particle for easier access
                const auto& part = sr->common.ixn.dlp[ixn].part.dlp[ipart];

                //Put reco particle cuts here
                if((abs(part.pdg) == 11) and part.primary) {
                    reco_ixn_electrons++;
                    if (part.contained) reco_ixn_electrons_contained++;
                    if (part.E > maxShowerEnergy) maxShowerEnergy=part.E;
                }

                if((abs(part.pdg) == 22) and part.primary) {
                    reco_ixn_gammas++;
                    if (part.contained) reco_ixn_gammas_contained++;
                    if (part.E > maxShowerEnergy) maxShowerEnergy=part.E;
                }

                
                if((abs(part.pdg) == 13) and part.primary){
                    reco_ixn_muons++;
                } 
                if(((abs(part.pdg) == 13) || (abs(part.pdg) == 2212) || (abs(part.pdg) == 211) || (abs(part.pdg) == 321)) and part.primary){
                //if(((abs(part.pdg) == 13)) and part.primary){
                    auto* muon_start = &part.start;
                    auto* muon_end = &part.end;
                    double diffVertexdZ=abs(muon_start->z-sr->common.ixn.dlp[ixn].vtx.z);
	                double diffVertexdX=abs(muon_start->x-sr->common.ixn.dlp[ixn].vtx.x);
	                double diffVertexdY=abs(muon_start->y-sr->common.ixn.dlp[ixn].vtx.y);
	                double diffVertex=TMath::Sqrt(diffVertexdZ*diffVertexdZ+diffVertexdY*diffVertexdY+diffVertexdX*diffVertexdX);
	                //if (diffVertex>5) continue;
                    double dX=(muon_end->x-muon_start->x);
                    double dY=(muon_end->y-muon_start->y);
                    double dZ=(muon_end->z-muon_start->z);
                    double length=TMath::Sqrt(dX*dX+dY*dY+dZ*dZ);
                    double dirX=dX/length; double dirY=dY/length; double dirZ=dZ/length;

                    if (dirZ<0){ dirZ=-dirZ; dirX=-dirX; dirY=-dirY;     
                        muon_start = &part.end;
                        muon_end = &part.start;
                    }
                    if (std::isnan(muon_start->z)) length=-999;
                    int maxPartMinerva=-999; int maxTypeMinerva=-999; int maxIxnMinerva=-999; int maxPartMinervaUS=-999; int maxTypeMinervaUS=-999;
                    // Check if muon goes through downstream Mx2
                    int minervaPass=0;
                    double dotProductDS=-999; double deltaExtrapYUS=-999; double deltaExtrapY=-999; 
                    double dotProductUS=-999; double deltaExtrapX=-999; double deltaExtrapXUS=-999;	
                    if(abs(muon_start->z)>fv_abs_z_max || abs(muon_end->z)>fv_abs_z_max){
                        //std::cout << "DEBUG: Muon start Z: " << muon_start->z << std::endl;
                        //std::cout << "DEBUG: Muon end Z: " << muon_end->z << std::endl;
                        //std::cout << "DEBUG: Mx2 ixn size: " << sr->nd.minerva.ixn.size() << std::endl;
                	 	for(int i=0; i<sr->nd.minerva.ixn.size(); i++){
                        
                		    for (int j=0; j<sr->nd.minerva.ixn[i].ntracks; j++){
                		        double dir_z=sr->nd.minerva.ixn[i].tracks[j].dir.z;
                		        double end_z=sr->nd.minerva.ixn[i].tracks[j].end.z;
                		        double start_z=sr->nd.minerva.ixn[i].tracks[j].start.z;
                		        double end_x=sr->nd.minerva.ixn[i].tracks[j].end.x;
                		        double start_x=sr->nd.minerva.ixn[i].tracks[j].start.x;
                		        double end_y=sr->nd.minerva.ixn[i].tracks[j].end.y;
                		        double start_y=sr->nd.minerva.ixn[i].tracks[j].start.y;
                                //std::cout << "DEBUG: Track start Z: " << start_z << std::endl;
                        
                		        if (start_z>0 && ((muon_start->z)>fv_abs_z_max || (muon_end->z)>fv_abs_z_max) ){
                		            int truthPart=sr->nd.minerva.ixn[i].tracks[j].truth[0].part;
                		            double dXMnv=(sr->nd.minerva.ixn[i].tracks[j].end.x-sr->nd.minerva.ixn[i].tracks[j].start.x);
                		            double dYMnv=(sr->nd.minerva.ixn[i].tracks[j].end.y-sr->nd.minerva.ixn[i].tracks[j].start.y);
                		            double dZMnv=(sr->nd.minerva.ixn[i].tracks[j].end.z-sr->nd.minerva.ixn[i].tracks[j].start.z);
                		            double lengthMinerva=TMath::Sqrt(dXMnv*dXMnv+dYMnv*dYMnv+dZMnv*dZMnv);
                                    //std::cout << "DEBUG: Length Minerva: " << lengthMinerva << std::endl;
                	                if (lengthMinerva<10) continue;
                          	        double dirXMinerva=dXMnv/lengthMinerva;
                		            double dirYMinerva=dYMnv/lengthMinerva;
                		            double dirZMinerva=dZMnv/lengthMinerva;
                		            double dotProduct=dirXMinerva*dirX+dirYMinerva*dirY+dirZ*dirZMinerva;
                		            double extrapdZ=start_z-muon_end->z+mnvOffsetZDS;
                		            double extrapY=dirY/dirZ*(extrapdZ)+muon_end->y-start_y;
                		            double extrapX=dirX/dirZ*(extrapdZ)+muon_end->x-start_x;
                		            double diffExtrap=TMath::Sqrt(TMath::Power(extrapY-start_y,2));
                        
                        
                		            if (dotProductDS<dotProduct && abs(extrapY-mnvOffsetY)<mnvDistanceAllowance  && abs(TMath::ATan(dirXMinerva/dirZMinerva)-TMath::ATan(dirX/dirZ))<mnvAngleAllowance && abs(TMath::ATan(dirYMinerva/dirZMinerva)-TMath::ATan(dirY/dirZ))<mnvAngleAllowance && abs(extrapX-mnvOffsetX)<mnvDistanceAllowance){ 
                                        dotProductDS=dotProduct;
                		                deltaExtrapY=extrapY;
                		                deltaExtrapX=extrapX;
                		                dirZExiting=dirZ;
                		                if (is_mc){
                                            maxPartMinerva=sr->nd.minerva.ixn[i].tracks[j].truth[0].part;
                		                    maxTypeMinerva=sr->nd.minerva.ixn[i].tracks[j].truth[0].type;
                                            maxIxnMinerva=sr->nd.minerva.ixn[i].tracks[j].truth[0].ixn;
                                        }	
                		                //if (end_z>300){ minervaPass=1;} if(dirZExiting<dirZ){ dirZExiting=dirZ;}
                                    }
                                }
                            
                		        if (start_z<0 && end_z>0 && ( (muon_start->z<det_z_min && muon_end->z>det_z_max) || (muon_start->z>det_z_max && muon_end->z<det_z_min))){
                		            int truthPart=sr->nd.minerva.ixn[i].tracks[j].truth[0].part;
                		            double dXMnv=(sr->nd.minerva.ixn[i].tracks[j].end.x-sr->nd.minerva.ixn[i].tracks[j].start.x);
                		            double dYMnv=(sr->nd.minerva.ixn[i].tracks[j].end.y-sr->nd.minerva.ixn[i].tracks[j].start.y);
                		            double dZMnv=(sr->nd.minerva.ixn[i].tracks[j].end.z-sr->nd.minerva.ixn[i].tracks[j].start.z);
                		            double lengthMinerva=TMath::Sqrt(dXMnv*dXMnv+dYMnv*dYMnv+dZMnv*dZMnv);
                		            double dirXMinerva=dXMnv/lengthMinerva;
                		            double dirYMinerva=dYMnv/lengthMinerva;
                		            double dirZMinerva=dZMnv/lengthMinerva;
                		            double dotProduct=dirXMinerva*dirX+dirYMinerva*dirY+dirZ*dirZMinerva;
                                            
                		            double extrapdZUS=end_z-muon_end->z;
                		            double extrapYUS=dirY/dirZ*(extrapdZUS)+muon_end->y-end_y;
                		            double extrapXUS=dirX/dirZ*(extrapdZUS)+muon_end->x-end_x;
                        
                                    //double diffExtrap=TMath::Sqrt(TMath::Power(extrapY-end_y,2));
                		            if (dotProductUS<dotProduct  && abs(extrapYUS-mnvOffsetY)<mnvDistanceAllowance && abs(extrapXUS-mnvOffsetX)<mnvDistanceAllowance && abs(TMath::ATan(dirXMinerva/dirZMinerva)-TMath::ATan(dirX/dirZ))<mnvAngleAllowance && abs(TMath::ATan(dirYMinerva/dirZMinerva)-TMath::ATan(dirY/dirZ))<mnvAngleAllowance) dotProductUS=dotProduct;
                		            deltaExtrapYUS=extrapYUS;
                                    deltaExtrapXUS=extrapXUS;
                		            if (is_mc){
                                        maxPartMinervaUS=sr->nd.minerva.ixn[i].tracks[j].truth[0].part;
                		                maxTypeMinervaUS=sr->nd.minerva.ixn[i].tracks[j].truth[0].type;	
                                    }
                		        }
                            } // Loop over minerva tracks
                        } // Minerva event loop
                	    if (dotProductDS>maxDotProductDS){ maxDotProductDS=dotProductDS;
                            maxEventPar=maxPartMinerva;
                            maxEventTyp=maxTypeMinerva;
                            maxEventIxn=maxIxnMinerva;
                	    }
                	    if (dotProductDS>mnvMatchDotProdCut){  minervaTracks++; 
                            if (minervaPass==1){ minervaThrough++;
                                startZMuonCand=muon_start->z; if (muon_start->z>muon_end->z) startZMuonCand=muon_end->z;
                            }
                	    } 
                    } // end of if muon goes through minerva
                
                    if (dotProductUS>maxDotProductUS) maxDotProductUS=dotProductUS;
                //if (dotProductUS>0.99){   
           /*		if (maxPartNumber==maxPartMinervaUS && maxTypeMinervaUS==maxTypeNumber){ deltaYGoodUS->Fill(deltaExtrapYUS); deltaXGoodUS->Fill(deltaExtrapXUS); goodMINERvAMatchUS++;} 
                         else{ deltaYBadUS->Fill(deltaExtrapYUS); deltaXBadUS->Fill(deltaExtrapXUS); }
                  totalMINERvAMatchUS++;     
                  */	//}

                   // Get True Direction Muon Angle

                }
                if((abs(part.pdg) == 211) and part.primary)
                    reco_ixn_chpi++;

                if((abs(part.pdg) == 2212) and part.primary)
                    reco_ixn_proton++;

                if((abs(part.pdg) == 321) and part.primary)
                    reco_ixn_chkaon++;

            } // END OF BIG FOR LOOP FOR NOW (loop over reco particles)

            // Requirements on reco particles
            // Only save interactions with exactly 1 primary muon

            //if (reco_ixn_muons != 1) continue;
            //ixnsMuonCut++;
            //std::cout << "DEBUG: Passed muon cut" << std::endl;
            //if (is_mc && is_signal) signalMuonCut++;
            // Only save if muon match to Mx2 throughgoing track
            if ( /*minervaThrough<1  ||*/  maxDotProductDS<mnvMatchDotProdCut) continue;
            ixnsMx2Cut++;
            if (is_mc && is_signal) signalMx2Cut++;
            //std::cout << "DEBUG: Passed minerva cut" << std::endl;
            
            // Require two showers (photon or electron bc bad PID)
            //std::cout << "DEBUG: Reco ixn gammas + electrons: " << reco_ixn_gammas << " // " << reco_ixn_electrons << std::endl;
            if((reco_ixn_gammas + reco_ixn_electrons) != 2) continue;

            ixnsShowerCut++;
            if (is_mc && is_signal) signalShowerCut++;
            //std::cout << "DEBUG: Passed shower cut" << std::endl;

            // Require no charged pions or kaons
            /*if(reco_ixn_chpi != 0)
                continue;

            ixnsPionCut++;
            if (is_mc && is_signal) signalPionCut++;


            if(reco_ixn_chkaon != 0)
                continue;

            ixnsKaonCut++;
            if (is_mc && is_signal) signalKaonCut++;*/

            auto muonAngle = 0.; auto trueMuonAngle = 0.;
            auto recoLeadShowerEnergy = 0.; auto recoLeadShowerPX = 0.; auto recoLeadShowerPY = 0.; auto recoLeadShowerPZ = 0.; auto recoLeadShowerPDG = 0.; 
            auto recoSubleadShowerEnergy = 0.; auto recoSubleadShowerPX = 0.; auto recoSubleadShowerPY = 0.; auto recoSubleadShowerPZ = 0.; auto recoSubleadShowerPDG = 0.; 
            auto recoLeadShowerStartX = 0.; auto recoLeadShowerStartY = 0.; auto recoLeadShowerStartZ = 0.;
            auto recoSubleadShowerStartX = 0.; auto recoSubleadShowerStartY = 0.; auto recoSubleadShowerStartZ = 0.;
            auto trueLeadShowerEnergy = 0.; auto trueLeadShowerPX = 0.; auto trueLeadShowerPY = 0.; auto trueLeadShowerPZ = 0.; auto trueLeadShowerPDG = 0.; auto trueLeadShowerOvlp = 0;
            auto trueLeadShowerStartX = 0.; auto trueLeadShowerStartY = 0.; auto trueLeadShowerStartZ = 0.;
            auto trueSubleadShowerStartX = 0.; auto trueSubleadShowerStartY = 0.; auto trueSubleadShowerStartZ = 0.;
            auto trueSubleadShowerEnergy = 0.; auto trueSubleadShowerPX = 0.; auto trueSubleadShowerPY = 0.; auto trueSubleadShowerPZ = 0.; auto trueSubleadShowerPDG = 0.; auto trueSubleadShowerOvlp = 0;
            bool isLeadShower = false; 
            auto true_pvec = TVector3(-99999, -99999, -99999);
            auto true_dir = TVector3(-99999, -99999, -99999);
            //auto true_shower_to_true_vtx = -99999;
            // Loop over particles again to get particle kinematics info
            for(unsigned long ipart = 0; ipart < sr->common.ixn.dlp[ixn].part.dlp.size(); ++ipart)
            {
                //Store current reco particle for easier access
                const auto& part = sr->common.ixn.dlp[ixn].part.dlp[ipart];

                // Reset isLeadShower
                isLeadShower = false;
                if (maxShowerEnergy == part.E) {
                    isLeadShower = true;
                }

                // Get Reco Direction Muon Angle
                if (abs(part.pdg) == 13 && part.primary){

                    auto pvec = TVector3(part.p.x, part.p.y, part.p.z);
                    auto pvec_unit = pvec.Unit();
                    muonAngle = pvec_unit.Angle(beam_dir);

                }

                if ((abs(part.pdg) == 22 || abs(part.pdg) == 11) && part.primary) {
                    
                    if (isLeadShower) {
                        recoLeadShowerEnergy = part.E;
                        recoLeadShowerPX = part.p.x;
                        recoLeadShowerPY = part.p.y;
                        recoLeadShowerPZ = part.p.z;
                        recoLeadShowerPDG = part.pdg;
                        recoLeadShowerStartX = part.start.x;
                        recoLeadShowerStartY = part.start.y;
                        recoLeadShowerStartZ = part.start.z;
                    }
                    else {
                        recoSubleadShowerEnergy = part.E;
                        recoSubleadShowerPX = part.p.x;
                        recoSubleadShowerPY = part.p.y;
                        recoSubleadShowerPZ = part.p.z;
                        recoSubleadShowerPDG = part.pdg;
                        recoSubleadShowerStartX = part.start.x;
                        recoSubleadShowerStartY = part.start.y;
                        recoSubleadShowerStartZ = part.start.z;
                    }
                }

                // Get Shower information 
                // Get truth match(es) for this reco particle
                // Variables need to be wrapped in the Proxy object (sometimes)
                if (is_mc){
                    caf::Proxy<caf::SRTrueParticle>* truth_match = nullptr;
                    const auto& vec_truth_id = part.truth;
                    const auto& vec_overlap  = part.truthOverlap;
/*
                    ////If the truth overlap vector is empty, then assume no truth match and skip
                    //if(vec_overlap.empty())
                    //{
                    //    //std::cout << "No truth match... skipping reco particle..." << std::endl;
                    //    continue;
                    //}
*/  
                    //auto true_length_val = -99999;
                    double part_current_max = -1;
                    if (!(vec_overlap.empty()))
                    {//Find the truth particle with the largest overlap
                        //auto result = std::max_element(vec_overlap.begin(), vec_overlap.end());
                        //auto max_overlap = std::distance(vec_overlap.begin(), result);
                        //auto truth_id = vec_truth_id.at(max_overlap);
                        part_current_max = 0;
                        unsigned int max_overlap = 0;
                        for(unsigned int i = 0; i < vec_overlap.size(); ++i)
                        {
                            auto val = vec_overlap.at(i);
                            if(val > part_current_max)
                            {
                                part_current_max = val;
                                max_overlap = i;
                            }
                        }
//  
                    ////if(current_max < 0.25)
                    ////    continue;
                        const auto& truth_id = vec_truth_id.at(max_overlap);
//  
                        //Get pointer to the corresponding truth particle
                        if(truth_id.type == 1)
                            truth_match = &(sr->mc.nu[truth_id.ixn].prim[truth_id.part]);
                        else if(truth_id.type == 3)
                            truth_match = &(sr->mc.nu[truth_id.ixn].sec[truth_id.part]);
                        else
                        {
                            std::cout << "Invalid truth id type!" << std::endl;
                            continue;
                        }
                    
                    
                        true_pvec = TVector3(truth_match->p.px, truth_match->p.py, truth_match->p.pz);
                        true_dir = true_pvec.Unit();

                        if (abs(part.pdg) == 13 && part.primary){
                            trueMuonAngle = true_dir.Angle(beam_dir);
                        }

                        if ((abs(part.pdg) == 22 || abs(part.pdg) == 11) && part.primary) {
                            if (isLeadShower) {
                                trueLeadShowerEnergy = truth_match->p.E;
                                trueLeadShowerPX = truth_match->p.px;
                                trueLeadShowerPY = truth_match->p.py;
                                trueLeadShowerPZ = truth_match->p.pz;
                                trueLeadShowerPDG = truth_match->pdg;
                                trueLeadShowerOvlp = part_current_max;
                                trueLeadShowerStartX = truth_match->start_pos.x;
                                trueLeadShowerStartY = truth_match->start_pos.y;
                                trueLeadShowerStartZ = truth_match->start_pos.z;
                            }
                            else {
                                trueSubleadShowerEnergy = truth_match->p.E;
                                trueSubleadShowerPX = truth_match->p.px;
                                trueSubleadShowerPY = truth_match->p.py;
                                trueSubleadShowerPZ = truth_match->p.pz;
                                trueSubleadShowerPDG = truth_match->pdg;
                                trueSubleadShowerOvlp = part_current_max;
                                trueSubleadShowerStartX = truth_match->start_pos.x;
                                trueSubleadShowerStartY = truth_match->start_pos.y;
                                trueSubleadShowerStartZ = truth_match->start_pos.z;
                            }
                        }
                    }
                        //true_shower_to_true_vtx = (TVector3(truth_match->start_pos.x, truth_match->start_pos.y, truth_match->start_pos.z) - TVector3(true_ixn_vtx_x, true_ixn_vtx_y, true_ixn_vtx_z)).Mag();
                }
            } // End of second loop over reco particles  

            
            true_match_lead_shower_energy.push_back(trueLeadShowerEnergy);
            true_match_lead_shower_p_x.push_back(trueLeadShowerPX); 
            true_match_lead_shower_p_y.push_back(trueLeadShowerPY); 
            true_match_lead_shower_p_z.push_back(trueLeadShowerPZ);
            true_match_sublead_shower_energy.push_back(trueSubleadShowerEnergy);
            true_match_sublead_shower_p_x.push_back(trueSubleadShowerPX);  
            true_match_sublead_shower_p_y.push_back(trueSubleadShowerPY); 
            true_match_sublead_shower_p_z.push_back(trueSubleadShowerPZ);
            true_match_lead_shower_start_x.push_back(trueLeadShowerStartX);
            true_match_lead_shower_start_y.push_back(trueLeadShowerStartY);
            true_match_lead_shower_start_z.push_back(trueLeadShowerStartZ);
            true_match_sublead_shower_start_x.push_back(trueSubleadShowerStartX);
            true_match_sublead_shower_start_y.push_back(trueSubleadShowerStartY);
            true_match_sublead_shower_start_z.push_back(trueSubleadShowerStartZ);
            true_match_muon_angle.push_back(trueMuonAngle);
            true_match_lead_shower_pdg.push_back(trueLeadShowerPDG);
            true_match_sublead_shower_pdg.push_back(trueSubleadShowerPDG);
            true_match_lead_ovlp.push_back(trueLeadShowerOvlp);
            true_match_sublead_ovlp.push_back(trueSubleadShowerOvlp);
            true_ixn_pi0_mult.push_back(true_ixn_pi0s);
            true_ixn_e_mult.push_back(true_ixn_electrons);
            true_ixn_gamma_mult.push_back(true_ixn_gammas);
            true_ixn_cont_pi0_mult.push_back(true_ixn_pi0s_contained);
            true_ixn_muon_mult.push_back(true_ixn_muons);
            true_ixn_chpi_mult.push_back(true_ixn_chpi);
            true_ixn_proton_mult.push_back(true_ixn_proton);
            true_ixn_chkaon_mult.push_back(true_ixn_chkaon);
            true_ixn_vtx_x_pos.push_back(true_ixn_vtx_x);
            true_ixn_vtx_y_pos.push_back(true_ixn_vtx_y);
            true_ixn_vtx_z_pos.push_back(true_ixn_vtx_z);
            true_ixn_is_cc.push_back(true_ixn_isCC);
            true_ixn_is_in_fv.push_back(true_ixn_isInFV);
            true_ixn_target.push_back(true_ixn_targetPDGValue);
            true_ixn_mode.push_back(true_ixn_ixnMode);
            true_ixn_nu_pdg.push_back(true_ixn_nuPDG);
            true_ixn_is_signal.push_back(is_signal);
           
            //std::cout << "DEBUG: Filling Reco TTree values" << std::endl;
            reco_lead_shower_energy.push_back(recoLeadShowerEnergy);
            reco_lead_shower_p_x.push_back(recoLeadShowerPX); 
            reco_lead_shower_p_y.push_back(recoLeadShowerPY); 
            reco_lead_shower_p_z.push_back(recoLeadShowerPZ);
            reco_sublead_shower_energy.push_back(recoSubleadShowerEnergy);
            reco_sublead_shower_p_x.push_back(recoSubleadShowerPX); 
            reco_sublead_shower_p_y.push_back(recoSubleadShowerPY); 
            reco_sublead_shower_p_z.push_back(recoSubleadShowerPZ);
            reco_lead_shower_start_x.push_back(recoLeadShowerStartX);
            reco_lead_shower_start_y.push_back(recoLeadShowerStartY);
            reco_lead_shower_start_z.push_back(recoLeadShowerStartZ);
            reco_sublead_shower_start_x.push_back(recoSubleadShowerStartX);
            reco_sublead_shower_start_y.push_back(recoSubleadShowerStartY);
            reco_sublead_shower_start_z.push_back(recoSubleadShowerStartZ);
            reco_muon_angle.push_back(muonAngle);
            reco_muon_angle_Mx2_match.push_back(dirZExiting);
            reco_lead_shower_pdg.push_back(recoLeadShowerPDG);
            reco_sublead_shower_pdg.push_back(recoSubleadShowerPDG);
            reco_ixn_e_mult.push_back(reco_ixn_electrons);
            reco_ixn_gamma_mult.push_back(reco_ixn_gammas);
            reco_ixn_e_cont_mult.push_back(reco_ixn_electrons_contained);
            reco_ixn_cont_gamma_mult.push_back(reco_ixn_gammas_contained);
            reco_ixn_muon_mult.push_back(reco_ixn_muons);
            reco_ixn_chpi_mult.push_back(reco_ixn_chpi);
            reco_ixn_proton_mult.push_back(reco_ixn_proton);
            reco_ixn_chkaon_mult.push_back(reco_ixn_chkaon);
            reco_ixn_vtx_x_pos.push_back(vtx.x);
            reco_ixn_vtx_y_pos.push_back(vtx.y);
            reco_ixn_vtx_z_pos.push_back(vtx.z);
            overlap.push_back(current_max);
            true_ixn_index.push_back(truth_idx);
            reco_ixn_index.push_back(ixn);
            spill_index.push_back(sr->meta.lar2x2.event);
            file_index.push_back(file_num);
            event.push_back(sr->meta.nd_lar.event);
            run.push_back(sr->meta.lar2x2.run);
            subrun.push_back(sr->meta.lar2x2.subrun);
            caf_file_name.push_back(current_file.erase(0, current_file.find_last_of("/")+1).c_str());
            in_truth_dict.push_back(save_ixn);
            truth_dict_key.push_back(truth_ixn_dict_key_str);
        } // end of loop over interactions
    } //end of reco spills

    } //end of file loop

    // Push back PurEff values
    if(is_mc){
        total_ixns.push_back(trueInteractions);
        total_signal.push_back(trueSignal);
        signal_reco_no_cuts.push_back(signalNoCuts);
        all_reco_no_cuts.push_back(ixnsNoCuts);
        signal_post_vtx_cut.push_back(signalVtxCut);
        all_post_vtx_cut.push_back(ixnsVtxCut);
        signal_post_muon_cut.push_back(signalMuonCut);
        all_post_muon_cut.push_back(ixnsMuonCut);
        signal_post_mx2_cut.push_back(signalMx2Cut);
        all_post_mx2_cut.push_back(ixnsMx2Cut);
        signal_post_shower_cut.push_back(signalShowerCut);
        all_post_shower_cut.push_back(ixnsShowerCut);
        signal_post_pion_cut.push_back(signalPionCut);
        all_post_pion_cut.push_back(ixnsPionCut);
        signal_post_kaon_cut.push_back(signalKaonCut);
        all_post_kaon_cut.push_back(ixnsKaonCut);
    }

    const auto t_end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> t_elapsed{t_end - t_start};

        // Output TTree file name
    std::string file_name = "first_pass_general_CC1pi0_selection_SANDBOX_v11beta_with_mesons_fv_cut_xy2cm_z3cm_mx2_any_track_PID_no_single_muon_cut";
    //std::string file_name = "first_pass_general_CC1pi0_selection_SANDBOX_v6_with_mesons_fv_cut_xy2cm_z3cm_cathode2cm_mx2_any_track_PID_no_single_muon_cut";
    //std::string file_name = "first_pass_general_CC1pi0_selection_MR6p4_1000_files_with_mesons_fv_cut_xy2cm_z3cm_cathode2cm_mx2_any_track_PID_no_single_muon_cut";
    //std::string file_name = "MR6p4_debug_multiple_truth_entries";
    // DEFINE: Output TFile
    TFile *f=new TFile(Form("%s.root", file_name.c_str()),"RECREATE");


    // POPULATE: Fill Truth TTree and write to output ROOT file
    if(is_mc){
        fTruthTree->Fill();
        fTruthTree->Write();
        std::cout << "Filled and wrote Truth TTree." << std::endl;

        fPurEffTree->Fill();
        fPurEffTree->Write();
        std::cout << "Filled and wrote PurEff TTree." << std::endl;
    }

    // POPULATE: Fill Reco TTree and write to output ROOT file
    fRecoTree->Fill();
    fRecoTree->Write();
    
    std::cout << "Filled and wrote Reco TTree." << std::endl;
    std::cout << "Number of entries in Reco TTree: " << fRecoTree->GetEntries() << std::endl;

    // CLOSE: Output ROOT file
    f->Close();

    std::cout << "Time elapsed: " << t_elapsed.count() << std::endl;
    std::cout << "Finished." << std::endl;
    return 0;
}
