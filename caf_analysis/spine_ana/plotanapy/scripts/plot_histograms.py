
from plotanapy.plotting import HistogramLoader, RecoPlotter, TruthPlotter
from plotanapy.plotting.styles import apply_theme
from plotanapy.plotting.utils import ParseYAML
import argparse


if __name__ == '__main__':
    # parsers for any cmd line args
    parser = argparse.ArgumentParser()
    parser.add_argument("--root_file", default='./outputs/testCC1pi0SelSPINE.root', type=str)
    parser.add_argument("--config_file", default='./include/config_data/CC1pi0.yaml', type=str)
    parser.add_argument("--theme", default='default', type=str)
    parser.add_argument("--all_plots_dir", default='./outputs', type=str, help='directory to store plots')
    parser.add_argument("--mcOnly", default='1', type=str, help='Analysis MC only? 1=True, 0=False')
    args = parser.parse_args()

    # Get config parameters from yaml file
    sel_config = ParseYAML(args.config_file, config_name="selection")
    beam_config = ParseYAML(args.config_file, config_name="beam")
    det_config = ParseYAML(args.config_file, config_name="detector")

    # Apply theme
    apply_theme(args.theme)

    # Load histograms and create plotters
    print(f"Loading histograms from {args.root_file}...")
    loader = HistogramLoader(args.root_file)
    reco_plotter = RecoPlotter(loader.get_reco_histograms(), \
                               sel_config=sel_config, \
                               beam_config=beam_config, \
                               det_config=det_config)
    # cuflow_plotter = CuflowPlotter(loader.get_cuflow_histograms())
    if (args.mcOnly == '1'):
        truth_plotter = TruthPlotter(loader.get_truth_histograms(), \
                                     sel_config=sel_config, \
                                     beam_config=beam_config, \
                                     det_config=det_config)
    #    truth_matched_plotter = TruthMatchedPlotter(loader.get_truth_matched_histograms())
    
    # Generate all plots at once
    print("Plotting reco histograms...")
    reco_plotter.plot_all(output_dir=args.all_plots_dir+"reco_plots")
    #print("Plotting cuflow histograms...")
    #cuflow_plotter.plot_all(output_dir=args.all_plots_dir+"cuflow_plots")
    if (args.mcOnly == '1'):
        print("Plotting truth histograms...")
        truth_plotter.plot_all(output_dir=args.all_plots_dir+"truth_plots")
    #    print("Plotting truth matched histograms...")
    #    truth_matched_plotter.plot_all(output_dir=args.all_plots_dir+"truth_matched_plots")
    
