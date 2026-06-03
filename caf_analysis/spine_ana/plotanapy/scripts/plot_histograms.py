
from plotanapy.plotting import HistogramLoader, RecoPlotter
from plotanapy.plotting.styles import apply_theme
import argparse


if __name__ == '__main__':
    # parsers for any cmd line args
    parser = argparse.ArgumentParser()
    parser.add_argument("--root_file", default='./outputs/testCC1pi0SelSPINE.root', type=str)
    parser.add_argument("--theme", default='default', type=str)
    parser.add_argument("--all_plots_dir", default='./outputs', type=str, help='directory to store plots')
    parser.add_argument("--mcOnly", default='1', type=str, help='Analysis MC only? 1=True, 0=False')
    args = parser.parse_args()

    # Apply theme
    apply_theme(args.theme)

    # Load histograms and create plotter
    print(f"Loading histograms from {args.root_file}...")
    loader = HistogramLoader(args.root_file)
    reco_plotter = RecoPlotter(loader.get_reco_histograms())
    # cuflow_plotter = CuflowPlotter(loader.get_cuflow_histograms())
    #if (args.mcOnly == '1'):
    #    truth_plotter = TruthPlotter(loader.get_truth_histograms())
    #    truth_matched_plotter = TruthMatchedPlotter(loader.get_truth_matched_histograms())
    
    # Generate all plots at once
    print("Plotting reco histograms...")
    reco_plotter.plot_all(output_dir=args.all_plots_dir+"reco_plots")
