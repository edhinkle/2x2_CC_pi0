
from plotdqmpy.plotting import HistogramLoader, BeamPlotter
from plotdqmpy.plotting.styles import apply_theme
from plotdqmpy.plotting.utils import ParseYAML
import argparse


if __name__ == '__main__':
    # parsers for any cmd line args
    parser = argparse.ArgumentParser()
    parser.add_argument("--root_file", default='./outputs/Physics_DQM_Run1_Sandboxv11.root', type=str)
    parser.add_argument("--sample_name", default='Sandboxv11', type=str, help='Name of the sample being analyzed')
    parser.add_argument("--theme", default='default', type=str)
    parser.add_argument("--all_plots_dir", default='./outputs', type=str, help='directory to store plots')
    parser.add_argument("--mcOnly", default='0', type=str, help='Analysis MC only? 1=True, 0=False')
    args = parser.parse_args()

    # Apply theme
    apply_theme(args.theme)

    # Load histograms and create plotters
    print(f"Loading histograms from {args.root_file}...")
    loader = HistogramLoader(args.root_file)
    beam_plotter = BeamPlotter(loader.get_beam_histograms(), sample_name=args.sample_name)

    # Generate all plots at once
    print("Plotting beam histograms...")
    beam_plotter.plot_all(output_dir=args.all_plots_dir+"beam_plots")

