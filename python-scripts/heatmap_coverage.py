import argparse
import os.path
import matplotlib.pyplot as plt
import pandas as pd
import pathlib
import seaborn as sns


def heatmap(input_file: pathlib.PosixPath, output_dir: pathlib.PosixPath, figure_width: float, figure_height: float):
    """
    :param input_file: The path to input file
    :param output_dir: The path to output file directory
    :param figure_width: The width of the heatmap
    :param figure_height: The height of the heatmap
    """
    df = pd.read_csv(input_file, sep="\t", header=0, index_col=0)
    plt.figure(figsize=(figure_width, figure_height))
    plot = sns.heatmap(df, annot=True, fmt='.3g', cmap="RdYlGn")
    plot.set_xlabel('Reads per sample', fontsize=15)
    plot.set_ylabel('Reads per amplicon', fontsize=15)
    plot.set_title('Proportion of amplicons (%) with the target coverage', fontsize=15)
    plt.savefig(os.path.join(output_dir, "heatmap_coverage.png"))


def main(input_file, output_dir, figure_width, figure_height):
    heatmap(input_file, output_dir, figure_width, figure_height)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script for coverage heatmap")
    parser.add_argument("-i", "--input_file", type=lambda p: pathlib.Path(p).absolute(),
                        help="The path to input file")
    parser.add_argument("-o", "--output_dir", type=lambda p: pathlib.Path(p).absolute(),
                        help="The path to output files")
    parser.add_argument("-w", "--figure_width", type=float, default=10,
                        help="The width of the heatmap")
    parser.add_argument("-e", "--figure_height", type=float, default=6,
                        help="The height of the heatmap")
    args = parser.parse_args()
    main(input_file=args.input_file, output_dir=args.output_dir, figure_width=args.figure_width,
         figure_height=args.figure_height)
