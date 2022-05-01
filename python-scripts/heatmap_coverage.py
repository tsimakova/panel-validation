import argparse
import matplotlib.pyplot as plt
import pandas as pd
import pathlib
import seaborn as sns


def heatmap(input_file: pathlib.PosixPath, output_dir: pathlib.PosixPath):
    """
    :param input_file: The path to input file
    :param output_dir: The path to output file directory
    """
    df = pd.read_csv(input_file, sep="\t", header=0, index_col=0)
    plt.figure(figsize=(16, 5))
    plot = sns.heatmap(df, annot=True, fmt='.3g', cmap="RdYlGn")
    plot.set_xlabel('Reads per sample', fontsize=15)
    plot.set_ylabel('Reads per amplicon', fontsize=15)
    plot.set_title('Proportion of amplicons (%) with the target coverage', fontsize=15)
    plt.savefig(f"{output_dir}/heatmap_coverage.png")


def parse_args():
    parser = argparse.ArgumentParser(description="Script for coverage heatmap")
    parser.add_argument("-i", "--input_file", type=lambda p: pathlib.Path(p).absolute(),
                        help="The path to input file")
    parser.add_argument("-o", "--output_dir", type=lambda p: pathlib.Path(p).absolute(),
                        help="The path to output files")
    return parser.parse_args()


def main():
    args = parse_args()
    heatmap(args.input_file, args.output_dir)


if __name__ == "__main__":
    main()
