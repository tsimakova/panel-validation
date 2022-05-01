import argparse
import matplotlib.pyplot as plt
import pandas as pd
import pathlib
import seaborn as sns


def table_for_plot(first_point: int, last_point: int, points: int, input_file: pathlib.PosixPath) -> pd.DataFrame:
    """
    :param first_point: The first point among numbers of reads per amplicon
    :param last_point: The last point among numbers of reads per amplicon
    :param points: The number of points
    :param input_file: The path to an input file
    :return pd.DataFrame for subsequent plotting
    """
    df = pd.DataFrame()
    lqr_prop = pd.read_csv(input_file, sep="\t", header=None)
    step = (last_point - first_point) // (points - 1)
    read_per_amp = []
    for i in range(first_point, last_point + step, step):
        read_per_amp.append(i)
    df["Reads per amplicon"] = read_per_amp
    df["Proportion of LQRs"] = lqr_prop[0]
    return df


def lqr_plot(df: pd.DataFrame, output_dir: pathlib.PosixPath):
    """
    :param df: pd.DataFrame containing the proportion of LQRs for each point
    :param output_dir: The path to an output file directory
    """
    sns.set(rc={'figure.figsize': (15, 8)})
    lqr = sns.lineplot(data=df, x="Reads per amplicon", y="Proportion of LQRs", color="purple")
    lqr.set(title="Percentage of LQR in a panel target regions")
    lqr.set(ylabel="Percentage", xlabel="Reads per amplicon")
    plt.savefig(f"{output_dir}/LQR_proportion_plot.png")


def parse_args():
    parser = argparse.ArgumentParser(description="Script for plotting the LQR proportion for each point")
    parser.add_argument("-f", "--first_point", type=int, help="The first point among numbers of reads per amplicon")
    parser.add_argument("-l", "--last_point", type=int, help="The last point among numbers of reads per amplicon")
    parser.add_argument("-p", "--points", type=int, help="The number of points")
    parser.add_argument("-i", "--input_file", help="The path to input file")
    parser.add_argument("-o", "--output_dir", type=lambda p: pathlib.Path(p).absolute(),
                        help="The path to output file directory")
    return parser.parse_args()


def main():
    args = parse_args()
    df_for_plot = table_for_plot(args.first_point, args.last_point, args.points, args.input_file)
    lqr_plot(df_for_plot, args.output_dir)


if __name__ == "__main__":
    main()
