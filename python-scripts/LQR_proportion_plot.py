import argparse
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os.path
from pathlib import Path


def table_for_plot(first_point: int, last_point: int, points: int, amp_number: int, mapped_reads: int, input_file:
                   Path) -> pd.DataFrame:
    """
    :param first_point: The first point among numbers of reads per amplicon
    :param last_point: The last point among numbers of reads per amplicon
    :param points: The number of points
    :param amp_number: The number of amplicons in a panel
    :param mapped_reads: The number of mapped reads in a BAM file
    :param input_file: The path to an input file
    :return pd.DataFrame for subsequent plotting
    """
    # Check the first and last points
    if last_point - first_point <= 0:
        raise ValueError("The last point is less than or equal to the first point")
    # Check that the multiplication of the maximum number of reads per amplicon and the number of amplicons in a panel
    # does not exceed the number of reads in the bam file
    if mapped_reads // amp_number < last_point:
        last_point = mapped_reads // amp_number
        print(f"The last point (maximum number of reads per amplicon) was changed to {last_point}")
    df = pd.DataFrame()
    lqr_prop = pd.read_csv(input_file, sep="\t", header=None)
    step = (last_point - first_point) // (points - 1)
    read_per_amp = []
    for i in range(first_point, last_point + 1, step):
        read_per_amp.append(i)
    df["Reads_per_amplicon"] = read_per_amp
    df["Proportion_of_LQRs"] = lqr_prop[0]
    return df


def lqr_plot(df: pd.DataFrame, output_dir: Path, figure_width: float, figure_height: float):
    """
    :param figure_width: The width of the LQR plot
    :param figure_height: The height of the LQR plot
    :param df: pd.DataFrame containing the proportion of LQRs for each point
    :param output_dir: The path to an output file directory
    """
    sns.set(rc={'figure.figsize': (figure_width, figure_height)})
    lqr = sns.lineplot(data=df, x="Reads_per_amplicon", y="Proportion_of_LQRs", color="purple")
    lqr.set_title("Percentage of LQR in panel target regions", fontsize=20)
    lqr.set_ylabel("Percentage", fontsize=15)
    lqr.set_xlabel("Number of reads per amplicon", fontsize=15)
    plt.savefig(os.path.join(output_dir, "LQR_proportion_plot.png"))


def main(first_point, last_point, points, amp_number, mapped_reads, input_file, output_dir, figure_width, figure_height):
    df_for_plot = table_for_plot(first_point, last_point, points, amp_number, mapped_reads, input_file)
    lqr_plot(df_for_plot, output_dir, figure_width, figure_height)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script for plotting the LQR proportion for each point")
    parser.add_argument("-f", "--first_point", type=int, help="The first point among numbers of reads per amplicon")
    parser.add_argument("-l", "--last_point", type=int, help="The last point among numbers of reads per amplicon")
    parser.add_argument("-p", "--points", type=int, help="The number of points")
    parser.add_argument("-a", "--amp_number", type=int, help="The number of amplicons in a panel")
    parser.add_argument("-m", "--mapped_reads", type=int, help="The number of mapped reads in a .bam file")
    parser.add_argument("-i", "--input_file", help="The path to input file")
    parser.add_argument("-o", "--output_dir", type=lambda p: Path(p).absolute(),
                        help="The path to output file directory")
    parser.add_argument("-w", "--figure_width", type=float, default=10, help="The width of the LQR plot")
    parser.add_argument("-e", "--figure_height", type=float, default=6, help="The height of the LQR plot")
    args = parser.parse_args()
    main(first_point=args.first_point, last_point=args.last_point, points=args.points, amp_number=args.amp_number,
         mapped_reads=args.mapped_reads, input_file=args.input_file, output_dir=args.output_dir,
         figure_width=args.figure_width, figure_height=args.figure_height)
