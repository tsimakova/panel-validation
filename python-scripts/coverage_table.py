import argparse
import numpy as np
import pandas as pd
import pathlib


def cov_table(first_point: int, last_point: int, points: int, amp_number: int, output_dir: pathlib.PosixPath):
    """
    :param first_point: The first point among numbers of reads per amplicon
    :param last_point: The last point among numbers of reads per amplicon
    :param points: The number of points
    :param amp_number: The number of amplicons in a panel
    :param output_dir: The path to output files
    """
    # Create an empty pd.Dataframe:
    df = pd.DataFrame(np.zeros((points, points)))
    # Define the row names (number of reads per amplicon):
    step = (last_point - first_point) // (points - 1)
    read_per_amp = []
    for i in range(first_point, last_point + step, step):
        read_per_amp.append(i)
    df.index = read_per_amp
    # Define the column names (number of reads per sample):
    first_col = df.index[0] * amp_number
    last_col = df.index[-1] * amp_number
    step = (last_col - first_col) // (points - 1)
    df.columns = range(first_col, last_col + step, step)
    # Fill the pd.Dataframe:
    for i in range(len(df.index)):
        for j in range(len(df.columns)):
            if df.columns[j] / (amp_number * df.index[i]) * 100 >= 100:
                df.iat[i, j] = round(100, 2)
            else:
                df.iat[i, j] = round(df.columns[j] / (amp_number * df.index[i]) * 100, 2)
    # Save coverage table:
    df.to_csv(f"{output_dir}/coverage_table.txt", sep="\t", index=True, header=True)


def parse_args():
    parser = argparse.ArgumentParser(description="Script for coverage table creation")
    parser.add_argument("-f", "--first_point", type=int, help="The first point among numbers of reads per amplicon")
    parser.add_argument("-l", "--last_point", type=int, help="The last point among numbers of reads per amplicon")
    parser.add_argument("-p", "--points", type=int, help="The number of points")
    parser.add_argument("-a", "--amp_number", type=int, help="The number of amplicons in a panel")
    parser.add_argument("-o", "--output_dir", type=lambda p: pathlib.Path(p).absolute(),
                        help="The path to output files")
    return parser.parse_args()


def main():
    args = parse_args()
    cov_table(args.first_point, args.last_point, args.points, args.amp_number, args.output_dir)


if __name__ == "__main__":
    main()
