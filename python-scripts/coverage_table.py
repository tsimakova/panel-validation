import argparse
import numpy as np
import pandas as pd
import pathlib
import os.path


def cov_table(first_point: int, last_point: int, points: int, amp_number: int, correction: int,
              output_dir: pathlib.PosixPath):
    """
    :param first_point: The first point among numbers of reads per amplicon
    :param last_point: The last point among numbers of reads per amplicon
    :param points: The number of points
    :param amp_number: The number of amplicons in a panel
    :param correction: Correction coefficient
    :param output_dir: The path to output files
    """
    # Check the first and last points
    if last_point - first_point <= 0:
        raise ValueError("The last point is less than or equal to the first point")
    # Create an empty pd.Dataframe:
    df = pd.DataFrame(np.zeros((points, points)))
    # Define the row names (number of reads per amplicon):
    step = (last_point - first_point) // (points - 1)
    read_per_amp = []
    for i in range(first_point, last_point + 1, step):
        read_per_amp.append(i)
    df.index = read_per_amp
    # Define the column names (number of reads per sample):
    first_col = int(df.index[0] * amp_number * (1 + correction/100))
    last_col = int(df.index[-1] * amp_number * (1 + correction/100))
    step = (last_col - first_col) // (points - 1)
    df.columns = range(first_col, last_col + 1, step)
    # Fill the pd.Dataframe:
    for i in range(len(df.index)):
        for j in range(len(df.columns)):
            if (df.columns[j] / (1 + correction/100)) / (amp_number * df.index[i]) * 100 >= 100:
                df.iat[i, j] = round(100, 2)
            else:
                df.iat[i, j] = round((df.columns[j] / (1 + correction/100)) / (amp_number * df.index[i]) * 100, 2)
    # Save coverage table:
    df.to_csv(os.path.join(output_dir, "coverage_table.txt"), sep="\t", index=True, header=True)


def main(first_point, last_point, points, amp_number, correction, output_dir):
    cov_table(first_point, last_point, points, amp_number, correction, output_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script for coverage table creation")
    parser.add_argument("-f", "--first_point", type=int, help="The first point among numbers of reads per amplicon")
    parser.add_argument("-l", "--last_point", type=int, help="The last point among numbers of reads per amplicon")
    parser.add_argument("-p", "--points", type=int, help="The number of points")
    parser.add_argument("-a", "--amp_number", type=int, help="The number of amplicons in a panel")
    parser.add_argument("-c", "--correction", type=int, help="Correction coefficient")
    parser.add_argument("-o", "--output_dir", type=lambda p: pathlib.Path(p).absolute(),
                        help="The path to output files")
    args = parser.parse_args()
    main(first_point=args.first_point, last_point=args.last_point, points=args.points, amp_number=args.amp_number,
         correction=args.correction, output_dir=args.output_dir)
