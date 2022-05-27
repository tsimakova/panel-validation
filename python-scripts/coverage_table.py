import argparse
import numpy as np
import pandas as pd
import pathlib
import os.path


def cov_table(tsv_files: list, first_point: int, last_point: int, points: int, amp_number: int, mapped_reads: int,
              correction: int, output_dir: pathlib.PosixPath):
    """
    :param tsv_files: The list of tsv files, obtained after coverage analysis
    :param first_point: The first point among numbers of reads per amplicon
    :param last_point: The last point among numbers of reads per amplicon
    :param points: The number of points
    :param amp_number: The number of amplicons in a panel
    :param mapped_reads: The number of mapped reads in a BAM file
    :param correction: Correction coefficient
    :param output_dir: The path to output files
    """
    # Check the first and last points
    if last_point - first_point <= 0:
        raise ValueError("The last point is less than or equal to the first point")
    # Check that the multiplication of the maximum number of reads per amplicon and the number of amplicons in a panel
    # does not exceed the number of reads in the bam file
    if mapped_reads // amp_number < last_point:
        last_point = mapped_reads // amp_number
        print(f"The last point (maximum number of reads per amplicon) was changed to {last_point}")
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

    for i in range(len(tsv_files)):
        temp_data = pd.read_csv(tsv_files[i], sep="\t", header=0)
        new_column = []
        for elem in df.index:
            sum = 0
            for el in temp_data[["total_reads"]].to_numpy():
                if int(el) >= elem:
                    sum += 1
            new_column.append(round(sum / amp_number * 100))
        se = pd.Series(new_column)
        df[df.columns[i]] = se.values

    df.to_csv(os.path.join(output_dir, "coverage_table.txt"), sep="\t", index=True, header=True)


def main(tsv_files, first_point, last_point, points, amp_number, mapped_reads, correction, output_dir):
    cov_table(tsv_files, first_point, last_point, points, amp_number, mapped_reads, correction, output_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script for coverage table creation")
    parser.add_argument("-t", "--tsv_files", nargs="+", type=lambda p: pathlib.Path(p).absolute(),
                        help="The list of tsv files, obtained after coverage analysis")
    parser.add_argument("-f", "--first_point", type=int, help="The first point among numbers of reads per amplicon")
    parser.add_argument("-l", "--last_point", type=int, help="The last point among numbers of reads per amplicon")
    parser.add_argument("-p", "--points", type=int, help="The number of points")
    parser.add_argument("-a", "--amp_number", type=int, help="The number of amplicons in a panel")
    parser.add_argument("-m", "--mapped_reads", type=int, help="The number of mapped reads in a .bam file")
    parser.add_argument("-c", "--correction", type=int, help="Correction coefficient")
    parser.add_argument("-o", "--output_dir", type=lambda p: pathlib.Path(p).absolute(),
                        help="The path to output files")
    args = parser.parse_args()
    main(tsv_files=args.tsv_files, first_point=args.first_point, last_point=args.last_point, points=args.points,
         amp_number=args.amp_number, mapped_reads=args.mapped_reads, correction=args.correction,
         output_dir=args.output_dir)
