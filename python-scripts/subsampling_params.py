import argparse
import json
import os.path
from pathlib import Path


def get_params(first_point: int, last_point: int, points: int, amp_number: int, mapped_reads: int,
               output_dir: Path):
    """
    :param first_point: The first point among numbers of reads per amplicon
    :param last_point: The last point among numbers of reads per amplicon
    :param points: The number of points
    :param amp_number: The number of amplicons in a panel
    :param mapped_reads: The number of mapped reads in a BAM file
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
    # Count the number of reads per amplicon for each point and create a dictionary:
    step = (last_point - first_point) // (points - 1)
    read_per_amp = []
    for i in range(first_point, last_point + 1, step):
        read_per_amp.append(i)
    amplicon_dict = dict.fromkeys(read_per_amp, None)
    # Count the number of reads per sample for each point:
    for el in amplicon_dict.keys():
        amplicon_dict[el] = el * amp_number
    # Count the percentage of reads per sample for each point:
    for el in amplicon_dict.keys():
        amplicon_dict[el] = amplicon_dict[el] / mapped_reads
    # Make a .json file with parameters for a subsequent subsampling
    with open(os.path.join(output_dir, "subsampling_params.json"), 'w') as f:
        json.dump(amplicon_dict, f)


def main(first_point, last_point, points, amp_number, mapped_reads, output_dir):
    get_params(first_point, last_point, points, amp_number, mapped_reads, output_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script for subsampling parameters")
    parser.add_argument("-f", "--first_point", type=int, help="The first point among numbers of reads per amplicon")
    parser.add_argument("-l", "--last_point", type=int, help="The last point among numbers of reads per amplicon")
    parser.add_argument("-p", "--points", type=int, help="The number of points")
    parser.add_argument("-a", "--amp_number", type=int, help="The number of amplicons in a panel")
    parser.add_argument("-m", "--mapped_reads", type=int, help="The number of mapped reads in a .bam file")
    parser.add_argument("-o", "--output_dir", type=lambda p: Path(p).absolute(),
                        help="The path to output files")
    args = parser.parse_args()
    main(first_point=args.first_point, last_point=args.last_point, points=args.points, amp_number=args.amp_number,
         mapped_reads=args.mapped_reads, output_dir=args.output_dir)
