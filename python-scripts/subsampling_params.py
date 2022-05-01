import argparse
import json
import pathlib


def get_params(first_point: int, last_point: int, points: int, amp_number: int, mapped_reads: int,
               output_dir: pathlib.PosixPath):
    """
    :param first_point: The first point among numbers of reads per amplicon
    :param last_point: The last point among numbers of reads per amplicon
    :param points: The number of points
    :param amp_number: The number of amplicons in a panel
    :param mapped_reads: The number of mapped reads in a BAM file
    :param output_dir: The path to output files
    """
    # Count the number of reads per amplicon for each point and create a dictionary:
    step = (last_point - first_point) // (points - 1)
    read_per_amp = []
    for i in range(first_point, last_point + step, step):
        read_per_amp.append(i)
    amplicon_dict = dict.fromkeys(read_per_amp, None)
    # Count the number of reads per sample for each point:
    for el in amplicon_dict.keys():
        amplicon_dict[el] = el * amp_number
    # Count the percentage of reads per sample for each point:
    for el in amplicon_dict.keys():
        amplicon_dict[el] = amplicon_dict[el] / mapped_reads
    # Make a .json file with parameters for a subsequent subsampling
    with open(f"{output_dir}/subsampling_params.json", 'w') as f:
        json.dump(amplicon_dict, f)


def parse_args():
    parser = argparse.ArgumentParser(description="Script for subsampling parameters")
    parser.add_argument("-f", "--first_point", type=int, help="The first point among numbers of reads per amplicon")
    parser.add_argument("-l", "--last_point", type=int, help="The last point among numbers of reads per amplicon")
    parser.add_argument("-p", "--points", type=int, help="The number of points")
    parser.add_argument("-a", "--amp_number", type=int, help="The number of amplicons in a panel")
    parser.add_argument("-m", "--mapped_reads", type=int, help="The number of mapped reads in a .bam file")
    parser.add_argument("-o", "--output_dir", type=lambda p: pathlib.Path(p).absolute(),
                        help="The path to output files")
    return parser.parse_args()


def main():
    args = parse_args()
    get_params(args.first_point, args.last_point, args.points, args.amp_number, args.mapped_reads, args.output_dir)


if __name__ == "__main__":
    main()
