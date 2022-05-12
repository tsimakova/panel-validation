import argparse
import os
import pathlib
import re


def parse_bed(qv: int, input_dir: pathlib.PosixPath, output_dir: pathlib.PosixPath):
    """
    :param qv: The quality threshold
    :param input_dir: The path to input files directory
    :param output_dir: The path to output files directory
    """
    files = os.listdir(input_dir)
    filename = []
    for file in files:
        if file.endswith('sequtils.bed'):
            filename.append(file)
    # Sort the "filename" list according to the number of reads per amplicon:
    filename.sort(key=lambda test_string: list(map(int, re.findall(r'\d+', test_string)))[-1])
    for input_filename in filename:
        dict_from_bam = {"contig": [], "start": [], "stop": [], "fwd_cov": [], "rev_cov": []}
        with open(f"{input_dir}/{input_filename}") as f:
            table_data = [line.split() for line in f]
            for i in table_data:
                if (float(i[4]) < qv) or ((float(i[4]) >= qv) and (float(i[5]) < qv)):
                    dict_from_bam["contig"].append(i[0])
                    dict_from_bam["start"].append(i[1])
                    dict_from_bam["stop"].append(str(int(i[2])+1))
                    dict_from_bam["fwd_cov"].append(i[4])
                    dict_from_bam["rev_cov"].append(i[5] + "\n")
        # If the directory for input and output files contains "sequtils" -> error
        with open(f"{output_dir}/{input_filename}".replace("sequtils", "LQR"), "w") as lqr_file:
            for value in zip(dict_from_bam["contig"], dict_from_bam["start"], dict_from_bam["stop"],
                             dict_from_bam["fwd_cov"], dict_from_bam["rev_cov"]):
                lqr_file.write('\t'.join(value))


def main(quality_threshold, input_dir, output_dir):
    parse_bed(quality_threshold, input_dir, output_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script for LQR counting")
    parser.add_argument("-q", "--quality_threshold", type=int, help="The quality threshold")
    parser.add_argument("-i", "--input_dir", type=lambda p: pathlib.Path(p).absolute(),
                        help="The path to input files directory")
    parser.add_argument("-o", "--output_dir", type=lambda p: pathlib.Path(p).absolute(),
                        help="The path to output files directory")
    args = parser.parse_args()
    main(quality_threshold=args.quality_threshold, input_dir=args.input_dir, output_dir=args.output_dir)
