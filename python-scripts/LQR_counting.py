import argparse
import os
import re
from dataclasses import dataclass
from pathlib import Path


@dataclass
class BedParse:
    contig: str
    start: str
    stop: str
    fwd_cov: str
    rev_cov: str

    def __str__(self) -> str:
        return '\t'.join([self.contig, self.start, self.stop, self.fwd_cov, self.rev_cov])


def parse_bed(qv: int, input_dir: Path, output_dir: Path):
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
        inp = os.path.join(input_dir, input_filename)
        out = os.path.join(output_dir, input_filename).replace("sequtils", "LQR")
        with open(inp, "r") as f, open(out, "w") as lqr_file:
            table_data = [line.split() for line in f]
            for i in table_data:
                if (float(i[4]) < qv) or ((float(i[4]) >= qv) and (float(i[5]) < qv)):
                    bed_string = BedParse(i[0], i[1], str(int(i[2])+1), i[4], i[5])
                    # If the directory for input and output files contains "sequtils" -> error
                    lqr_file.write(str(bed_string) + '\n')


def main(quality_threshold, input_dir, output_dir):
    parse_bed(quality_threshold, input_dir, output_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script for LQR counting")
    parser.add_argument("-q", "--quality_threshold", type=int, help="The quality threshold")
    parser.add_argument("-i", "--input_dir", type=lambda p: Path(p).absolute(),
                        help="The path to input files directory")
    parser.add_argument("-o", "--output_dir", type=lambda p: Path(p).absolute(),
                        help="The path to output files directory")
    args = parser.parse_args()
    main(quality_threshold=args.quality_threshold, input_dir=args.input_dir, output_dir=args.output_dir)
