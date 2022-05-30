import subprocess
import os
import time


# 1. Enter the range and points number
def reads_per_amplicon(first, last, points):
    """
    :param first: the first point
    :param last: the last point
    :param points: number of points
    :return: an array containing the number of reads per amplicon
    """
    step = (last - first) // (points - 1)
    read_per_amp = []
    for i in range(first, last + step, step):
        read_per_amp.append(i)
    return read_per_amp


# 2. Count the number of reads per sample for each point
def reads_per_sample(read_per_amp, amp_num):
    """
    :param read_per_amp: an array containing the number of reads per amplicon == range_and_point_number() output
    :param amp_num: the number of amplicon in a panel
    :return: the list containing the number of reads per sample
    """
    read_per_sample = []
    for el in read_per_amp:
        read_per_sample.append(amp_num * el)
    return read_per_sample


# 3. Count the percentage of reads in a .BAM file
def percentage(read_per_sample):
    """
    :param read_per_sample: an array containing the number of reads per sample == reads_per_sample() output
    :return: an array containing the percentage of reads (for subsampling)
    """
    files = os.listdir(".")
    bam_files = []
    for file in files:
        if file.endswith('genome_alignment.bam'):
            bam_files.append(file)
    my_stdin = open(bam_files[0])
    command = "samtools view -c -F 4"
    cmd = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stdin=my_stdin)
    mapped_reads = int(cmd.stdout.decode())
    perc = []
    for el in read_per_sample:
        perc.append(el / mapped_reads)
    return perc


# 4. BAM subsampling
def bam_subsampling(bam, read_per_amp, perc):
    """
    :param bam: .bam file for subsampling
    :param read_per_amp: an array containing the number of reads per amplicon
    :param perc: an array containing the percentage for subsampling
    """
    split_bam_file_name = str(bam).split(".bam")[0]
    # for el in range(len(perc)):
    #     perc[el] = str(perc[el])
    for el in range(len(perc)):
        subprocess.Popen(["samtools", "view", "-s", str(perc[el]), "-b", bam, "-o",
                          f"{split_bam_file_name}_sub_{read_per_amp[el]}.bam"])


# 5. Run Sequtils.jar
def sequtils(target_region):
    """
    :param target_region: .bed file containing the target panel regions
    """
    files = os.listdir(".")
    bam_files = []
    split_bam_file_name = []
    for file in files:
        if "genome_alignment_sub_" in file:
            bam_files.append(file)
            split_bam_file_name.append(str(file).split(".bam")[0])
    for el in range(len(bam_files)):
        subprocess.Popen(["java", "-jar", "sequtils.jar", "regions", "-t", target_region, "-b",
                          f"{bam_files[el]}", "-o", f"{split_bam_file_name[el]}_sequtils.bed"])


# 6. Count LQR
def count_lqr(qv):
    """
    :param qv: quality threshold
    """
    files = os.listdir(".")
    filename = []
    for file in files:
        if file.endswith('sequtils.bed'):
            filename.append(file)
    for input_filename in filename:
        chrom = []
        start = []
        end = []
        fwd_cov = []
        rev_cov = []
        with open(input_filename) as f:
            table_data = [line.split() for line in f]
            for i in table_data:
                if float(i[4]) < qv:
                    chrom.append(i[0])
                    start.append(i[1])
                    end.append(str(int(i[2]) + 1))
                    fwd_cov.append(i[4])
                    rev_cov.append(i[5] + "\n")
                if (float(i[4]) >= qv) and (float(i[5]) < qv):
                    chrom.append(i[0])
                    start.append(i[1])
                    end.append(str(int(i[2]) + 1))
                    fwd_cov.append(i[4])
                    rev_cov.append(i[5] + "\n")
        lqr_file = open(input_filename.replace("sequtils", "LQR"), "w")
        for row in zip(chrom, start, end, fwd_cov, rev_cov):
            lqr_file.write('\t'.join(row))
        lqr_file.close()


if __name__ == "__main__":
    reads_amp = reads_per_amplicon(100, 1000, 6)
    reads_sample = reads_per_sample(reads_amp, 181)
    per = percentage(reads_sample)
    bam_subsampling("sub_5000_NSCLC-HR-TruQ4_S23.genome_alignment.bam", reads_amp, per)
    # без задержки у меня не успевали появиться файлы после subsampling и sequtils
    time.sleep(1)
    sequtils("nsclc__target_regions.bed")
    time.sleep(35)
    count_lqr(64)
