import json
import os
import subprocess
import time


rule all:
    input:
        expand(os.path.join("{run_dir}", "OUTPUT_FILES", "{sample}_{type}_amplicons.txt"), sample = config["cov_analysis_result"],
                            run_dir=config["run_dir"], type=config["amplicon_type"]) ,
        os.path.join(config["run_dir"], "temporal_files", f'{config["bam_sample"]}_number_of_mapped_reads.txt') ,
        os.path.join(config["run_dir"], "temporal_files", "subsampling_params.json") , #
        expand(os.path.join("{run_dir}","temporal_files", "subsampling_results", "{bam_sample}_sub{index}.bam"),
                            bam_sample=config["bam_sample"], run_dir=config["run_dir"],index=range(int(config["points"]))) ,
        expand(os.path.join("{run_dir}","temporal_files", "sequtils_results", "{bam_sample}_sub{index}_sequtils.bed"),
                            bam_sample=config["bam_sample"], run_dir=config["run_dir"],
                            index=range(int(config["points"]))) ,
        expand(os.path.join("{run_dir}","temporal_files", "LQR_count_results", "{bam_sample}_sub{index}_LQR.bed"),
                            bam_sample=config["bam_sample"], run_dir=config["run_dir"],
                            index=range(int(config["points"]))) ,
        os.path.join(config["run_dir"], "temporal_files", "total_number_of_lqr_positions.txt") ,
        os.path.join(config["run_dir"], "temporal_files", "lqr_proportions.txt") ,
        os.path.join(config["run_dir"], "OUTPUT_FILES", "LQR_proportion_plot.png") ,
        expand(os.path.join("{run_dir}", "temporal_files", "subsampling_results", "{bam_sample}_sub{index}.bam.bai"),
                            run_dir=config["run_dir"], bam_sample=config["bam_sample"], index=range(int(config["points"]))) ,
        expand(os.path.join("{run_dir}", "temporal_files", "coverage_analysis_results",
                            "{bam_sample}_sub{index}.amplicon.cov.tsv"), run_dir=config["run_dir"],
                            bam_sample=config["bam_sample"], index=range(int(config["points"]))) ,
        os.path.join(config["run_dir"], "OUTPUT_FILES", "coverage_table.txt") ,
        os.path.join(config["run_dir"], "OUTPUT_FILES", "heatmap_coverage.png")


rule amplicon_coverage:
    input:
        expand(os.path.join("{run_dir}", "{sample}.tsv"), sample = config["cov_analysis_result"],
                            run_dir=config["run_dir"])
    output:
        expand(os.path.join("{run_dir}", "OUTPUT_FILES", "{sample}_{type}_amplicons.txt"),
                            sample = config["cov_analysis_result"], run_dir=config["run_dir"], type=config["amplicon_type"])
    message:
        "Look for under- and overcovered amplicons and create a linear regression plot"
    params:
        script_path=os.path.join(config["scripts_dir"], "amplicon_coverage.py"),
        threshold=config["threshold"],
        under_ratio=config["under_ratio"],
        over_ratio=config["over_ratio"],
        output_dir=os.path.join(config["run_dir"], "OUTPUT_FILES"),
        run_dir=config["run_dir"],
        width=config["lin_reg_width"],
        height=config["lin_reg_height"]
    shell:
        """
        python3 {params.script_path} -t {params.threshold} -u {params.under_ratio} -o {params.over_ratio} -i {input} \
        -d {params.output_dir} -w {params.width} -e {params.height}
        """


rule count_mapped_reads:
    input:
        os.path.join(config["run_dir"], f'{config["bam_sample"]}.bam')
    output:
        os.path.join(config["run_dir"], "temporal_files", f'{config["bam_sample"]}_number_of_mapped_reads.txt')
    message:
        "Run samtools to count the number of mapped reads in a BAM file"
    shell:
        """
        samtools view -c -F 4 {input} > {output}
        """


def return_number_of_mapped_reads(wildcards):
    """
    :return: the number of mapped reads
    """
    with open(os.path.join(config["run_dir"], "temporal_files", f'{config["bam_sample"]}_number_of_mapped_reads.txt'),
              "r") as bam_reads:
        for line in bam_reads:
            mapped_reads = int(line)
    return mapped_reads


rule params_for_subsampling:
    input:
        os.path.join(config["run_dir"], "temporal_files", f'{config["bam_sample"]}_number_of_mapped_reads.txt')
    output:
        os.path.join(config["run_dir"], "temporal_files", "subsampling_params.json")
    message:
        "Count the percentage of reads in files after subsampling"
    params:
        script_path=os.path.join(config["scripts_dir"], "subsampling_params.py"),
        f_point=config["first_point"],
        l_point=config["last_point"],
        points=config["points"],
        amp_number=config["amp_number"],
        m_reads= lambda wildcards: return_number_of_mapped_reads(wildcards),
        output_dir=os.path.join(config["run_dir"], "temporal_files")
    shell:
        """
        python3 {params.script_path} -f {params.f_point} -l {params.l_point} -p {params.points} -a {params.amp_number} \
        -m {params.m_reads} -o {params.output_dir}
        """


def return_json_data(wildcards):
    """
    :return: the list of percentage for subsampling
    """
    with open(os.path.join(config["run_dir"], "temporal_files", "subsampling_params.json"), "r") as js_data:
        a = json.load(js_data)
        json_values = []
        for elem in a:
            json_values.append(a[elem])
        return json_values


rule bam_subsampling:
    input:
        os.path.join(config["run_dir"], f'{config["bam_sample"]}.bam') ,
        rules.params_for_subsampling.output
    output:
        expand(os.path.join("{run_dir}", "temporal_files", "subsampling_results", "{bam_sample}_sub{index}.bam"), run_dir=config["run_dir"],
                            bam_sample=config["bam_sample"], index=range(int(config["points"])))
    message:
        "Use samtools for BAM file subsampling"
    params:
        path = return_json_data
    run:
        for el in range(len(params.path)):
            subprocess.Popen(["samtools", "view", "-s", str(params.path[el]), "-b", input[0], "-o", output[el]])
        time.sleep(5)


rule run_sequtils:
    input:
        expand(os.path.join("{run_dir}", "temporal_files", "subsampling_results", "{bam_sample}_sub{index}.bam"),
                            run_dir=config["run_dir"], bam_sample=config["bam_sample"], index=range(int(config["points"])))
    output:
        expand(os.path.join("{run_dir}", "temporal_files", "sequtils_results", "{bam_sample}_sub{index}_sequtils.bed"),
                            bam_sample=config["bam_sample"], run_dir=config["run_dir"], index=range(int(config["points"])))
    message:
        "Run sequtils"
    params:
        sequtils_path = config["path_to_sequtils"],
        target_path = os.path.join(config["run_dir"], config["tagret_regions"]),
        points = config["points"]
    run:
        for el in range(params.points):
            subprocess.Popen(["java", "-jar", params.sequtils_path, "regions", "-t", params.target_path, "-b",
                              input[el], "-o", output[el]])
            time.sleep(20)
        time.sleep(100)


rule count_lqr:
    input:
        expand(os.path.join("{run_dir}", "temporal_files", "sequtils_results", "{bam_sample}_sub{index}_sequtils.bed"),
                            bam_sample=config["bam_sample"], run_dir=config["run_dir"], index=range(int(config["points"])))
    output:
        expand(os.path.join("{run_dir}", "temporal_files", "LQR_count_results", "{bam_sample}_sub{index}_LQR.bed"),
                             run_dir=config["run_dir"], bam_sample=config["bam_sample"], index=range(int(config["points"])))
    message:
        "Count the positions with low sequence quality"
    params:
        qv = config["qv"],
        script_path = os.path.join(config["scripts_dir"], "LQR_counting.py"),
        run_dir= os.path.join(config["run_dir"], "temporal_files", "sequtils_results"),
        out_dir= os.path.join(config["run_dir"], "temporal_files",  "LQR_count_results")
    shell:
        """
        python3 {params.script_path} -q {params.qv} -i {params.run_dir} -o {params.out_dir}
        """


rule count_total_number_positions:
    input:
        os.path.join(config["run_dir"], config["tagret_regions"])
    output:
        os.path.join(config["run_dir"], "temporal_files", "total_number_of_lqr_positions.txt")
    message:
        "Count the number of positions (nucleotides) in a BED file with target regions"
    shell:
        """
        cat {input} | awk -F"\t" "BEGIN{{SUM=0}}{{ SUM+=\$3-\$2 }}END{{print SUM}}" > {output}
        """


def return_total_number_positions(wildcards):
    """
    :return: the number of analysed positions
    """
    with open(os.path.join(config["run_dir"], "temporal_files", "total_number_of_lqr_positions.txt"), "r") as total_position_number:
        for line in total_position_number:
            num = int(line)
    return num


rule lqr_proportion:
    input:
        expand(os.path.join("{run_dir}", "temporal_files", "LQR_count_results", "{bam_sample}_sub{index}_LQR.bed"),
                             run_dir=config["run_dir"], bam_sample=config["bam_sample"], index=range(int(config["points"])))
    output:
        os.path.join(config["run_dir"], "temporal_files", "lqr_proportions.txt")
    message:
        "Count the proportion of positions with low sequence quality"
    params:
        num = lambda wildcards: return_total_number_positions(wildcards),
        points = config["points"]
    run:
        for el in input:
            c = """awk -F"\t" "BEGIN{{SUM=0}}{{ SUM+=\$3-\$2 }}END{{print SUM}}" {el} | awk "{{print \$1/{params.num}}}" >> {output}"""
            shell(c)


rule create_LQRs_plot:
    input:
        os.path.join(config["run_dir"], "temporal_files", "lqr_proportions.txt")
    output:
        os.path.join(config["run_dir"], "OUTPUT_FILES", "LQR_proportion_plot.png")
    message:
        "Create a LQR proportion lineplot"
    params:
        script_path = os.path.join(config["scripts_dir"], "LQR_proportion_plot.py"),
        f_point=config["first_point"],
        l_point=config["last_point"],
        amp_number=config["amp_number"],
        m_reads= lambda wildcards: return_number_of_mapped_reads(wildcards),
        points=config["points"],
        output_dir=os.path.join(config["run_dir"], "OUTPUT_FILES"),
        width = config["lqr_plot_width"],
        height = config["lqr_plot_height"]
    shell:
        """
        python3 {params.script_path} -f {params.f_point} -l {params.l_point} -p {params.points} -a \
        {params.amp_number} -m {params.m_reads} -i {input} -o {params.output_dir} -w {params.width} -e {params.height}
        """


rule make_bai:
    input:
        expand(os.path.join("{run_dir}", "temporal_files", "subsampling_results", "{bam_sample}_sub{index}.bam"),
                             run_dir=config["run_dir"], bam_sample=config["bam_sample"], index=range(int(config["points"])))
    output:
        expand(os.path.join("{run_dir}", "temporal_files", "subsampling_results", "{bam_sample}_sub{index}.bam.bai"),
                             run_dir=config["run_dir"], bam_sample=config["bam_sample"], index=range(int(config["points"])))
    message:
        "Create a BAI index"
    run:
        for el in input:
            command = """samtools index -b {el}"""
            shell(command)


rule coverage_analysis:
    input:
        bam = expand(os.path.join("{run_dir}", "temporal_files", "subsampling_results", "{bam_sample}_sub{index}.bam"),
                                   run_dir=config["run_dir"], bam_sample=config["bam_sample"], index=range(int(config["points"]))) ,
        bai = expand(os.path.join("{run_dir}", "temporal_files", "subsampling_results", "{bam_sample}_sub{index}.bam.bai"),
                                   run_dir=config["run_dir"], bam_sample=config["bam_sample"], index=range(int(config["points"])))
    output:
        expand(os.path.join("{run_dir}", "temporal_files", "coverage_analysis_results", "{bam_sample}_sub{index}.amplicon.cov.tsv"),
                                        run_dir=config["run_dir"], bam_sample=config["bam_sample"],
                                        index=range(int(config["points"])))
    message:
        "Run coverage analysis"
    params:
        ref = config["ref_genome"] ,
        cov = config["path_to_analyze_cov"] ,
        amp = os.path.join(config["run_dir"], config["designed_amplicons"]) ,
        out_dir = os.path.join(config["run_dir"], "temporal_files", "coverage_analysis_results")
    run:
        for el in input.bam:
            command = """cwltool --outdir {params.out_dir} {params.cov} --aligned_to_genome_indexed_reads {el} \
                         --designed_bed {params.amp} --reference_genome {params.ref}"""
            shell(command)


rule create_coverage_table:
    input:
        map_reads = os.path.join(config["run_dir"], "temporal_files", f'{config["bam_sample"]}_number_of_mapped_reads.txt') ,
        tsv_files = expand(os.path.join("{run_dir}", "temporal_files", "coverage_analysis_results",
                           "{bam_sample}_sub{index}.amplicon.cov.tsv"), run_dir=config["run_dir"],
                           bam_sample=config["bam_sample"], index=range(int(config["points"])))
    output:
        os.path.join(config["run_dir"], "OUTPUT_FILES", "coverage_table.txt")
    message:
        "Create a coverage table"
    params:
        script_path = os.path.join(config["scripts_dir"], "coverage_table.py"),
        f_point=config["first_point"],
        l_point = config["last_point"],
        points = config["points"],
        amp =  config["amp_number"],
        m_reads= lambda wildcards: return_number_of_mapped_reads(wildcards),
        corr=config["correction_coeff"],
        output_dir=os.path.join(config["run_dir"], "OUTPUT_FILES"),
    shell:
        """
        python3 {params.script_path} -t {input.tsv_files} -f {params.f_point} -l {params.l_point} -p {params.points} -a\
        {params.amp} -m {params.m_reads} -c {params.corr} -o {params.output_dir}
        """


rule create_coverage_heatmap:
    input:
        os.path.join(config["run_dir"], "OUTPUT_FILES", "coverage_table.txt")
    output:
        os.path.join(config["run_dir"], "OUTPUT_FILES", "heatmap_coverage.png")
    message:
        "Create a coverage heatmap"
    params:
        script_path = os.path.join(config["scripts_dir"],"heatmap_coverage.py"),
        output_dir = os.path.join(config["run_dir"], "OUTPUT_FILES"),
        width=config["heatmap_width"],
        height=config["heatmap_height"]
    shell:
        """
        python3 {params.script_path} -i {input} -o {params.output_dir} -w {params.width} -e {params.height}
        """