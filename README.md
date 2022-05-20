# Pipeline for a targeted gene sequencing panel validation

## Requirements

* awk==5.1.0
* cat==8.32
* openjdk==11.0.13
* python==3.9.12
* samtools==1.11
* sequtils
* snakemake==7.3.8
* matplotlib==3.5.1
* seaborn==0.11.2
* pandas==1.4.0
* scikit-learn==1.0.2
* numpy==1.22.3

## Installation 

Install dependencies 

```pip install -r requirements.txt```

## Usage

### 1. Script for searching for under- and overcovered amplicons

This script is looking for under- and overcovered amplicons. It needs the coverage analysis results (VariFind or/and OncoScope). Using linear regression, the script predicts the relative amplicon coverage. Then the ratio of actual and predicted coverage is counting. In the case the ratio is less then 0.5 (set as a parameter), the amplicon is considered to be undercovered.
When running the script you will be requested to select the path to VariFind or/and OncoScope coverage analysis results, output files directory, threshold for linear regression coefficient, threshold ratio for under- and overcovered amplicons, and width and height of received plot.
When the script completed, we received a file with under- and overcovered amplicons, as well as a linear regression plot for our amplicons (undercovered amplicons are marked in red, overcovered - in green).


### Input

* -i, \--input_files: The path to TSV file(s), containing coverage analysis results (VariFind, OncoScope)
* -d, \--output_dir: The path to the output files directory
* -t, \--threshold: Threshold for R2 in linear regression
* -u, \--under_ratio: The ratio of observed to predicted number of reads for undercovered amplicons
* -o, \--over_ratio: The ratio of observed to predicted number of reads for overcovered amplicons
* -w, \--figure_width: The width of the linear regression plot
* -e, \--figure_height: The height of the linear regression plot


### Run script

```commandline
python3  amplicon_coverage.py 
```

### Output

* ```commandline 
  <TSV-file_prefix>_undercovered_amplicons.txt 
  ```
* ```commandline 
  <TSV-file_prefix>_overrcovered_amplicons.txt
  ```
* ```commandline 
  <TSV-file_prefix>_amplicon_coverage_scatterplot.png
  ```

### 2. Script for selecting the subsampling parameters

This script calculates the parameters for BAM subsampling - the percentage of reads that should remain in the initial BAM file after subsampling, and saves the result in a JSON file. When running the script you will be requested to select the range of the number of reads per amplicon (the first and last points), the number of points, the number of amplicons in the panel, and the number of mapped reads in a BAM file. 

### Input

* -f, \--first_point: The first point among numbers of reads per amplicon
* -l, \--last_point: The last point among numbers of reads per amplicon
* -p, \--points: The number of points (numbers of reads per amplicon)
* -a, \--amp_number: The number of amplicons in a panel
* -m, \--mapped_reads: The number of mapped reads in a BAM file
* -o, \--output_dir: The path to output files

### Run script

```commandline
python3  subsampling_params.py
```

### Output

```commandline 
subsampling_params.json 
```


### 3. Script for counting low quality regions (LQRs)

This script extracts low quality regions (LQR) from 'sequtils regions' results and writes these regions into a new BED-file with 'LQR' suffix. Quality value (QV) corresponds to the sequtils quality corrected coverage (QCC) threshold, i.e. 2^n. When running the script you will be requested to select quality value, input directory with BED files and the output folder.

### Input

* -q, \--quality_threshold: The sequencing quality threshold
* -i, \--input_dir: The path to input BED files directory containing the output of sequtils.jar
* -o, \--output_dir: The path to output files directory

### Run script

```commandline
python3 LQR_counting.py
```

### Output

```commandline 
<BAM_file_prefix>_sub<subsampling_index>_LQR.bed 
```


### 4. Script for plotting the percentage of LQRs for each point (number of reads per amplicon)

This script plots the percentage of LQRs for each number of selected points (the number of reads per amplicon). To run, it needs a TXT file that contains the proportion of positions that belong to the region with low sequencing quality for each point. When running the script you will be requested to select the path to input TXT file, output file directory, the range of the number of reads per amplicon (the first and last points), the number of points, and the width and height of received plot. When the script completed, we received a PNG file, containing the plot of LQR proportion in a target region. 


### Input

* -f, \--first_point: The first point among numbers of reads per amplicon
* -l, \--last_point: The last point among numbers of reads per amplicon
* -p, \--points: The number of points (numbers of reads per amplicon)
* -i, \--input_file: The path to input TXT file, proportion of positions that belong to the region with low sequencing quality for each point
* -o, \--output_dir: The path to output files
* -w, \--figure_width: The width of the LQR plot
* -e, \--figure_height: The height of the LQR plot


### Run script

```commandline
python3 LQR_proportion_plot.py
```

### Output

```commandline 
LQR_proportion_plot.png
```


### 5. Script for creating the table for calculating the number of reads per sample

This script creates the table for calculating the proportion of amplicons with the target coverage. Row names correspond to the number of reads per amplicon, column names - to the number of reads per sample (taking into account the correction coefficient, which is set as a parameter). This table allows the user to chose an appropriate number of reads per sample and calculate how many samples could be multiplexing on a sequencing run. When running the script you will be requested to select the output file directory, the range of the number of reads per amplicon (the first and last points), the number of points, the number of amplicons in a panel, and a correction coefficient. When the script completed, we received a TXT file, containing the proportion of amplicons (%) with the target coverage.


* -f, \--first_point: The first point among numbers of reads per amplicon
* -l, \--last_point: The last point among numbers of reads per amplicon
* -p, \--points: The number of points (numbers of reads per amplicon)
* -a, \--amp_number: The number of amplicons in a panel
* -с, \--correction: The correction coefficient for the number of reads per sample
* -o, \--output_dir: The path to output files

### Run script

```commandline
python3 coverage_table.py
```

### Output

```commandline 
coverage_table.txt
```


### 6. Script for visualization of the table for calculating the number of reads per sample

This script creates the heatmap of proportions of amplicons (%) with the target coverage. To run, it needs a TXT file that contains table for calculating the number of reads per sample. When running the script you will be requested to select the path to input TXT file, output file directory, and the width and height of a heatmap. When the script completed, we received a PNG file, containing the plot of LQR proportion in a target region. 


* -i, \--input_file: The path to input TXT file, containing the table for calculating the number of reads per sample
* -o, \--output_dir: The path to output files
* -w, \--figure_width: The width of the heatmap
* -e, \--figure_height: The height of the heatmap

### Run script

```commandline
python3 heatmap_coverage.py
```

### Output

```commandline 
heatmap_coverage.png
```


## Snakemake pipeline

To run the snakemake pipeline, you need to put the BAM file, BED file with target regions and TSV file with coverage analysis results in a working directory and specify the path to this folder in the configuration file. You also need to enter the prefix of the BAM and TSV files, and the name of BED file (without extension), specify the path to sequtils.jar and other params. 

### Pipeline input:

* \--snakefile: The path to snakefile
* \--configfile: The path to YAML configuration file
* \--cores: Specify the maximum number of CPU cores to be used at the same time (enter a number of cores "--cores N" or do not enter anything "--cores") 

### Run pipeline

```commandline
snakemake --snakefile Snakefile --configfile config.yaml --cores
```

### Pipeline output

* ```commandline 
  <TSV-file_prefix>_undercovered_amplicons.txt 
  ```
* ```commandline 
  <TSV-file_prefix>_overrcovered_amplicons.txt
  ```
* ```commandline 
  <TSV-file_prefix>_amplicon_coverage_scatterplot.png
  ```
* ```commandline 
  LQR_proportion_plot.png 
  ```
* ```commandline 
  coverage_table.txt 
  ```
* ```commandline 
  heatmap_coverage.png 
  ```
