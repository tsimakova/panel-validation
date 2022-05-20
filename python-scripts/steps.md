Step 1: Search for undercovered amplicons (input: TXT file with the results of a coverage analysis (VariFind or Oncoscope)) -> TXT file containing tab-separated table with undercovered amplicons + linear regression plot (script: **amplicon_coverage.py**)

Step 2: Run samtools to count the number of mapped reads in a BAM file (input: BAM file) -> number of mapped reads (shell: **samtools**)

Step 3: Count the percentage of reads (input: first and last point of number of reads per amplicon, number of points, number of amplicons in a panel) -> JSON file with the subsampling parameters (script: **subsampling_params.py**)

Step 4: Run samtools for subsampling (input: BAM file, JSON file) -> several BAM files (shell: **samtools)**

Step 5: Run sequtils (input: files after BAM subsampling, BED file with target regions) -> several BED files (shell: **sequtils.jar**)

Step 6: Count LQRs (input: BED files after sequtils) -> several BED files (script: **LQR_counter.py**)

Step 7: Count the proportion of LQRs for each point (input: BED files after LQR_counter.py) -> TXT file with proportion of LQRs for each point (shell: **awk**)

Step 8: Create a LQR proportion lineplot (input: TXT file with proportions of LQR) -> PNG file with LQR lineplot (script: **LQR_proportion_plot.py**)

Step 9: Create a coverage table (input: first and last point of number of reads per amplicon, number of points, number of amplicons in a panel) -> TXT file containing tab-separated coverage table (script: **coverage_table.py**)

Step 10: Create a heatmap for previous step (input: TXT file with coverage table) -> PNG file with coverage heatmap (script: **heatmap_coverage.py**)

