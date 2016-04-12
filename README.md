# NBGW Plant Illumina Pipeline
This is a pipeline for processing of plant DNA data from an illumina run as employed by the National Botanic Garden of Wales. It is based upon Dan Smith's (Aberystwyth University) original plant pipeline.

# Overview
The plant illumina pipeline is a Python script that aids in processing illumina sequences. It has been successfully used to process DNA barcode rbcL sequences that have been sequenced using an illumina machine. The script is currently specilised to run on a SLURM base HPC, specifically the Wales HPC as the SLURM resources and queue names are hardcoded at the moment.

The pipeline has the folling flow:

 1. Copy *.fastq.gz files in to a directory.
 2. Run bin/plant-pipeline.py
 3. Set work directory
 4. Run adaptor check to look for contanimation
 5. Run fastq validator to valid fastq files
 6. Run fastqc to get a measure of quality
 7. Run trim and pair to get good quality paired reads
 8. Run merge on trim and pair output
 9. Run length selection on merged output
 10. Conver and collapse to merged fastq to fasta format
 11. Run BLAST on fasta output against a given BLAST database
 
# Dependencies
The script expects the following to be avaliable on the HPC. 

 1. Python
 2. FastQValidator
 2. FastQC/0.11.2
 3. Java
 4. Trimmomatic/0.33
 5. FLASH/1.2.11
 6. fastx_toolkit/0.0.13.2
 7. BLAST+/2.2.31

# How to run the Plant Pipeline
First copy the *.fastq.gz files that were output by the illumina processing in to a directory. Run the plant-pipeline.py program and point it at the directory that contains the *.fastq.gz files. The program is menu driven and takes you through each step. The script will setup the directory with the *.fastq.gz files in the following manner:
 
```
  <working-directory>
     original-fastq/                           <- the script will copy all the *.fastq.gz files in to here
     slrum-files/                              <- the script will output slrum files to be run by the user
     original-qc/                              <- QC output on original files
         adaptor_check.txt                     <- Adaptor check output
         fastq_validator_check.log             <- Fastq validator output
         fastq/                                <- Fastqc output directory, zip, html, sub-dirs for each zip
                                                  contaning histograms, reports html and summary
     trim-paired-fastq/
         paired/                               <- Sequences that survived pair and trim
         unpaired/                             <- Sequences that didn't survive pair and trim
         paired-merged/
             merged/                           <- Sequences that could be merged
             merged-length-selected/           <- Merged sequences selected for length (fastq)
             merged-length-secected-fasta/     <- Merged sequences output as fasta
             merged-length-selected-calapsed/  <- Merged sequences collapsed from fastas 
                                                  i.e. identical sequences combined
             not-merged/                       <- Sequences that could not be merged
             qc/                               <- Quality histograms on merged data
```

At each stage the pipeline will create SLURM files and place them in the slrum-files directory. These should be run manually by the user using the appropriate SLURM commands e.g. sbatch ```<file>.slrum```. The running SLURM processes will first touch a "running" file to indicate it has started e.g. "adaptor_check_running" when the process has finished it will mv the running file to a "done" file e.g. "adaptor_check_done". The running process will also redirect stdout and stderr to the process name plus their job number e.g. "adaptor_check_964605.out" and "adaptor_check_964605.err". Some of the processes output to stdout and some to stderr so it's best to check both.

## Helper scripts
There are two helper scripts that extract data out of some of the output files in to a summary table for easy digestion. These are outlined below.

### process_trim_stats.py
This script works upon trim_and_pair_original_files_<jobid>.err file. Given the file as input it will output a CSV file containing the sample, total number of reads, successful trim and paired number, forwards only, reverse only and dropped.

### process_merge_stats.py
This script works upon merge_trim_and_paired_files_<jobid>.out file. Given the file as input it will output a CSV file containing the sample, number of reads, combined, uncombined and percent combined. 

# Running BLAST on the output
Once you have run through the pipeline and you have your collapsed FASTA files then you might want to run them through BLAST. To run the BLAST follow these instructions. Note that the scripts are mostly hardcoded at the moment and will require editing for your particular paths and needs.

 1. cd BLAST-Tools
 2. Create a directory called "fasta/DNA" and copy the collapsed FASTA output from the pipeline in to it.
 3. Create a directory called "blast-db" and copy your BLAST database in to it. You can download the whole BLAST database from the NCBI or create your own.
 4. Create a directory called "blast_summary" and "blast_results" and "blast_summary_table"
 5. cd jobs
 6. Edit the blast-file-template.txt directory paths and edit the database names in the blastn command line.
 7. Run the make-jobs.sh script. This will create the jobs in the slrum/ directory and start running them.
 8. The csv output for each BLAST file will be output in "blast_results"
 9. Once the BLAST has finished run bin/blast_summary.py to summerise the BLAST results for each file.
 10. Once the BLAST summary has finished run the bin/create_blast_results_table.py to summarise the summaries.
