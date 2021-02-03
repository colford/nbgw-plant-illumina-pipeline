# NBGW Plant Illumina Pipeline

This is a pipeline for processing of plant DNA data from an Illumina run as employed by the National Botanic Garden of Wales. It is based upon Dan Smith's (Aberystwyth University) original plant pipeline.

# Overview
The plant illumina pipeline is a Python script that aids in processing Illumina sequences. It has been successfully used to process DNA barcode *rbcL* and ITS2 sequences that have been sequenced using an Illumina machine. The script is currently specialised to run on a SLURM base HPC, specifically the Wales HPC as the SLURM resources and queue names are hardcoded.

The pipeline has the following flow:

 1. Copy \*.fastq.gz files in to a directory.
 2. Run bin/plant-pipeline.py
 3. Set work directory
 4. Set project quality score
 5. Run adaptor check to look for contamination
 6. Run fastq validator to valid fastq files
 7. Run fastqc to get a measure of quality
 8. Run trim and pair to get good quality paired reads
 9. Run merge on trim and pair output  
10. (a) Demultiplex the merged output by primer, length select, convert to fasta and dereplicate the sequences within each sample
11. (b) Concatenate all fasta files and cluster sequences and remove singletons.
12. (c) BLAST the sequences and summarise results

The final output creates an excel file which can be used to manually check the identifications. When the IDs are confirmed, the results are summarised using *create_summary_matrix.py* to output a matrix of samples versus taxa.
 
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
 8. VSERACH

# How to run the Plant Pipeline
First copy the *.fastq.gz files that were output by the illumina processing in to a directory. Run the plant-pipeline.py program and point it at the directory that contains the *.fastq.gz files. The program is menu driven and takes you through each step. The script will setup the directory with the *.fastq.gz files in the following manner:
 
```
  <working-directory>
     original-fastq/                           <- the script will copy all the *.fastq.gz files into here
     slrum-files/                              <- the script will output slrum files to be run by the user
     original-qc/                              <- QC output on original files
         adaptor_check.txt                     <- Adaptor check output
         fastq_validator_check.log             <- Fastq validator output
         fastq/                                <- Fastqc output directory, zip, html, sub-dirs for each zip
                                                  contaning histograms, reports html and summary
     trim-paired-fastq/
         paired/                               <- Sequences that survived pair and trim (go forward to merge)
         unpaired/                             <- Sequences that didn't survive pair and trim
         paired-merged/
             merged/                           <- Sequences that could be merged (go forward to demultiplex)
             not-merged/                       <- Sequences that could not be merged
     demultiplexed/
         unknown/                              <- Sequences which couldn't be identified by primer
             fastq/
         PRIMER/                               
             fastq/                            <- Demultiplexed sequences by primer
             length-selected-fastq/            <- Short sequences removed
             length-selected-fasta/            <- Converted to fasta
             dereplicated-fasta/               <- Sequences are dereplicated within a sample (go forward to cluster)
     clustered/                                <- Dereplicated fastas are concatenated and sequences clustered across samples
     blasted/                                  <- Clustered sequences are blasted against reference database
     manual-checking/                          <- Summarised format of blast results ready for manual checking

```

At each stage the pipeline will create SLURM files and place them in the slrum-files directory. These should be run manually by the user using the appropriate SLURM commands e.g. ```sbatch <file>.slrum```. The running SLURM processes will first touch a "running" file to indicate it has started e.g. "adaptor_check_running" and when the process has finished it will ```mv``` the running file to a "done" file e.g. "adaptor_check_done". The running process will also redirect stdout and stderr to the process name plus their job number e.g. "adaptor_check_964605.out" and "adaptor_check_964605.err". Some of the processes output to stdout and some to stderr so it's best to check both.

## Final output with *create_summary_matrix.py*

After manual checking of the summary blast output, the *-centroid-ids.csv and *-now-manual-edit.xlsx file can be converted into a matrix summarising the total number of sequence reads found in each sample for each identified taxon.

## Helper scripts
There are two helper scripts that extract data out of some of the output files in to a summary table for easy digestion. These are outlined below.

### *process_trim_stats.py*
This script works upon trim_and_pair_original_files_<jobid>.err file. Given the file as input it will output a CSV file containing the sample, total number of reads, successful trim and paired number, forwards only, reverse only and dropped.

### *process_merge_stats.py*
This script works upon merge_trim_and_paired_files_<jobid>.out file. Given the file as input it will output a CSV file containing the sample, number of reads, combined, uncombined and percent combined. 

### Pipeline flow
<img src="https://github.com/colford/nbgw-plant-illumina-pipeline/blob/master/images/barcode-pipeline-flow.png" width="400" height="600" />

# Papers that the pipeline was used in
1. [Using DNA metabarcoding to investigate honey bee foraging reveals limited flower use despite high floral availability](http://www.nature.com/articles/srep42838)
2. [Using DNA Metabarcoding to Identify the Floral Composition of Honey: A New Tool for Investigating Honey Bee Foraging Preferences](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0134735)
3. [Temperate airborne grass pollen defined by spatio-temporal shifts in community composition](https://www.nature.com/articles/s41559-019-0849-7?platform=hootsuite)
4. [Generalisation and specialisation in hoverfly (Syrphidae) grassland pollen transport networks revealed by DNA metabarcoding](https://doi.org/10.1111/1365-2656.12828)
5. [Pollen metabarcoding reveals broad and species-specific resource use by urban bees](https://doi.org/10.7717/peerj.5999)
6. [Floral resource partitioning by individuals within generalised hoverfy pollination networks revealed by DNA metabarcoding](https://doi.org/10.1038/s41598-018-23103-0)
7. [Shifts in honeybee foraging reveal historical changes in floral resources](https://www.nature.com/articles/s42003-020-01562-4)
