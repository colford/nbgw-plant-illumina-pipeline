#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse
import os
import shutil
import glob

def find_files(directory,pattern):
    return glob.glob(os.path.join(directory,pattern))

parser = argparse.ArgumentParser(description='blast the multis cluster fasta file')
parser.add_argument('primer', choices=['its2', 'rbcl'],
            help='specify primer: its2 or rbcl')
primer_name = parser.parse_args()

cluster_file = os.path.abspath('vsearch_concat_clustered.fasta')
cluster_file_multis = os.path.abspath('vsearch_concat_clustered_multis.fasta')
concat_file = os.path.abspath('vsearch_concat.fasta')
clustered_results = os.path.abspath('vsearch_concat.fasta.results')

blast_database = vars(primer_name)['primer'] + "-2019-11-restricted"

with open('blast-concat-clustered-file.slurm', 'w') as sl:
    sl.write("#!/bin/bash --login\n")
    sl.write("#SBATCH --job-name=vsearch\n")
    sl.write("#SBATCH --output=blast_cluster_%j.out\n")
    sl.write("#SBATCH --error=blast_cluster_%j.err\n")
    sl.write("#SBATCH --exclusive\n")
    sl.write("#SBATCH --ntasks=32\n")
    sl.write("#SBATCH --time=0-24:00\n")
    sl.write("#SBATCH --mem-per-cpu=4000\n")

    sl.write("# Modules\n")
    sl.write("module add vsearch/2.3.2\n")
    sl.write('module add compiler/intel/16.0\n')
    sl.write('module add BLAST/2.2.31+\n')

    sl.write('# Move to the BLAST db directory\n')
    sl.write('pushd .\n')
    sl.write('cd /home/b.bsp801/BLAST-Tools/blast-db\n')
    sl.write('# Run the BLAST\n')

    parts = os.path.split(cluster_file_multis)
    clustered_blast_fname = os.path.join(parts[0],parts[1].split('.')[0] + '_blast.csv')
    blast_log = os.path.join(parts[0],parts[1].split('.')[0] + '_blast_log.txt')
    
    sl.write("blastn -query '%s' -db '%s' -out '%s' -outfmt '10 std score qcovs stitle' -max_target_seqs 20 -num_threads 16 >& %s\n"
             % (cluster_file_multis, blast_database, clustered_blast_fname, blast_log))
	sl.write('echo "Finished BLAST"\n')
			 
