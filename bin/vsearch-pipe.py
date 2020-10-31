#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import shutil
import glob

def find_files(directory,pattern):
    return glob.glob(os.path.join(directory,pattern))

def create_dir(path):
    try:
        os.makedirs(path)
    except:
        pass

def concat_files(fname, files):
    with open(fname,'wb') as wfd:
        for f in files:
            with open(f,'rb') as fd:
                shutil.copyfileobj(fd, wfd, 1024*1024*10)

fasta_fnames = find_files('.', '*.fasta')
derep_dir = os.path.abspath('./derep')
cluster_dir = os.path.abspath('./cluster')
concat_file = 'vsearch_concat.fasta'

create_dir(derep_dir)
create_dir(cluster_dir)

with open('vsearch.slurm', 'w') as sl:
    sl.write("#!/bin/bash --login\n")
    sl.write("#SBATCH --job-name=vsearch\n")
    sl.write("#SBATCH --output=vsearch_%j.out\n")
    sl.write("#SBATCH --error=vsearch_%j.err\n")
    sl.write("#SBATCH --exclusive\n")
    sl.write("#SBATCH --ntasks=32\n")
    sl.write("#SBATCH --time=0-24:00\n")
    sl.write("#SBATCH --mem-per-cpu=4000\n")
    sl.write("touch /home/colin.ford/Honey2016_17_LiverpoolRun_201804/slrum-files/vsearch\n")

    sl.write("# Modules\n")
    sl.write("module add vsearch/2.3.2\n")
    sl.write('module add compiler/intel/16.0\n')
    sl.write('module add BLAST/2.2.31+\n')

    sl.write("# Code to run the procedure\n")
    sl.write('echo "Running vsearch on the length selected fasta files"\n')

    derep_files = []

    for fname in fasta_fnames:
        full_fname = os.path.abspath(fname)
        tag = os.path.basename(fname).split('_')[0]
        out_fname = os.path.join(derep_dir, tag + ".fasta")
        sl.write('vsearch --derep_fulllength=%s --sizeout --output=%s --relabel=%s-\n' 
                  % (full_fname, out_fname, tag))
        derep_files.append(out_fname)

    concat_fname = os.path.join(cluster_dir, concat_file)

    sl.write('cat ')
    for fname in derep_files:
        sl.write('%s ' % fname)
    sl.write(' > %s\n' % concat_fname)

    parts = os.path.split(concat_fname)
    clustered_fname = os.path.join(parts[0],parts[1].split('.')[0] + '_clustered.fasta')
    clustered_results = concat_fname + ".results"
    sl.write('vsearch --cluster_size=%s --sizein --sizeout --uc=%s --id=1.0 --consout=%s\n' 
             % (concat_fname,clustered_results,clustered_fname))

    sl.write('# Move to the BLAST db directory\n')
    sl.write('pushd .\n')
    sl.write('cd /home/b.bsp801/BLAST-Tools/blast-db\n')
    sl.write('# Run the BLAST\n')

    parts = os.path.split(clustered_fname)
    clustered_blast_fname = os.path.join(parts[0],parts[1].split('.')[0] + '_blast.csv')
    blast_log = os.path.join(parts[0],parts[1].split('.')[0] + '_blast_log.txt')
    
    sl.write("blastn -query '%s' -db 'restrictedblastdb' -out '%s' -outfmt '10 std score stitle' -max_target_seqs 20 -num_threads 16 >& %s\n"
             % (clustered_fname, clustered_blast_fname, blast_log))

    sl.write('popd\n')

    sl.write('module rm BLAST/2.2.31+\n')
    sl.write('module rm compiler/intel/16.0\n')
    sl.write('module add python/2.6.7\n')

    # Process the BLAST results
    parts = os.path.split(clustered_fname)
    blast_summary = os.path.join(parts[0],parts[1].split('.')[0] + '_blast_summary.csv')
    sl.write("python /home/b.bsp801/bin/vsearch-blast-summary.py --blastinput=%s --otuinput=%s --outsummary=%s\n" %
             (clustered_blast_fname, clustered_results, blast_summary))

    # Sort the file in to largest first to smallest
    parts = os.path.split(blast_summary)
    blast_summary_sorted = os.path.join(parts[0],parts[1].split('.')[0] + '_sorted.csv')
    sl.write("sort -t $',' -k 2 -r -g -o %s %s\n" %
             (blast_summary_sorted, blast_summary))
