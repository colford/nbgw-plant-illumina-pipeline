#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import shutil
import glob
import argparse

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


def cluster(project_dir, input_dir):
    fasta_fnames = find_files(input_dir, '*.fasta')
    bin_dir = '%s/../bin' % project_dir
    cluster_dir = os.path.abspath('%s/clustered' % project_dir)
    primer = os.path.basename(os.path.split(input_dir)[0])
    clustered_fname = '%s/clustered-%s.fasta' % (cluster_dir, primer)
    concat_file = '%s/concatenated-%s.fasta' % (cluster_dir, primer)
    results_file = '%s/results-%s.txt' % (cluster_dir, primer)
    vsearch_slurm = 'vsearch-%s' % primer
    vsearch_slurm_fname = '%s.slurm' % vsearch_slurm
    vsearch_path_name = '%s/slrum-files/%s' % (project_dir, vsearch_slurm_fname)

    create_dir(cluster_dir)

    with open(vsearch_path_name, 'w') as sl:
        sl.write("#!/bin/bash --login\n")
        sl.write("#SBATCH --job-name=vsearch\n")
        sl.write("#SBATCH --output=%s_%%j.out\n" % vsearch_slurm)
        sl.write("#SBATCH --error=%s_%%j.err\n" % vsearch_slurm)
        sl.write("#SBATCH --exclusive\n")
        sl.write("#SBATCH --ntasks=32\n")
        sl.write("#SBATCH --time=0-24:00\n")
        sl.write("#SBATCH --mem-per-cpu=4000\n")
        sl.write("touch %s_running\n" % vsearch_slurm)

        sl.write("# Modules\n")
        sl.write("module add vsearch/2.3.2\n")

        sl.write("# Code to run the procedure\n")
        sl.write('echo "Running vsearch on the dereplicated fasta files"\n')

        sl.write('cat ')
        for fname in fasta_fnames:
            sl.write('%s ' % fname)
        sl.write(' > %s\n' % concat_file)

        sl.write('vsearch --cluster_size=%s --sizein --sizeout --uc=%s --id=1.0 --consout=%s\n' 
                 % (concat_file,results_file,clustered_fname))

        sl.write('%s/splitting_concat_clustered_file.py --input_file=%s\n' 
                % (bin_dir, clustered_fname))

        sl.write('mv %s_running %s_done\n' % (vsearch_slurm,vsearch_slurm))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description="""Concatenates and clusters dereplicated files""")
    parser.add_argument(
            '--input_dir',
            dest='input_dir',
            required=True,
            help='Input directory containing dereplicated fasta files')
    parser.add_argument(
            '--project_dir',
            dest='project_dir',
            required=True,
            help='Project directory')

    args = parser.parse_args()
    cluster(args.project_dir, args.input_dir)
