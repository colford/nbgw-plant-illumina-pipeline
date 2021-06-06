#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import shutil
import glob
import argparse

def script_dir():
    return os.path.dirname(os.path.realpath(__file__))

def create_dir(path):
    try:
        os.makedirs(path)
    except:
        pass

def blast(project_dir, primer):

    cluster_dir = '%s/clustered' % project_dir
    blast_dir = '%s/blasted' % project_dir

    results_file = '%s/results-%s.txt' % (cluster_dir, primer)

    blast_slurm = 'blast-%s' % primer
    blast_slurm_fname = '%s.slurm' % blast_slurm
    blast_path_name = '%s/slrum-files/%s' % (project_dir, blast_slurm_fname)

    clusterd_file = 'clustered-%s-multis' % primer
    clusterd_path = '%s/%s.fasta' % (cluster_dir, clusterd_file)

    create_dir(blast_dir)

    with open(blast_path_name, 'w') as sl:
        sl.write("#!/bin/bash --login\n")
        sl.write("#SBATCH --job-name=vsearch\n")
        sl.write("#SBATCH --output=%s_%%j.out\n" % blast_slurm)
        sl.write("#SBATCH --error=%s_%%j.err\n" % blast_slurm)
        sl.write("#SBATCH --exclusive\n")
        sl.write("#SBATCH --ntasks=32\n")
        sl.write("#SBATCH --time=0-24:00\n")
        sl.write("#SBATCH --mem-per-cpu=4000\n")
        sl.write("touch %s_running\n" % blast_slurm)

        sl.write("# Modules\n")
        sl.write('module add compiler/intel/16.0\n')
        sl.write('module add BLAST/2.2.31+\n')

        sl.write("# Code to run the procedure\n")
        sl.write('echo "Running BLAST on clustered file"\n')

        sl.write('# Move to the BLAST db directory\n')
        sl.write('pushd .\n')
        sl.write('cd /home/b.bsp801/BLAST-Tools/blast-db\n')
        sl.write('# Run the BLAST\n')

        blast_csv = '%s/blast-%s.csv' % (blast_dir, clusterd_file)
        blast_log = '%s/blast-log-%s.txt' % (blast_dir, clusterd_file)
        
        sl.write("blastn -query '%s' -db 'restrictedblastdb' -out '%s' -outfmt '10 std score qcovs stitle' -max_target_seqs 20 -num_threads 16 >& %s\n"
                 % (clusterd_path, blast_csv, blast_log))

        sl.write('popd\n')

        sl.write('module rm BLAST/2.2.31+\n')
        sl.write('module rm compiler/intel/16.0\n')
        sl.write('module add python/2.6.7\n')

        # Process the BLAST results
        blast_summary = '%s/blast-summary-%s.csv' % (blast_dir, clusterd_file)
        sl.write("python %s/vsearch-blast-summary.py --blastinput=%s --otuinput=%s --outsummary=%s\n" %
                 (script_dir(), blast_csv, results_file, blast_summary))

        # Cur header from the top
        blast_summary_without_header = '%s.wh' % blast_summary
        sl.write("sed -e '1,1d' < %s > %s\n" % (blast_summary, blast_summary_without_header))
        sl.write('mv -f %s %s\n' % (blast_summary_without_header, blast_summary)) 

        # Process blast summary for manual excel checking
        sl.write('module unload python/2.6.7\n')
        sl.write('module add python/2.7.9\n')
        sl.write('python %s/process_vsearch_blast_output_for_excel.py --csv=%s --project_dir=%s\n' %
                    (script_dir(), blast_summary, project_dir)) 

        # Sort the file in to largest first to smallest
        #blast_summary_sorted = '%s/blast-sorted-%s.csv' % (blast_dir, clusterd_file)
        #sl.write("sort -t $',' -k 2 -r -g -o %s %s\n" %
        #         (blast_summary_sorted, blast_summary))

        sl.write('mv %s_running %s_done\n' % (blast_slurm,blast_slurm))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description="""BLASTS dereplicated data""")
    parser.add_argument(
            '--project_dir',
            dest='project_dir',
            required=True,
            help='Project directory')
    parser.add_argument(
            '--primer',
            dest='primer',
            required=True,
            help='Primer name')

    args = parser.parse_args()
    blast(args.project_dir, args.primer)
