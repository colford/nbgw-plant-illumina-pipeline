#!/usr/bin/python
#
# Dereplicates fasta files using vsearch
#

import os
import subprocess
import argparse
import glob

def find_files(directory,pattern):
    return glob.glob(os.path.join(directory,pattern))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description="""Converts fastq files to fasta""")
    parser.add_argument(
            '--input_dir',
            dest='input_dir',
            required=True,
            help='Input directory containing fastq files')
    parser.add_argument(
            '--output_dir',
            dest='output_dir',
            required=True,
            help='Output directory to hold fasta files')
    args = parser.parse_args()

    fasta_files = find_files(args.input_dir, '*.fasta')

    for fasta in fasta_files:
    	tag = os.path.basename(fasta).split('_')[0]
        out_fname = os.path.join(args.output_dir, tag + ".fasta")
        try:
    	    subprocess.check_output(['vsearch', 
    		    					 '--derep_fulllength=%s' % fasta,
    			    				 '--sizeout',
    				    			 '--output=%s' % out_fname,
    					    		 '--relabel=%s-' % tag])
        except:
            print("ERROR: Failed to run vsearch dereplication on %s" % fasta)
