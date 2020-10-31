# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 13:14:39 2018

@author: Laura.Jones
"""

from Bio import SeqIO
import os
import fnmatch
import gzip
import argparse

def fastq_gz_files(start_dir='.'):
    matches = []
    for root, dirnames, filenames in os.walk(start_dir):
        for filename in fnmatch.filter(filenames, '*.fastq.gz'):
            matches.append(os.path.join(root, filename))
    return matches

def compressed_fastq_to_fasta(fastq_path, fasta_dir):
   path, fname = os.path.split(fastq_path)
   fasta_path = os.path.join(fasta_dir, '%s.fasta' % fname.split('.')[0])
   with gzip.open(fastq_path, 'rt') as handle:
       with open(fasta_path, 'wt') as fasta_fh:
           SeqIO.write(
                   sequences=list(SeqIO.parse(handle, 'fastq')),
                   handle=fasta_fh,
                   format='fasta')


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

    fastq_files = fastq_gz_files(args.input_dir)
    for path in fastq_files:
        compressed_fastq_to_fasta(path, args.output_dir)
