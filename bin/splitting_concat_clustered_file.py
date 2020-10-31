#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 10:52:15 2019

@author: Laura.Jones

splits the vsearch_concat_clustered.fasta file into two, separating out the 
singletons

"""
import os
import argparse
from Bio import SeqIO

def split_cluster_file(path_and_file_fastq):
    pathdir, fname = os.path.split(path_and_file_fastq)
    multis = '%s/%s-multis.fasta' % (pathdir, fname.split('.')[0]) 
    singles = '%s/%s-singles.fasta' % (pathdir, fname.split('.')[0]) 

    with open(multis, 'w') as mfd:
        with open(singles, 'w') as sfd:
            for index, record in enumerate(SeqIO.parse(path_and_file_fastq, 'fasta')):
                centroid = record.description
                sequence = record.seq
                seq_n = centroid.split(';')[1]
                cent_n =  centroid.split(';')[2]
                if seq_n == "seqs=1" and cent_n == "size=1":
                    sfd.write('>%s\n' % centroid)
                    sfd.write('%s\n' % sequence)
                else:
                    mfd.write('>%s\n' % centroid)
                    mfd.write('%s\n' % sequence)
                

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description="""Splits clustered sequences into singles and multis""")
    parser.add_argument(
            '--input_file',
            dest='input_file',
            required=True,
            help='Input file containing clusterd fasta sequences')
    args = parser.parse_args()
    split_cluster_file(args.input_file)

