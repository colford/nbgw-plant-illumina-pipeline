#!/usr/bin/python
# -*- coding: utf-8 -*-
###############################################################################
# Should be run in the directory that you want to process, could be run in the
# top level directory or in the directory you want to extract from.
#
# The script will read the forward (_R1_) and reverse (_R2_) reads and look
# for specific primers so that it can then seperate these output.
###############################################################################

from Bio import SeqIO
import os
import fnmatch
import errno
import gzip
import re


primer_seq = {
    '_R1_': {
        'rbcLaF': '^ATGTCACCACAAACAGAGACTAAAGC',
        'ITS2F': '^ATGCGATACTTGGTGTGAAT'
    },
    '_R2_': {
        'rbcLr506': '^AGGGGACGACCATACTTGTTCA',
        'ITS3R': '^GACGCTTCTCCAGACTACAAT',
        'uniplantR': '^CCCG[A|C|T][C|T]TGA[C|T][C|T]TG[A|G]GGTC[A|G|T]C',
    }
}


def compressed_fastq_to_fasta(fastq_path):
    path, fname = os.path.split(fastq_path)
    fasta_path = os.path.join(path, '%s.fasta' % fname.split('.')[0])
    with gzip.open(fastq_path, "rt") as handle:
        with open(fasta_path, 'wt') as fasta_fh:
            SeqIO.write(
                sequences=list(SeqIO.parse(handle, "fastq")),
                handle=fasta_fh,
                format="fasta")


def fastq_gz_files(start_dir='.'):
    matches = []
    for root, dirnames, filenames in os.walk(start_dir):
        for filename in fnmatch.filter(filenames, '*.fastq.gz'):
            matches.append(os.path.join(root, filename))
    return matches


def mkdirs():
    try:
        os.makedirs('unknown')
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    for direction in primer_seq:
        for primer in primer_seq[direction]:
            try:
                os.makedirs(primer)
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise


def process(fname, direction):
    print('Processing:', fname)
    metrics = {}
    with gzip.open(fname, "rt") as handle:
        metrics['filename'] = fname
        metrics['direction'] = direction
        metrics['seq'] = {}
        metrics['seq']['unknown'] = []
        for primer in primer_seq[direction]:
            metrics['seq'][primer] = []
        for rseq in SeqIO.parse(handle, "fastq"):
            found = False
            for primer in primer_seq[direction]:
                bcode = primer_seq[direction][primer]
                if bcode in str(rseq.seq):
                    metrics['seq'][primer].append(rseq)
                    found = True
            if not found:
                metrics['seq']['unknown'].append(rseq)
        for prime in metrics['seq']:
            newfname = './%s/%s' % (prime, os.path.basename(fname))
            with gzip.open(newfname, 'wt') as fh:
                SeqIO.write(
                    sequences=metrics['seq'][prime],
                    handle=fh,
                    format="fastq")
    return metrics


def num_seq(meta, primer):
    if primer in meta['seq']:
        return len(meta['seq'][primer])
    return 0


def demultiplex(fastq_files):
    with open('meta_data.csv', 'w') as md:
        md.write('Filename,direction,rbcLaF,rbcLr506,ITS2F,ITS3R,uniplantR,unknown\n')
        for fname in fastq_files:
            for direction in primer_seq:
                if direction in fname:
                    meta = process(fname, direction)
                    md.write('%s,%s,%d,%d,%d,%d,%d\n' % (
                        meta['filename'],
                        meta['direction'],
                        num_seq(meta, 'rbcLaF'),
                        num_seq(meta, 'rbcLr506'),
                        num_seq(meta, 'ITS2F'),
                        num_seq(meta, 'ITS3R'),
                        num_seq(meta, 'uniplantR'),
                        num_seq(meta, 'unknown')))


if __name__ == "__main__":
    mkdirs()
    fastq_files = fastq_gz_files()
    #test = ['./Sample_WHS135-NEXT_TSHT-173/WHS135-NEXT_TSHT-173_TAGGTAGG-CCTACCTA_L001_R1_001.fastq.gz']
    demultiplex(fastq_files)
