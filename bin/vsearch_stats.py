#!/usr/bin/python

import numpy as np
import os

singletons = []
multis = []
multis_list = []

totalsingletons = 0
totalimultis = 0
multinumberofclusters = 0
root_dir = '/home/b.bsp801/Abi/All-fastq/trim-paired-fastq/paired-merged/merged/rbcL-length-selected/cluster/'

with open(os.path.join(root_dir, 'rbcl_blast_summary.csv'), 'r') as fd:
    single_fd = open(os.path.join(root_dir, 'vsearch_concat_clustered_blast_summary_singles.csv'), 'w') 
    multi_fd = open(os.path.join(root_dir, 'vsearch_concat_clustered_blast_summary_multis.csv'), 'w') 
    fd.readline()
    for line in fd:
        parts = line.split(',')
        bitscore = float(parts[2])
        totalseq = int(parts[1])
        clusterno = int(parts[0].split(';')[1].split('=')[1])
        if clusterno == 1 and totalseq == 1:
            singletons.append(bitscore)
            totalsingletons += 1
            single_fd.write(line)
        else:
            multis.append(bitscore)
            multis_list.append(totalseq) 
            totalimultis += totalseq
            multinumberofclusters += 1
            multi_fd.write(line)

    single_fd.close()
    multi_fd.close()

singletons = np.array(singletons)
multis = np.array(multis)
multis_list = np.array(multis_list)

print('============')
print('Total reads                  : ', totalsingletons + totalimultis)
print('Total singles                : ', totalsingletons)
print('Total multies                : ', totalimultis)
print('Total clusters               : ', multinumberofclusters)
print('Singles bitscore min,max,mean: ', singletons.min(), singletons.max(), singletons.mean())
print('Multies bitscore min,max,mean: ', multis.min(), multis.max(), multis.mean())
print('Multies list min,max,mean    : ', multis_list.min(), multis_list.max(), multis_list.mean())
