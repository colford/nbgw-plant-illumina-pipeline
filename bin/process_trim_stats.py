#!/usr/bin/python

import argparse
import os

# Function to find n'th instance in a string
def findnth(haystack, needle, n):
    parts= haystack.split(needle, n+1)
    if len(parts)<=n+1:
        return -1
    return len(haystack)-len(parts[-1])-len(needle)

# Should be called with the trim file
parser = argparse.ArgumentParser(description='Creates summary of trim statistics')
parser.add_argument('--trimfile', 
            help='Trim output file to process',
            type=argparse.FileType('r'),
            required=True)
args = parser.parse_args()

# Open the csv file
fd = open('trim_stats.csv','w')
fd.write('Sample,Reads,Both,Forward,Reverse,Dropped\n')

# Work our way through the file. 
# First look for the name of the processed pair and then the status
for line in args.trimfile:
    parts = line.split()
    if len(parts) == 0:
        continue
    # Look for the filename
    if parts[0] == 'TrimmomaticPE:':
        if len(parts) > 3:
            fname_part = parts[7]
            head, tail = os.path.split(fname_part)
            sample_name = tail[0:findnth(tail,'_',2)]
    # Look for the stats
    if parts[0] == 'Input':
        reads = parts[3]
        both = parts[6]
        forwards = parts[11]
        reverse = parts[16]
        dropped = parts[19]
        fd.write('%s,%s,%s,%s,%s,%s\n' %(sample_name,reads,both,forwards,reverse,dropped))

fd.close()
