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
parser = argparse.ArgumentParser(description='Creates summary of merge statistics')
parser.add_argument('--mergefile', 
			help='Merge file to process',
			type=argparse.FileType('r'),
			required=True)
args = parser.parse_args()

# Open the csv file
fd = open('merge_stats.csv','w')
fd.write('Sample,Reads,Combined,Uncombined,Percent\n')

# Keep state
state = 'find_input'

# Work our way through the file. 
# First look for the name of the processed pair and then the status
for line in args.mergefile:
	parts = line.split()
	# Look for the look for the input files section
	if state == 'find_input' and len(parts) > 2 and parts[1] == 'Input' and parts[2] == 'files:':
		state = 'grab_sample_type'
		sample_name = ''
		continue

	# Now grab the input file
	if state == 'grab_sample_type':
		sample_name = line.split('/')[-1].split('.')[0]
		sample_name = sample_name[0:findnth(sample_name,'_',2)]
		state = 'find_stats'
		continue

	# Look for the stats
	if state == 'find_stats' and len(parts) > 2 and parts[1] == 'Total' and parts[2] == 'pairs:':
		total_pairs = parts[3]
		continue

	if state == 'find_stats' and len(parts) > 2 and parts[1] == 'Combined' and parts[2] == 'pairs:':
		combined_pairs = parts[3]
		continue

	if state == 'find_stats' and len(parts) > 2 and parts[1] == 'Uncombined' and parts[2] == 'pairs:':
		uncombined_pairs = parts[3]
		continue

	if state == 'find_stats' and len(parts) > 2 and parts[1] == 'Percent' and parts[2] == 'combined:':
		percent_combined = parts[3]
		fd.write('%s,%s,%s,%s,%s\n' % (sample_name,total_pairs,combined_pairs,uncombined_pairs,percent_combined))
		state = 'find_input'
		continue

fd.close()
