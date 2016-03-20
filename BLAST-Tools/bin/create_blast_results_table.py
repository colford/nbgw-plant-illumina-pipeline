# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 12:51:31 2015
@author: cford
"""

import os
import sys
import argparse
import logging
import csv
import glob

# Work dir where this script is
workdir = os.path.dirname(os.path.realpath(__file__))

# Where the csv summaries are
csv_summary_dir = '%s/../blast_summary' % workdir
table_directory = '%s/../blast_summary_table' % workdir

#
# Makes a directory if does not already exist
def mkdir(dname):
    if not os.path.exists(dname):
        os.makedirs(dname)

#
# Given a directory this returns a list of ext files
def ls(dir,ext):
    return glob.glob('%s/*.%s' % (dir,ext))

#
# Read the CSV file and create summarise by match
def summarise_csv_file(fname,master_table):
    logging.debug('Starting: %s' % fname)
    fd = open(fname,'rb')
    reader = csv.reader(fd)
    reader.next()
    fdata = {}
    for row in reader:
        nos = int(row[1])
        match = row[3].strip(' %*')
        if match in fdata: 
            fdata[match] += nos
        else:
            fdata[match] = nos
    fd.close()
    master_table.append({os.path.basename(fname):fdata})
    outfname = '%s/%s' % (table_directory,os.path.basename(fname))
    logging.debug('CSV output file: %s' % outfname)
    fd = open(outfname,'w')
    for spp in fdata:
        fd.write('%s,%d\n' % (spp,fdata[spp]))
    fd.close()
    logging.debug('Number of unique species for file: %d' % len(fdata)) 
    logging.debug('Done: %s' % fname)
    return

#
# Output the summary table of all the files
def output_summary_table(fname,table):
    # Hold the species data infomation
    spp = []
    tags = ['']
    # Get all the unique species
    for dict in table:
        for tag in dict:
            tags.append(tag)
	    for species in dict[tag]:
                spp.append(species) 
    spp = list(set(spp))
    spp.sort()

    fd = open(fname,'wb')
    writer = csv.writer(fd)
    writer.writerow(tags) 
    for species in spp:
        spp_row = [species]
        for dict in table:
            for tag in dict:
                if species in dict[tag]:
                    spp_row.append(dict[tag][species])
                else:
                    spp_row.append('0')
        writer.writerow(spp_row)
    fd.close()

# Main
if __name__ == "__main__":
    # Setup the logging
    logging.basicConfig(
           level=logging.DEBUG,
           format='[%(asctime)s %(levelname)s] %(message)s' )

    # Parse the command line for the file we are instructed to process
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-dir', '-d', nargs=1, help='full directory path to BLAST summary CSV files',
                          required=False, dest='csv_dir')
    parser.add_argument('--input-file', '-i', nargs=1, help='full path to a single BLAST summary CSV file',
                          required=False, dest='csv_file')
    parser.add_argument('--output-file', '-o', nargs=1, help='full path to a summary output table file',
                          required=True, dest='out_file')
    args = parser.parse_args()

    if args.csv_dir == None and args.csv_file == None:
        parser.print_help()
        sys.exit(1)

    if args.csv_file != None:
        csv_files = [args.csv_file[0]]
    	logging.debug('BLAST summary input file: %s' % csv_files[0])
    else:
    	csv_dir = args.csv_dir[0]
    	csv_files = ls(csv_dir,'csv')
    	logging.debug('BLAST summary input directory: %s' % csv_dir)

    csv_out = args.out_file[0]

    # Make the results table directory if not already
    mkdir(table_directory)

    # Run the summary
    master_table = []
    for fname in csv_files:
    	summary = summarise_csv_file(fname, master_table)

    # Output the summary table
    output_summary_table(csv_out,master_table)

    # End of process
    logging.debug('Done')
