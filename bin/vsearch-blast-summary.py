#!/usr/bin/python
# -*- coding: utf-8 -*-
###########################################################################
# Given the output of the vsearch-pipe after the BLAST this will create
# a SLRUM file to parse the BLAST data so as to create a single entry for
# each sequence. Based upon the BLAST data it will try and classify the
# sequence to species or genus or family. If it can't it will classify
# the sequence as Various.
#
# Command line arguments:
#
# --blastinput <blast input file i.e. output of vsearch-pipe>
# --otuinput   <OTU input file "results" file>
# --outsummary <Output summary CSV file>
#
###########################################################################

#
# Imports
import os
import re
import glob
import csv
from collections import defaultdict
import argparse

# Regular expression to allow us to "mark" results that come from
# our created UK barcodes
uk_re = re.compile('NMW|NBGW|RBGE')

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
# Returns the path to the CSV results file
def mkcsv(dir,file):
    return '%s/%s.csv' % (dir,os.path.splitext(os.path.basename(file))[0])

#
# Returns if this record is an NCBI record or a FPUK record
def is_ncbi(rec):
    ncbi = rec[5].rsplit('|') 
    return len(ncbi) > 2

#
# Gets the list of species and bit scores
def get_species(matched):
    species_list = []
    for match in matched:
        bits = match[3].strip()
        perident = match[2].strip()
        quercov = match[6].strip()
        parts = match[4].split('|')
        species_list.append('%s %s %s (%s; PID:%s%%; QC:%s%%)' % (parts[3],parts[1],parts[2],bits,perident,quercov))

    return species_list

#
# Extracts the species from the matched list
def extract_species(matched):
    species_list = []

    for rows in matched:
        #print(rows)
        parts = rows[4].split('|')

        genus = parts[1]
        species = parts[2]
        if len(parts) < 4:
            print('No family data for this line?', rows)
            family = ''
        else:
            family = parts[3]
        species_list.append('%s %s %s' % (family,genus,species))

    unique = list(set(species_list))
    species_list = []

    for spp in unique:
        tmp = spp.split(' ')
        species_list.append([tmp[0],tmp[1],tmp[2]])

    return species_list

#
# Look through the genus in the list and if one of them
# has a grater influence of 60% then pick that one
def genus_percentage(spp):
    total = len(spp)
    genus = defaultdict(list)

    for sp in spp:
        gen = sp[1]
        if not gen in genus:
            genus[gen] = 1
        else:
            genus[gen] += 1

    for gen in genus:
        if (float(genus[gen])/float(total))*100.0 >= 60:
            return '%s %%' % gen

        families = []
        for sp in spp:
                families.append(sp[0])

        if len(list(set(families))) == 1:
                return families[0]

    return 'Various'

#
# Returns the match type for the species it can be:
# Zero
# Species
# Genus
# Various
def get_match_type(species):
    if len(species) == 0:
        return '----'
    if len(species) == 1:
        return '%s %s' % (species[0][1],species[0][2])

    genus = species[0][1]
    for item in species:
        if genus != item[1]:
            return genus_percentage(species)

    return genus

#
# Outputs the header
def header(fd):
    fd.write('SID,Number-Of-Sequences,Score,Match,Top Species\n')

#
# Output a row in the file
def row(fd,sid,numberof,score,type,species,blast_results,otu):
    try:
        # ID
        fd.write('%s,' % sid)
        # Number of sequences for this sequence
        fd.write('%d,' % int(numberof))
        #fd.write('%d,' % 1)
        # Top score
        fd.write('%f,' % float(score))
        # Type of match ... speices, genus, various
        fd.write('%s' % type)
    except:
        print(sid)
        bang        
        
    # Top bit score taxa matches
    for item in species:
        fd.write(',%s %s' % (item[0],item[1]))

    # We also output all the BLAST results on the same line with their bit score
    # So if something does not look correct we can take a look at the full results
    top_10 = get_species(blast_results)
    fd.write(',')
    
    for spp in top_10:
        fd.write(',%s' % spp)        

    # And now we output what this is made up of
    fd.write(',')
    for ids in otu:
        fd.write(',%s,%d' % (ids,otu[ids]))

    fd.write('\n')

#
# Process one set of blast results for each ID
def process_set(fd,blast_results,otu):
    top_bit_score = 0.0
    top_set = []

    # Work out the top bit score set
    for blast_id in blast_results:
        
        # Sort the entries by top bit score
        blast_results[blast_id].sort(key = lambda row:float(row[3]),reverse=True)

        # Go through the entries and pick out all the top bit scored ones
        for entry in blast_results[blast_id]:
            if float(entry[3]) >= float(top_bit_score):
                top_set.append(entry)
                top_bit_score = entry[3]
            else:
                break
    
    # Extract the species and the matched type
    # i.e. does the top match to a specific species, genus or various?
    species_list = extract_species(top_set)
    matched = get_match_type(species_list)

    # When our sequences where processed they where given the tags
    # id-number, where number was the number of sequences within the 
    # huge data-set that matched exactly with this sequence. i.e. the
    # sequences where merged. So we extract the number for the summaries
    numberof = 1
    if "size" in blast_id:
        parts = blast_id.split('=')
        if parts[-1].isdigit():
            numberof = parts[-1]    
    
    # Output a row in the summary file for this ID set
    row(fd,blast_id,numberof,top_bit_score,matched,species_list,blast_results[blast_id],otu)


#
# Reads the OTU file and creates a dictonary to assosiated clusters including
# the number of sequences that make them up.
def parse_otu(fname):
    print('Parsing:',fname)
    data = {}
    with open(fname, 'r') as otu:
        for line in otu:
            # S record is a single 
            # H record is sussumed by the cluster
            # C record is the total
            parts = line.split()
            if parts[0] == 'S':
                ids = parts[8].split(';')
                data[ids[0]] = {ids[0]:int(ids[1].split('=')[1])}
            elif parts[0] == 'H':
                ids1 = parts[8].split(';')
                ids2 = parts[9].split(';')
                data[ids2[0]][ids1[0]] = int(ids1[1].split('=')[1])
            elif parts[0] == 'C':
                pass
            else:
                pass

    print('Finished Parsing:',fname)
    return data

#
# Read the OTU file, then the BLAST file and process the summary
def blast_summary(blast_fname, otu_fname, out_fname):

    # Read the OTU file in to memory
    otu_dict = parse_otu(otu_fname)

    # Process the BLAST results file
    print('Processing: %s' % blast_fname)
    with open(blast_fname, 'r') as in_fd:
        with open(out_fname, 'w') as out_fd:
            records = csv.reader(in_fd,delimiter=",")

            header(out_fd)

            working_set = defaultdict(list)
            last_id = ""

            # Go through each line in the file
            for row in records:
                # Check to see if we have got to the end of a set for
                # a particular ID. There will be a set of results per ID
                if len(working_set) > 0 and last_id != row[0]:
                    # Process this ID's set
                    ids = last_id.split('=')[1].split(';')[0]
                    process_set(out_fd, working_set, otu_dict[ids])
                    working_set.clear()
            
                # Grab the data from the line
                sid = row[0]
                sid_description = row[1]
                percent_score = row[2]
                bit_score = row[11]
                query_cover = row[13] 
                desc = row[14]
                if len(row) > 15:
                    desc = row[14] + '|' + row[15].split('|')[1]
                description = row[1] + '|' + desc

                # Save off the data into the set
                last_id = sid
                working_set[row[0]].append([sid,1,percent_score,bit_score,description,sid_description,query_cover])

            # There might be one left to process
            if len(working_set) > 0:
                ids = last_id.split('=')[1].split(';')[0]
                process_set(out_fd, working_set, otu_dict[ids])
                working_set.clear()


#
# Main
print('')
print('Running')

# Grab the input parameters
parser = argparse.ArgumentParser(
            description="""Reads the CSV BLAST output to create a summary output""")
parser.add_argument(
            '--blastinput',
            dest='blast_input',
            nargs=1,
            required=True,
            help='BLAST input file as output by vsearch pipeline')
parser.add_argument(
            '--otuinput',
            dest='otu_input',
            nargs=1,
            required=True,
            help='OTU input file as generated by the vsearch cluster')
parser.add_argument(
            '--outsummary',
            dest='out_summary',
            nargs=1,
            required=True,
            help='Output summary CSV file')
args = parser.parse_args()

# Create the summaries for each BLAST result file
blast_summary(args.blast_input[0], args.otu_input[0], args.out_summary[0])
print('End')
