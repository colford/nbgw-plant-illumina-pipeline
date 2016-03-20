###########################################################################
# Given the output of "blast_with_ncbi.py" we process each tag file that
# contains all the BLAST results and output a CSV file with our best guess
# at what the BLAST result is informing us.
#
# Change the defines for your own specific requirements
###########################################################################

#
# Imports
import os
import re
import glob
import csv
from   collections import defaultdict

#
# Defines

# Work dir where this script is
workdir   = os.path.dirname(os.path.realpath(__file__))

# BLAST results directory that was created by "blast_with_ncbi.py"
file_dir  = '%s/../blast_results' % workdir

# Where to place our summaries
res_dir   = '%s/../blast_summary' % workdir

# Regular expression to allow us to "mark" results that come from
# our created UK barcodes
uk_re     = re.compile('NMW|NBGW|RBGE')

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
# Gets the list of species and bit scores
def get_species(matched):
    species_list = []
    for match in matched:
        bits = match[3].strip()
        species = match[4]
        is_ncbi = species.rsplit('|')
        is_uk   = species.split('-')

        # Test to see if this ncbi
        if len(is_ncbi) > 2:        
            # From NCBI
            str = is_ncbi[-1]
            parts = str.split(' ') # Always a space at the start
            genus = parts[1]
            spec  = parts[2]
            if uk_re.search(str):
                spec += '*'
            species_list.append('%s %s (%s)' % (genus,spec,bits))
        else:
            # From UK
            species_list.append('%s %s* (%s)' % (is_uk[1],is_uk[2],bits))

    return species_list

#
# Extracts the species from the matched list
def extract_species(matched):
    species_list = []
    for match in matched:
        species = match[4]
        is_ncbi = species.rsplit('|')
        is_uk   = species.split('-')

        # Test to see if this ncbi
        if len(is_ncbi) > 2:        
            # From NCBI
            str = is_ncbi[-1]
            parts = str.split(' ') # Always a space at the start
            genus = parts[1]
            spec  = parts[2]
            if uk_re.search(str):
                spec += '*'
            species_list.append('%s %s' % (genus,spec))
        else:
            # From UK we mark this with a *
            species_list.append('%s %s*' % (is_uk[1],is_uk[2]))

    unique = list(set(species_list))
    species_list = []

    for spp in unique:
        tmp = spp.split(' ')
        species_list.append([tmp[0],tmp[1]])

    return species_list

#
# Look through the genus in the list and if one of them
# has a grater influence of 60% then pick that one
def genus_percentage(spp):
    total = len(spp)
    genus = defaultdict(list)

    for sp in spp:
        gen = sp[0]
        if not gen in genus:
            genus[gen] = 1
        else:
            genus[gen] += 1

    for gen in genus:
        if (float(genus[gen])/float(total))*100.0 >= 60:
            return '%s %%' % gen

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
        return '%s %s' % (species[0][0],species[0][1])

    genus = species[0][0]
    for item in species:
        if genus != item[0]:
            return genus_percentage(species)

    if species[0][1][-1] == '*':
        genus = '%s*' % genus

    return genus

#
# Outputs the header
def header(fd):
    fd.write('SID,Number-Of-Sequences,Score,Match,Top Species\n')

#
# Output a row in the file
def row(fd,sid,numberof,score,type,species,blast_results):
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

    fd.write('\n')

#
# Process one set of blast results for each ID
def process_set(fd,blast_results):
    top_bit_score = 0.0
    top_set = []

    # Work out the top bit score set
    for blast_id in blast_results:
        
	   # Sort the entries by top bit score
        blast_results[blast_id].sort(key = lambda row: row[3],reverse=True)
		
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
    if "-" in blast_id:
        parts = blast_id.split('-')
        if parts[0].isdigit():
            numberof = parts[1]    
    
	# Output a row in the summary file for this ID set
    row(fd,blast_id,numberof,top_bit_score,matched,species_list,blast_results[blast_id])

#
# Goes through a list of BLAST CSV results and summaries each file as a CSV file
def blast_summary(dir,files):
    # Create the summary folder
    mkdir(dir)
	
    # Process the BLAST results file by file
    for file in files:
        print('Processing: %s' % file)
        in_fd = open(file,'r')
        out_fd = open(mkcsv(dir,file),'w')
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
                process_set(out_fd, working_set)
                working_set.clear()
        
		 # Grab the data from the line
            sid = row[0]
            percent_score = row[2]
            bit_score = row[11]
            description = row[13]
			
		 # Save off the data into the set
            last_id = sid
            working_set[row[0]].append([sid,1,percent_score,bit_score,description])

        # There might be one left to process
        if len(working_set) > 0:
            process_set(out_fd, working_set)
            working_set.clear()

        in_fd.close()
        out_fd.close()

#
# Main
print('')
print('Running')

# Get the BLAST results
files = ls(file_dir,'csv')
#files = ['/home/colin.ford/BLAST-Tools/blast_results/TREE6_TAAGGCGA-TATCCTCT_L001_merged.csv']

# Create the summaries for each BLAST result file
blast_summary(res_dir,files)
print('Enend')
