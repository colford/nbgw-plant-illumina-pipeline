#####################################################################################
# BLAST sequences against a local BLAST database
#
# We created a local BLAST database from the downloaded NCBI and then filtered the
# database for plant based GI's therefore cutting down our search area. We also
# have a local BLAST database for the extra UK species that we have barcoded but
# have yet to upload to GenBank.
#
# The output of this program is used as input for "blast_summary.py"
#
# Change the defines for your specific requirements.
# We assume that the sequences have already had their primers and tags removed
# and are stored in fasta files with each file being for one tag. 
#####################################################################################

#
# Imports
import os
import glob
from   Bio import SeqIO
from   Bio.Blast.Applications import NcbiblastnCommandline

#
# Defines
# Change these for your specific needs

# Main directory for work
workdir   = os.path.dirname(os.path.realpath(__file__))

# Local BLAST database and GI filter list
blast_db  = '%s/../blast-db' % workdir

# FASTA directory, where to find the sequences
fasta_dir = '%s/../fasta' % workdir

# Output of our BLAST results
outdir    = '%s/../blast_results' % workdir

#
# Given a directory this returns a list of fasta files
# Change if 'fa' is not the extension that you want to find.
def get_fasta_files(fasta_dir):
    abs_path = os.path.abspath(fasta_dir)
    print('Looking for fasta files in: ',abs_path)
    return glob.glob('%s/*.fasta' % abs_path)

#
# Returns the output file
def outfile(dir,file):
    return '%s/%s.csv' % (dir,os.path.splitext(os.path.basename(file))[0])

#
# BLASTS the sequence file against the local database
def ncbi_blast(in_file,out_file,dbp):
    cur_dir = os.getcwd()
    os.chdir(dbp)
	# Replace the database names with your own local databases
	# ...plus we are using 8 threads so change according to the resources available
    cmd_line = NcbiblastnCommandline(query=in_file, db="'nt_ncbi_plants fpuk'", out=out_file, outfmt="'10 std score stitle'", max_target_seqs=20,num_threads=8)
    cmd_line()
    os.chdir(cur_dir)

#
# Each fasta file contains a set of sequences that where
# matched for a given tag. They are either reverse of forward.
# We blast them to the NCBI database.
def blast_sequences( fastas, odir, db_dir ):
    for file in fastas:
        print( 'Processing: ', file )
        ncbi_blast(file,outfile(odir,file),db_dir)            

# Main
print('Running')
os.chdir(workdir)
fasta_files = get_fasta_files(fasta_dir)
blast_sequences(fasta_files,outdir,blast_db)
print('Done')