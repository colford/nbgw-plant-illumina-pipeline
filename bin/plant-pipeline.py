#!/usr/bin/python
# -*- coding: utf-8 -*-
################################################################################
# Plant illumina pipleline based on Dan Smith's (das50@aber.ac.uk) original
# plant pipeline. This pipeline is targeted for use in a SLRUM HPC system.
#
# Author  : Col R Ford
# Version : V1.00.00
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

import os
import shutil
import fnmatch
import sys
import ConfigParser
import glob
import subprocess
import merged_demultiplex_by_primer as demux_script

################################################################################
# Allow access to the configuration file
################################################################################
class Configuration:
    """Session configuration"""
    config = ConfigParser.ConfigParser()

    def __init__(self):
        self.config.add_section('general')
        self.config.read(self.name())

    def name(self):
        return ('%s/.plant-pipeline.ini' % os.path.expanduser('~'))

    def section(self,section):
        dict1 = {}

        try:
            options = self.config.options(section)
        except:
            return None

        for option in options:
            try:
                dict1[option] = self.config.get(section, option)
                if dict1[option] == -1:
                    DebugPrint("skip: %s" % option)
            except:
                print("exception on %s!" % option)
                dict1[option] = None
        return dict1

    def add_section(self,section):
        self.config.add_section(section)

    def section_add(self,section,key,value):
        self.config.set(section,key,value)

    def __del__(self):
        fd = open(self.name(),'w')
        self.config.write(fd)
        fd.close()


################################################################################
################################################################################
def script_dir():
    return os.path.dirname(os.path.realpath(__file__))

################################################################################
################################################################################
def trimmomatic_jar():
    return '/app/genomics/Trimmomatic/0.33/Trimmomatic-0.33.jar'

################################################################################
################################################################################
def find_files(directory,pattern):
    return glob.glob(os.path.join(directory,pattern))

################################################################################
################################################################################
def fastq_files_recursive(start_dir='.'):
    matches = []
    for root, dirnames, filenames in os.walk(start_dir):
        for filename in fnmatch.filter(filenames, '*.fastq.gz'):
            matches.append(os.path.join(root, filename))
    return matches

################################################################################
################################################################################
def create_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)

################################################################################
################################################################################
def move_files(todir,files):
    for fname in files:
        shutil.move(fname,todir)

################################################################################
################################################################################
def remove_dir(files):
    dirs = set()
    for f in files:
        dirs.add(os.path.dirname(f))
    for dirname in dirs:
        os.rmdir(dirname)

################################################################################
################################################################################
def set_state(config,state):
    project = config.section('general')['project_dir']
    if config.section(project) == None:
        config.add_section(project)
    config.section_add(project,'state',state)

################################################################################
################################################################################
def do_set_current_project_dir(config):
    print('Enter directory:'),
    project_dir = raw_input()   
    if os.path.isdir(project_dir):
        config.section_add('general','project_dir',project_dir)

################################################################################
################################################################################
def do_set_quality_score(config):
    project = config.section('general')['project_dir']  
    print('Enter quality score [20]:'),
    qscore  = raw_input()
    if project == '':
        print('Project not yet set - please set project dir first!')
        return
    if qscore == '':
        qscore = '20'
    config.section_add(project,'quality',qscore)

################################################################################
################################################################################
def do_set_adaptor_dir(config):
    print('Enter directory:'),
    new_adaptor_dir = raw_input()   
    if os.path.isdir(new_adaptor_dir):
        config.section_add('general','adaptor_dir',new_adaptor_dir)

################################################################################
################################################################################
def adaptor_dir(config):
    if 'adaptor_dir' in config.section('general'):
        return config.section('general')['adaptor_dir']
    return 'Not set'

################################################################################
################################################################################
def original_files_dir(config):
    project = config.section('general')['project_dir']
    return os.path.join(project, 'original-fastq')

################################################################################
################################################################################
def trim_paired_files_dir(config):
    project = config.section('general')['project_dir']
    return os.path.join(project, 'trim-paired-fastq')

################################################################################
################################################################################
def slrum_dir(config):
    project = config.section('general')['project_dir']
    return os.path.join(project, 'slrum-files')

################################################################################
################################################################################
def demux_dir(config):
    project = config.section('general')['project_dir']
    return '%s/%s' % (project,"demultiplexed")

def demux_fastq_dir(config):
    return "%s/%s" % (demux_dir(config), "fastq")

################################################################################
################################################################################
def running_file(config,filename):
    return  os.path.join(slrum_dir(config), filename)

################################################################################
################################################################################
def original_qc_dir(config):
    project = config.section('general')['project_dir']
    return  os.path.join(project, 'original-qc')

################################################################################
################################################################################
def do_setup_pipeline_in_current_project_dir(config):
    print('')
    project = config.section('general')['project_dir']
    fastq_gz_files = fastq_files_recursive(project)
    # We expect to find the zipped fastq files in here, do a check to make sure
    number_of_fastq_files = len(fastq_gz_files) 
    if number_of_fastq_files == 0:
        print('ERROR: No *.fastq.gz files found in %s' % project)
        return
    print('INFO: Found %d *.fastq.gz files in %s' % (number_of_fastq_files, project))   
    print('INFO: Creating new project structure...')

    # Create the following directories
    create_dir(original_files_dir(config))
    create_dir(slrum_dir(config))
    move_files(original_files_dir(config),fastq_gz_files)
    try:
        remove_dir(fastq_gz_files)
    except Exception:
        pass

    set_state(config,'Ready for processing')

################################################################################
################################################################################
def current_project_state(config):
    try:
        project = config.section('general')['project_dir']
    except Exception:
        return 'Not set-up'
    if config.section(project) is None:
        return 'Not set-up'
    return config.section(project)['state']

################################################################################
################################################################################
def current_project_quality_score(config):
    project = config.section('general')['project_dir']
    if config.section(project) == None:
        return '20'
    try:
        quality = config.section(project)['quality']
        return quality
    except:
        return '20'

################################################################################
################################################################################
def do_end(config):
    sys.exit(0)

################################################################################
################################################################################
def do_run_adaptor_check_on_original_fastq_files(config):
    # Quick check to see if we can find the illumina adaptors
    # make sure everything is in order
    project = config.section('general')['project_dir']
    slrum_adaptor_check = os.path.join(slrum_dir(config), 'adaptor_check.slrum')

    # Make the QC directory for the results of the adpator check
    create_dir(original_qc_dir(config))

    # Create the SLRUM file to run the adptor check
    print('Writing job file: %s' % slrum_adaptor_check)
    fd = open(slrum_adaptor_check, 'w')
    
    fd.write('#!/bin/bash --login\n')
    fd.write('#SBATCH --job-name=adaptor_check\n')
    fd.write('#SBATCH --output=adaptor_check_%j.out\n')
    fd.write('#SBATCH --error=adaptor_check_%j.err\n')
    fd.write('#SBATCH --exclusive\n')
    fd.write('#SBATCH --ntasks=32\n')
    fd.write('#SBATCH --ntasks-per-node=32\n')
    fd.write('#SBATCH --time=0-07:00\n')
    fd.write('#SBATCH --mem-per-cpu=8000\n')
    fd.write('touch %s\n' % running_file(config,'adaptor_check_running'))
    fd.write('# Adaptors to look for\n')

    fd.write('declare -a adaptor=()\n')
    fd.write('# Nextera\n')
    fd.write('adaptor[0]=Nextera_NX1-2:AGATGTGTATAAGAGACAG\n')
    fd.write('adaptor[1]=Nextera_NX1-2rc:CTGTCTCTTATACACATCT\n')
    fd.write('adaptor[2]=Nextera_Trans1:TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG\n')
    fd.write('adaptor[3]=Nextera_Trans1_rc:CTGTCTCTTATACACATCTGACGCTGCCGACGA\n')
    fd.write('adaptor[4]=Nextera_Trans2:GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG\n')
    fd.write('adaptor[5]=Nextera_Trans2_rc:CTGTCTCTTATACACATCTCCGAGCCCACGAGAC\n')
    fd.write('# Barcode flanking sequences (from Illumina documentation)\n')
    fd.write('adaptor[6]=Nextera_TruseqPCRi5:AATGATACGGCGACCACCGAGATCTACAC\n')
    fd.write('adaptor[7]=Nextera_TruseqPCRi5rc:GTGTAGATCTCGGTGGTCGCCGTATCATT\n')
    fd.write('adaptor[8]=Nextera_TruseqPCRi7:CAAGCAGAAGACGGCATACGAGAT\n')
    fd.write('adaptor[9]=Nextera_TruseqPCRi7rc:ATCTCGTATGCCGTCTTCTGCTTG\n')
    fd.write('adaptor[10]=D701–D712_adapters:GATCGGAAGAGCACACGTCTGAACTCCAGTCAC\n')
    fd.write('adaptor[11]=Truseq3_D501–D508_afterIndex:ACACTCTTTCCCTACACGACGCTCTTCCGATCT\n')
    fd.write('adaptor[12]=Truseq3_D501–D508_afterIndex_rc:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT\n')
    fd.write('# Truseq3\n')
    fd.write('adaptor[13]=Truseq3_PE1:TACACTCTTTCCCTACACGACGCTCTTCCGATCT\n')
    fd.write('adaptor[14]=Truseq3_PE1rc_SE-Universal:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA\n')
    fd.write('adaptor[15]=Truseq3_PE2:GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT\n')
    fd.write('adaptor[16]=Truseq3_SE-PE2SE-Universalrc:AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC\n')
    fd.write('# Truseq2\n')
    fd.write('adaptor[17]=Truseq2_PE1:AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT\n')
    fd.write('adaptor[18]=Truseq2_PE2:CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT\n')
    fd.write('adaptor[19]=Truseq2_PE-PCR1:AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT\n')
    fd.write('adaptor[20]=Truseq2_PE-PCR1rc:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT\n')
    fd.write('adaptor[21]=Truseq2_PE-PCR2:CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT\n')
    fd.write('adaptor[22]=Truseq2_PE-PCR2rc:AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG\n')
    fd.write('adaptor[23]=Truseq2_PE-FlowCell1:TTTTTTTTTTAATGATACGGCGACCACCGAGATCTACAC\n')
    fd.write('adaptor[24]=Truseq2_PE-FlowCell2:TTTTTTTTTTCAAGCAGAAGACGGCATACGA\n')
    fd.write('adaptor[25]=TruSeq2_SE:AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG\n')
    fd.write('adaptor[26]=TruSeq2_PE_f:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT\n')
    fd.write('adaptor[27]=TruSeq2_PE_r:AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG\n')

    fd.write('# Code to run the check\n')
    fd.write('echo "Processing: %s" > "$LOGFILE"\n' % original_files_dir(config))
    fd.write('FILES="%s/*.fastq.gz"\n' % original_files_dir(config))
    fd.write('LOGFILE="%s"\n' % os.path.join(original_qc_dir(config), 'adaptor_check.txt'))
    fd.write('echo "Adapter contamination (total/per million records)" >> "$LOGFILE"\n')
    fd.write('for f in $FILES\n')
    fd.write('do\n')
    fd.write('  echo "Processing $f"\n')
    fd.write('  echo "$f" >> "$LOGFILE"\n')
    fd.write('  for adapt in "${adaptor[@]}"\n')
    fd.write('  do\n')
    fd.write("    tmpHead=`echo ${adapt}| cut -d':' -f 1`\n")
    fd.write("    tmpStr=`echo ${adapt}| cut -d':' -f 2`\n")
    fd.write('    tmpCount=$(zcat $f| head -n 3000000 | grep ${tmpStr} -c)\n')
    fd.write('    echo "Sequence: ${tmpHead} ${tmpStr} Count: ${tmpCount}"\n')
    fd.write('    if test $tmpCount -gt 1\n')
    fd.write('    then\n')
    fd.write('       echo "$tmpHead $tmpStr $tmpCount" >> "$LOGFILE"\n')
    fd.write('    fi\n')
    fd.write('   done\n')
    fd.write('done\n')
    fd.write('mv %s %s\n' % (running_file(config,'adaptor_check_running'), running_file(config,'adaptor_check_done')))

    fd.close()

################################################################################
################################################################################
def had_adaptor_check_been_run(config):
    if os.path.isfile(running_file(config,'adaptor_check_done')):
        return 'Yes'
    return 'No'

################################################################################
################################################################################
def do_run_fastqc_validator_on_original_fastq_files(config):
    # Quick check to see if we can find the illumina adaptors
    # make sure everything is in order
    project = config.section('general')['project_dir']
    slrum_file = os.path.join(slrum_dir(config), 'fastq_validator_check.slrum')

    # Make the QC directory for the results 
    create_dir(original_qc_dir(config))

    # Create the SLRUM file to run the adptor check
    print('Writing job file: %s' % slrum_file)
    fd = open(slrum_file, 'w')
    
    fd.write('#!/bin/bash --login\n')
    fd.write('#SBATCH --job-name=fastq_valid_check\n')
    fd.write('#SBATCH --output=fastq_validator_check_%j.out\n')
    fd.write('#SBATCH --error=fastq_validator_check_%j.err\n')
    fd.write('#SBATCH --exclusive\n')
    fd.write('#SBATCH --ntasks=32\n')
    fd.write('#SBATCH --ntasks-per-node=32\n')
    fd.write('#SBATCH --time=0-05:00\n')
    fd.write('#SBATCH --mem-per-cpu=8000\n')
    fd.write('touch %s\n' % running_file(config,'fastq_validator_check_running'))
    fd.write('# Modules\n')
    fd.write('module add FastQValidator\n')
    fd.write('# Defines\n')
    fd.write('LOGFILE="%s"\n' % os.path.join(original_qc_dir(config), 'fastq_validator_check.log'))
    fd.write('# Code to run the validation\n')
    fd.write('echo "Running fastq validator check" > "$LOGFILE"\n')
    for fname in find_files(original_files_dir(config),'*.fastq.gz'):
        fd.write('"fastQValidator" --file "%s" >> "$LOGFILE"\n' % fname)

    fd.write('mv %s %s\n' % (running_file(config,'fastq_validator_check_running'), running_file(config,'fastq_validator_check_done')))
    fd.close()

################################################################################
################################################################################
def had_fastq_been_validated(config):
    if os.path.isfile(running_file(config,'fastq_validator_check_done')):
        return 'Yes'
    return 'No'

################################################################################
################################################################################
def do_run_fastqc_on_original_fastq_files(config):
    # Quick check to see if we can find the illumina adaptors
    # make sure everything is in order
    project = config.section('general')['project_dir']
    slrum_file = os.path.join(slrum_dir(config), 'fastq_original.slrum')

    # Make the QC directory for the results 
    create_dir(original_qc_dir(config))
    qc_dir = os.path.join(original_qc_dir(config),'fastq')
    create_dir(qc_dir)

    # Create the SLRUM file to run the adptor check
    print('Writing job file: %s' % slrum_file)
    fd = open(slrum_file, 'w')
    
    fd.write('#!/bin/bash --login\n')
    fd.write('#SBATCH --job-name=fastq_original_check\n')
    fd.write('#SBATCH --output=fastq_original_check_%j.out\n')
    fd.write('#SBATCH --error=fastq_original_check_%j.err\n')
    fd.write('#SBATCH --exclusive\n')
    fd.write('#SBATCH --ntasks=32\n')
    fd.write('#SBATCH --ntasks-per-node=32\n')
    fd.write('#SBATCH --time=0-05:00\n')
    fd.write('#SBATCH --mem-per-cpu=8000\n')
    fd.write('touch %s\n' % running_file(config,'fastq_original_running'))
    fd.write('# Modules\n')
    fd.write('module add FastQC/0.11.2\n')
    fd.write('# Code to run the validation\n')
    fd.write('echo "Running fastq on original files"\n')
    for fname in find_files(original_files_dir(config),'*.fastq.gz'):
        fd.write('fastqc "%s" -o "%s" --extract\n' % (fname,qc_dir))

    fd.write('mv %s %s\n' % (running_file(config,'fastq_original_running'), running_file(config,'fastq_original_done')))
    fd.close()

################################################################################
################################################################################
def had_fastq_original_been_run(config):
    if os.path.isfile(running_file(config,'fastq_original_done')):
        return 'Yes'
    return 'No'

################################################################################
################################################################################
def do_run_trim_and_pair_original_fastq_files(config):
    # Quick check to see if we can find the illumina adaptors
    # make sure everything is in order
    project = config.section('general')['project_dir']
    quality = current_project_quality_score(config)
    slrum_file = os.path.join(slrum_dir(config), 'trim_and_pair_original_files.slrum')

    # Create the SLRUM file to run the adptor check
    print('Writing job file: %s' % slrum_file)
    fd = open(slrum_file, 'w')
    
    fd.write('#!/bin/bash --login\n')
    fd.write('#SBATCH --job-name=trim_and_pair_original_files\n')
    fd.write('#SBATCH --output=trim_and_pair_original_files_%j.out\n')
    fd.write('#SBATCH --error=trim_and_pair_original_files_%j.err\n')
    fd.write('#SBATCH --exclusive\n')
    fd.write('#SBATCH --ntasks=32\n')
    fd.write('#SBATCH --ntasks-per-node=32\n')
    fd.write('#SBATCH --time=0-05:00\n')
    fd.write('#SBATCH --mem-per-cpu=8000\n')
    fd.write('touch %s\n' % running_file(config,'trim_and_pair_original_files_running'))
    fd.write('# Modules\n')
    fd.write('module add Java\n')
    fd.write('module add Trimmomatic/0.33\n')
    fd.write('# Code to run the validation\n')
    fd.write('echo "Running trim and pair on original files"\n')

    unique = []
    for fname in find_files(original_files_dir(config),'*.fastq.gz'):
        head, tail = os.path.split(fname)
        unique.append('%s_R' % '_'.join(tail.split('_')[0:-1]))
    idens = list(set(unique))
 
    paired_dir = os.path.join( trim_paired_files_dir(config), 'paired' ) 
    single_dir = os.path.join( trim_paired_files_dir(config), 'unpaired' ) 

    create_dir(paired_dir)  
    create_dir(single_dir)  

    for iden in idens:
        fd.write( 'java -Xms64m -Xmx2000m -jar "%s" PE -threads 8 -phred33 '
                  '%s %s %s %s %s %s ILLUMINACLIP:"%s":2:30:10 '
                  'HEADCROP:3 LEADING:20 SLIDINGWINDOW:4:%s MINLEN:200\n' % ( 
                    trimmomatic_jar(),
                    os.path.join( original_files_dir(config), '%s1.fastq.gz' % iden ),
                    os.path.join( original_files_dir(config), '%s2.fastq.gz' % iden ),
                    os.path.join( paired_dir, '%s1_001P.fastq.gz' % iden ), 
                    os.path.join( single_dir, '%s1_001U.fastq.gz' % iden ), 
                    os.path.join( paired_dir, '%s2_001P.fastq.gz' % iden ), 
                    os.path.join( single_dir, '%s2_001U.fastq.gz' % iden ),
                    os.path.join( adaptor_dir(config), 'NexteraPE-PE.fa' ),
                    quality
                ))
        fd.write('\n')

    fd.write('mv %s %s\n' % (running_file(config,'trim_and_pair_original_files_running'), running_file(config,'trim_and_pair_original_files_done')))
    fd.close()

################################################################################
################################################################################
def do_run_merge_on_trim_and_paired_files(config):
    # Quick check to see if we can find the illumina adaptors
    # make sure everything is in order
    project = config.section('general')['project_dir']
    slrum_file = os.path.join(slrum_dir(config), 'merge_trim_and_paired_files.slrum')

    # Create the SLRUM file to run the adptor check
    print('Writing job file: %s' % slrum_file)
    fd = open(slrum_file, 'w')

    fd.write('#!/bin/bash --login\n')
    fd.write('#SBATCH --job-name=merge_trim_and_paired_files\n')
    fd.write('#SBATCH --output=merge_trim_and_paired_files_%j.out\n')
    fd.write('#SBATCH --error=merge_trim_and_paired_files_%j.err\n')
    fd.write('#SBATCH --exclusive\n')
    fd.write('#SBATCH --ntasks=32\n')
    fd.write('#SBATCH --ntasks-per-node=32\n')
    fd.write('#SBATCH --time=0-05:00\n')
    fd.write('#SBATCH --mem-per-cpu=8000\n')
    fd.write('touch %s\n' % running_file(config,'merge_trim_and_paired_files_running'))
    fd.write('# Modules\n')
    fd.write('module add FLASH/1.2.11\n')
    fd.write('# Code to run the procedure\n')
    fd.write('echo "Running merge on the trim and paired files"\n')

    merged_dir = os.path.join( trim_paired_files_dir(config), 'paired-merged' ) 
    merged_dir_qc = os.path.join( merged_dir, 'qc' )
    merged_dir_merged = os.path.join( merged_dir, 'merged' )
    merged_dir_not_merged = os.path.join( merged_dir, 'not-merged' )
    paired_dir = os.path.join( trim_paired_files_dir(config), 'paired' ) 

    unique = []
    for fname in find_files(paired_dir,'*.fastq.gz'):
        head, tail = os.path.split(fname)
        unique.append('%s_R' % '_'.join(tail.split('_')[0:-2]))
    idens = list(set(unique))

    create_dir(merged_dir)
    create_dir(merged_dir_merged)
    create_dir(merged_dir_qc)
    create_dir(merged_dir_not_merged)

    for iden in idens:
        fd.write( 'flash -d "%s" --max-overlap 450 --min-overlap 10 --max-mismatch-density 0.25 -t 1 -z %s %s\n' % (
                    merged_dir,
                    os.path.join( paired_dir, '%s1_001P.fastq.gz' % iden ),
                    os.path.join( paired_dir, '%s2_001P.fastq.gz' % iden ),
        ))
        fd.write('mv %s %s\n' % ( os.path.join( merged_dir, 'out.notCombined_1.fastq.gz' ), os.path.join( merged_dir_not_merged, '%s1_001P.fastq.gz' % iden )))
        fd.write('mv %s %s\n' % ( os.path.join( merged_dir, 'out.notCombined_2.fastq.gz' ), os.path.join( merged_dir_not_merged, '%s2_001P.fastq.gz' % iden )))
        fd.write('mv %s %s\n' % ( os.path.join( merged_dir, 'out.extendedFrags.fastq.gz' ), os.path.join( merged_dir_merged, '%s_merged.fastq.gz' % iden )))
        fd.write('mv %s %s\n' % ( os.path.join( merged_dir, 'out.hist' ), os.path.join( merged_dir_qc, '%s_hist.txt' % iden )))
        fd.write('mv %s %s\n' % ( os.path.join( merged_dir, 'out.histogram' ), os.path.join( merged_dir_qc, '%s_histogram.txt' % iden ))) 
        fd.write('\n')

    demultiplex_dir = demux_dir(config)
    
    create_dir(demultiplex_dir)

    fd.write('python %s/merged_demultiplex_by_primer.py --merge_dir=%s --output_dir=%s\n' 
                % (script_dir(), merged_dir_merged, demultiplex_dir))

    fd.write('mv %s %s\n' % (running_file(config,'merge_trim_and_paired_files_running'), 
                running_file(config,'merge_trim_and_paired_files_done')))

    fd.close()

################################################################################
################################################################################
def do_demux_length_selection(config):
    # Quick check to see if we can find the illumina adaptors
    # make sure everything is in order
    project = config.section('general')['project_dir']
    slrum_file = os.path.join(slrum_dir(config), 'demux_length_selection.slrum')

    # Create the SLRUM file to run the adptor check
    print('Writing job file: %s' % slrum_file)
    fd = open(slrum_file, 'w')

    fd.write('#!/bin/bash --login\n')
    fd.write('#SBATCH --job-name=demux_length_selected\n')
    fd.write('#SBATCH --output=demux_length_selected_%j.out\n')
    fd.write('#SBATCH --error=demux_length_selected_%j.err\n')
    fd.write('#SBATCH --exclusive\n')
    fd.write('#SBATCH --ntasks=32\n')
    fd.write('#SBATCH --ntasks-per-node=32\n')
    fd.write('#SBATCH --time=0-05:00\n')
    fd.write('#SBATCH --mem-per-cpu=8000\n')
    fd.write('touch %s\n' % running_file(config,'demux_length_selected_running'))
    fd.write('# Modules\n')
    fd.write('module add Java\n')
    fd.write('module add Trimmomatic/0.33\n')
    fd.write("module add vsearch/2.3.2\n")
    fd.write('# Code to run the procedure\n')
    fd.write('echo "Running selection on demuxed files"\n')

    # Get list of directories under demux dir
    # These will be the primers we are using but we need to ignore the "unknown" dir
    demux_d = demux_dir(config)

    for primer in demux_script.start_primers:
        primer_dir = '%s/%s' % (demux_d, primer)

        if not os.path.isdir(primer_dir):
            continue

        primer_fastq_dir = '%s/fastq' % (primer_dir)
        fastq_files = find_files(primer_fastq_dir,'*.fastq.gz')

        if len(fastq_files) == 0:
            continue

        length_selected_fastq_dir = '%s/%s' % (primer_dir, 'length-selected-fastq')
        create_dir(length_selected_fastq_dir)

        length_selected_fasta_dir = '%s/%s' % (primer_dir, 'length-selected-fasta')
        create_dir(length_selected_fasta_dir)

        dereplicated_fasta_dir = '%s/%s' % (primer_dir, 'dereplicated-fasta')
        create_dir(dereplicated_fasta_dir)

        min_seq_length = demux_script.minimum_seq_lengths[primer]

        for fname in fastq_files:
            head, tail = os.path.split(fname)
            fd.write( 'java -Xms64m -Xmx2000m -jar "%s" SE -threads 8 -phred33 %s %s MINLEN:%d\n' % (
                        trimmomatic_jar(),
                        fname,
                        os.path.join( length_selected_fastq_dir, tail ),
                        min_seq_length))
            fd.write('\n')
        
        fd.write('python %s/fastq_to_fasta.py --input_dir=%s --output_dir=%s\n' % 
                    (script_dir(), length_selected_fastq_dir, length_selected_fasta_dir))

        fd.write('python %s/fasta_dereplication.py --input_dir=%s --output_dir=%s\n' % 
                    (script_dir(), length_selected_fasta_dir, dereplicated_fasta_dir))

    fd.write('mv %s %s\n' % (running_file(config,'demux_length_selected_running'), 
                running_file(config,'demux_length_selected_done')))
    fd.close()

################################################################################
################################################################################
def do_concat_and_cluster(config):
    project = config.section('general')['project_dir']
    
    demux_d = demux_dir(config)

    for primer in demux_script.start_primers:
        primer_dir = '%s/%s' % (demux_d, primer)

        if not os.path.isdir(primer_dir):
            continue

        primer_fasta_dir = '%s/dereplicated-fasta' % (primer_dir)
        fasta_files = find_files(primer_fasta_dir,'*.fasta')

        if len(fasta_files) == 0:
            continue

        # Produce the slurm files for each primer
        subprocess.check_output(['python',
                                 '%s/vsearch-concat-and-cluster.py' % script_dir(),
                                 '--input_dir=%s' % primer_fasta_dir,
                                 '--project_dir=%s' % project])


################################################################################
################################################################################
def do_blast(config):
    project = config.section('general')['project_dir']
    
    demux_d = demux_dir(config)

    for primer in demux_script.start_primers:
        primer_dir = '%s/%s' % (demux_d, primer)

        if not os.path.isdir(primer_dir):
            continue

        primer_fasta_dir = '%s/dereplicated-fasta' % (primer_dir)
        fasta_files = find_files(primer_fasta_dir,'*.fasta')

        if len(fasta_files) == 0:
            continue

        # Produce the slurm files for each primer
        subprocess.check_output(['python',
                                 '%s/blast.py' % script_dir(),
                                 '--primer=%s' % primer,
                                 '--project_dir=%s' % project])


################################################################################
################################################################################
def had_trim_and_pair_on_original_been_run(config):
    if os.path.isfile(running_file(config,'trim_and_pair_original_files_done')):
        return 'Yes'
    return 'No'

################################################################################
################################################################################
def had_merge_on_trim_and_paired(config):
    if os.path.isfile(running_file(config,'merge_trim_and_paired_files_done')):
        return 'Yes'
    return 'No'

################################################################################
################################################################################
def had_merged_files_been_length_selected(config):
    if os.path.isfile(running_file(config,'merged_files_length_selected_done')):
        return 'Yes'
    return 'No'

################################################################################
################################################################################
def had_selected_merged_files_been_converted_to_fasta(config):
    if os.path.isfile(running_file(config,'selected_merged_files_convert_and_calapse_done')):
        return 'Yes'
    return 'No'

################################################################################
################################################################################
def top_menu(config):
    key_dispatch =  { 'e': do_end, 
                      'E': do_end,
                      '1': do_set_current_project_dir,
                      '2': do_set_adaptor_dir,
                      '3': do_setup_pipeline_in_current_project_dir,
                      '4': do_set_quality_score,
                      '5': do_run_adaptor_check_on_original_fastq_files,
                      '6': do_run_fastqc_validator_on_original_fastq_files,
                      '7': do_run_fastqc_on_original_fastq_files,
                      '8': do_run_trim_and_pair_original_fastq_files,
                      '9': do_run_merge_on_trim_and_paired_files,
                      'a': do_demux_length_selection,
                      'b': do_concat_and_cluster,
                      'c': do_blast }

    while True:
        print('')
        print('Plant Illumina Pipeline')
        print('-----------------------')
        if 'project_dir' in config.section('general'):
            print('Current project dir                                             : %s' % config.section('general')['project_dir'])
            print('Current adaptor dir                                             : %s' % adaptor_dir(config))
            print('Project dir state                                               : %s' % current_project_state(config))
            print('Project quality score                                           : %s' % current_project_quality_score(config))
            print('Has adaptor check been run?                                     : %s' % had_adaptor_check_been_run(config))
            print('Has fastq been validated?                                       : %s' % had_fastq_been_validated(config))
            print('Has fastq been run on original?                                 : %s' % had_fastq_original_been_run(config))
            print('Has trim and pair been run on original?                         : %s' % had_trim_and_pair_on_original_been_run(config))
            print('Has merge been run on trim and paired files?                    : %s' % had_merge_on_trim_and_paired(config))
            print('Has merged files been length selected?                          : %s' % had_merged_files_been_length_selected(config))
            print('Has selected merged files been converted and collapse to fasta? : %s' % had_selected_merged_files_been_converted_to_fasta(config))
        print('')
        print('1. Set current project directory')
        print('2. Set current adaptor directory')
        if current_project_state(config) == 'Not set-up':
            print('3. Set-up pipeline in current project directory')
        else:
            print('4. Set project quality score')
            print('5. Run adaptor check on original fastq files')
            print('6. Run fastq validator on original fastq files')
            print('7. Run fastqc on original fastq files')
            print('8. Run trim and pair on original fastq files')
            print('9. Run merge and demultiplex on trim and paired fastq files from (8)')
            print('a. Run length selection & dereplicate sequences within a sample (9)')
            print('b. Run concat and cluster from (a)')
            print('c. Run BLAST on clustered data and summarise (b)')
        print('E. Exit')
        print('')
        print('Enter option:'),
        option = raw_input()

        if option in key_dispatch:
            key_dispatch[option](config)

################################################################################
# Main
################################################################################

# Read the current configuration
config = Configuration()

# Process the user selection
top_menu(config)
