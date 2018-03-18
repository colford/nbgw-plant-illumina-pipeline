#!/usr/bin/python
# -*- coding: utf-8 -*-
###############################################################################
# Plant illumina pipleline based on Dan Smith's (das50@aber.ac.uk) original
# plant pipeline. This pipeline is targeted for use in a SLRUM HPC system.
#
# Author  : Col R Ford
# Version : V2.00.00
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
###############################################################################

import os
import shutil
import sys
import ConfigParser
import glob


class Configuration:
    """Project configuration class

    Stores the configuration file for the pipeline and the projects. Allows
    one to see what state they are in within a project and also allows some
    parameter settings. The pipeline works on one project at a time so that
    the user can run the pipeline from any directory.

    """

    config = ConfigParser.ConfigParser()

    def __init__(self):
        self.config.add_section('general')
        self.config.read(self.name())

    def name(self):
        return ('%s/.plant-pipeline.ini' % os.path.expanduser('~'))

    def section(self, section):
        dict1 = {}

        try:
            options = self.config.options(section)
        except Exception as e:
            return None

        for option in options:
            try:
                dict1[option] = self.config.get(section, option)
                if dict1[option] == -1:
                    print("skip: %s" % option)
            except Exception as e:
                print("exception on %s!" % option)
                dict1[option] = None
        return dict1

    def add_section(self, section):
        self.config.add_section(section)

    def section_add(self, section, key, value):
        self.config.set(section, key, value)

    def __del__(self):
        fd = open(self.name(), 'w')
        self.config.write(fd)
        fd.close()


def trimmomatic_jar():
    ''' Returns the path to the traimmomatic java program
    '''
    return '/app/genomics/Trimmomatic/0.33/Trimmomatic-0.33.jar'


def find_files(directory, pattern):
    ''' Given a directory and a pattern return the found files
    '''
    return glob.glob(os.path.join(directory, pattern))


def create_dir(path):
    ''' Given a path, if it doesn't exist then create the directory
    '''
    if not os.path.exists(path):
        os.makedirs(path)


def move_files(todir, files):
    ''' Moves the list of files given to the directory given
    '''
    for fname in files:
        shutil.move(fname, todir)


def set_state(config, state):
    ''' Sets the given state of the project i.e. where the project has got
        to within the pipeline. State is stored in the config file.
    '''
    project = config.section('general')['project_dir']
    if config.section(project) is None:
        config.add_section(project)
    config.section_add(project, 'state', state)


def project_dir(config):
    ''' Returns the currently set project directory
    '''
    return config.section('general')['project_dir']


def do_set_current_project_dir(config):
    ''' Allows the user to set the current project directory. This will be
        added to the configuration file.
    '''
    print('Enter directory:'),
    project_dir = input()
    if os.path.isdir(project_dir):
        config.section_add('general', 'project_dir', project_dir)


def do_set_quality_score(config):
    ''' Allows the user to set the quality score to be used in the trim.
        This is stored within the configuration file. Default quality score
        is 20.
    '''
    project = project_dir(config)
    print('Enter quality score [20]:'),
    qscore = input()
    if project == '':
        print('Project not yet set - please set project dir first!')
        return
    if qscore == '':
        qscore = '20'
    config.section_add(project, 'quality', qscore)


def do_set_adaptor_dir(config):
    ''' Allows the user to set the directory where the apdaptor file resides.
        This is stored in the configuration file.
    '''
    print('Enter directory:'),
    new_adaptor_dir = input()
    if os.path.isdir(new_adaptor_dir):
        config.section_add('general', 'adaptor_dir', new_adaptor_dir)


def adaptor_dir(config):
    ''' Returns the adaptor directory or the string "Not set"
    '''
    if 'adaptor_dir' in config.section('general'):
        return config.section('general')['adaptor_dir']
    return 'Not set'


def original_files_dir(config):
    ''' Returns the directory containing the original fastq files. Note that
        these files are never chnanged or manipulated in any form.
    '''
    project = project_dir(config)
    return os.path.join(project, 'original-fastq')


def trim_paired_files_dir(config):
    ''' Returns the directory that contains the trim paired fastq files.
    '''
    project = project_dir(config)
    return os.path.join(project, 'trim-paired-fastq')


def slrum_dir(config):
    ''' Returns the directory that contains the slrum files for running on the
        HPC. Use sbatch to queue up these files.
    '''
    project = project_dir(config)
    return os.path.join(project, 'slrum-files')


def running_file(config, filename):
    ''' Returns the filename to run within the slrum directory
    '''
    return os.path.join(slrum_dir(config), filename)


def original_qc_dir(config):
    ''' Returns the quality check directory for the original fastq files.
    '''
    project = project_dir(config)
    return os.path.join(project, 'original-qc')


def do_setup_pipeline_in_current_project_dir(config):
    ''' This assumes that the project directory contains all the fastq files.
        It will take these files and place them in the directory called
        "original-fastq". It will also create the other directories that are
        required when running the pipeline such as "slrum-files" and it will
        create an entry for the project within the configuration file.
    '''
    print('')
    project = project_dir(config)
    fastq_gz_files = find_files(project, "*.fastq.gz")
    # We expect to find the zipped fastq files in here, do a check to make sure
    number_of_fastq_files = len(fastq_gz_files)
    if number_of_fastq_files == 0:
        print('ERROR: No *.fastq.gz files found in %s' % project)
        return
    print('INFO: Found %d *.fastq.gz files in %s'
          % (number_of_fastq_files, project))
    print('INFO: Creating new project structure...')

    # Create the following directories
    create_dir(original_files_dir(config))
    create_dir(slrum_dir(config))
    move_files(original_files_dir(config), fastq_gz_files)

    set_state(config, 'Ready for processing')


def current_project_state(config):
    ''' Returns the current project state from the configuration file.
    '''
    project = config.section('general')['project_dir']
    if config.section(project) is None:
        return 'Not set-up'
    return config.section(project)['state']


def current_project_quality_score(config):
    ''' Returns the trim quality score for the project. The default is set
        to 20.
    '''
    project = project_dir(config)
    if config.section(project) is None:
        return '20'
    try:
        quality = config.section(project)['quality']
        return quality
    except Exception as e:
        return '20'


def do_end(config):
    ''' Exits the pipeline back to the shell or OS
    '''
    sys.exit(0)


def do_run_adaptor_check_on_original_fastq_files(config):
    ''' Quick check to see if we can find the illumina adaptors make sure
        everything is in order. This needs to be revised so that the
        adaptors are read from the adaptor file instead of being hard
        coded.

        Writes a slurm file out to be run with sbatch.
    '''
    slrum_adaptor_check = os.path.join(slrum_dir(config),
                                       'adaptor_check.slrum')

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
    fd.write('#SBATCH --ntasks=1\n')
    fd.write('#SBATCH --time=0-07:00\n')
    fd.write('#SBATCH --mem-per-cpu=8000\n')
    fd.write('touch %s\n' % running_file(config, 'adaptor_check_running'))
    fd.write('# Adaptors to look for\n')

    fd.write('declare -a adaptor=()\n')
    fd.write('# Nextera\n')
    fd.write('adaptor[0]=Nextera_NX1-2:AGATGTGTATAAGAGACAG\n')
    fd.write('adaptor[1]=Nextera_NX1-2rc:CTGTCTCTTATACACATCT\n')
    fd.write('adaptor[2]=Nextera_Trans1:TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG\n')
    fd.write(
        'adaptor[3]=Nextera_Trans1_rc:CTGTCTCTTATACACATCTGACGCTGCCGACGA\n')
    fd.write('adaptor[4]=Nextera_Trans2:GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG\n')
    fd.write(
        'adaptor[5]=Nextera_Trans2_rc:CTGTCTCTTATACACATCTCCGAGCCCACGAGAC\n')
    fd.write('# Barcode flanking sequences (from Illumina documentation)\n')
    fd.write('adaptor[6]=Nextera_TruseqPCRi5:AATGATACGGCGACCACCGAGATCTACAC\n')
    fd.write(
        'adaptor[7]=Nextera_TruseqPCRi5rc:GTGTAGATCTCGGTGGTCGCCGTATCATT\n')
    fd.write('adaptor[8]=Nextera_TruseqPCRi7:CAAGCAGAAGACGGCATACGAGAT\n')
    fd.write('adaptor[9]=Nextera_TruseqPCRi7rc:ATCTCGTATGCCGTCTTCTGCTTG\n')
    fd.write(
        'adaptor[10]=D701–D712_adapters:GATCGGAAGAGCACACGTCTGAACTCCAGTCAC\n')
    fd.write(
        'adaptor[11]=Truseq3_D501–D508_afterIndex:ACACTCTTTCCCTACACGACG'
        'CTCTTCCGATCT\n')
    fd.write(
        'adaptor[12]=Truseq3_D501–D508_afterIndex_rc:AGATCGGAAGAGCGTCGTG'
        'TAGGGAAAGAGTGT\n')
    fd.write('# Truseq3\n')
    fd.write('adaptor[13]=Truseq3_PE1:TACACTCTTTCCCTACACGACGCTCTTCCGATCT\n')
    fd.write(
        'adaptor[14]=Truseq3_PE1rc_SE-Universal:AGATCGGAAGAGCGTCGTGTAGGGA'
        'AAGAGTGTA\n')
    fd.write('adaptor[15]=Truseq3_PE2:GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT\n')
    fd.write(
        'adaptor[16]=Truseq3_SE-PE2SE-Universalrc:AGATCGGAAGAGCACACGTCTGA'
        'ACTCCAGTCAC\n')
    fd.write('# Truseq2\n')
    fd.write(
        'adaptor[17]=Truseq2_PE1:AATGATACGGCGACCACCGAGATCTACACTCTT'
        'TCCCTACACGACGCTCTTCCGATCT\n')
    fd.write(
        'adaptor[18]=Truseq2_PE2:CAAGCAGAAGACGGCATACGAGATCGGTCTCGGC'
        'ATTCCTGCTGAACCGCTCTTCCGATCT\n')
    fd.write(
        'adaptor[19]=Truseq2_PE-PCR1:AATGATACGGCGACCACCGAGATCTACACTC'
        'TTTCCCTACACGACGCTCTTCCGATCT\n')
    fd.write(
        'adaptor[20]=Truseq2_PE-PCR1rc:AGATCGGAAGAGCGTCGTGTAGGGAAAGAG'
        'TGTAGATCTCGGTGGTCGCCGTATCATT\n')
    fd.write(
        'adaptor[21]=Truseq2_PE-PCR2:CAAGCAGAAGACGGCATACGAGATCGGTCTCG'
        'GCATTCCTGCTGAACCGCTCTTCCGATCT\n')
    fd.write(
        'adaptor[22]=Truseq2_PE-PCR2rc:AGATCGGAAGAGCGGTTCAGCAGGAATGCCG'
        'AGACCGATCTCGTATGCCGTCTTCTGCTTG\n')
    fd.write(
        'adaptor[23]=Truseq2_PE-FlowCell1:TTTTTTTTTTAATGATACGGCGACCACC'
        'GAGATCTACAC\n')
    fd.write(
        'adaptor[24]=Truseq2_PE-FlowCell2:TTTTTTTTTTCAAGCAGAAGACGGCATACGA\n')
    fd.write('adaptor[25]=TruSeq2_SE:AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG\n')
    fd.write('adaptor[26]=TruSeq2_PE_f:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT\n')
    fd.write('adaptor[27]=TruSeq2_PE_r:AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG\n')

    fd.write('# Code to run the check\n')
    fd.write(
        'echo "Processing: %s" > "$LOGFILE"\n' % original_files_dir(config))
    fd.write('FILES="%s/*.fastq.gz"\n' % original_files_dir(config))
    fd.write('LOGFILE="%s"\n' %
             os.path.join(original_qc_dir(config), 'adaptor_check.txt'))
    fd.write('echo "Adapter contamination (total/per million records)"'
             '>> "$LOGFILE"\n')
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
    fd.write('mv %s %s\n' %
             (running_file(config, 'adaptor_check_running'),
              running_file(config, 'adaptor_check_done')))

    fd.close()


def had_adaptor_check_been_run(config):
    ''' Looks for the file that indicates that the adaptor check has been
        run. If true it reutnrs 'Yes'
    '''
    if os.path.isfile(running_file(config, 'adaptor_check_done')):
        return 'Yes'
    return 'No'


def do_run_fastqc_validator_on_original_fastq_files(config):
    ''' Runs the fastq validator module on the original fastq files.
    '''
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
    fd.write('#SBATCH --ntasks=16\n')
    fd.write('#SBATCH --time=0-05:00\n')
    fd.write('#SBATCH --mem-per-cpu=8000\n')
    fd.write('touch %s\n' %
             running_file(config, 'fastq_validator_check_running'))
    fd.write('# Modules\n')
    fd.write('module add FastQValidator\n')
    fd.write('# Defines\n')
    fd.write('LOGFILE="%s"\n' %
             os.path.join(original_qc_dir(config),
                          'fastq_validator_check.log'))
    fd.write('# Code to run the validation\n')
    fd.write('echo "Running fastq validator check" > "$LOGFILE"\n')
    for fname in find_files(original_files_dir(config), '*.fastq.gz'):
        fd.write('"fastQValidator" --file "%s" >> "$LOGFILE"\n' % fname)

    fd.write('mv %s %s\n' %
             (running_file(config, 'fastq_validator_check_running'),
              running_file(config, 'fastq_validator_check_done')))
    fd.close()


def had_fastq_been_validated(config):
    ''' Returns "Yes" if the fastq validator has been run by detecting the
        done file.
    '''
    if os.path.isfile(running_file(config, 'fastq_validator_check_done')):
        return 'Yes'
    return 'No'


def do_run_fastqc_on_original_fastq_files(config):
    ''' Runs the fastqc on the original fastq files.
    '''
    slrum_file = os.path.join(slrum_dir(config), 'fastq_original.slrum')

    # Make the QC directory for the results
    create_dir(original_qc_dir(config))
    qc_dir = os.path.join(original_qc_dir(config), 'fastq')
    create_dir(qc_dir)

    # Create the SLRUM file to run the adptor check
    print('Writing job file: %s' % slrum_file)
    fd = open(slrum_file, 'w')

    fd.write('#!/bin/bash --login\n')
    fd.write('#SBATCH --job-name=fastq_original_check\n')
    fd.write('#SBATCH --output=fastq_original_check_%j.out\n')
    fd.write('#SBATCH --error=fastq_original_check_%j.err\n')
    fd.write('#SBATCH --exclusive\n')
    fd.write('#SBATCH --ntasks=16\n')
    fd.write('#SBATCH --time=0-05:00\n')
    fd.write('#SBATCH --mem-per-cpu=8000\n')
    fd.write('touch %s\n' % running_file(config, 'fastq_original_running'))
    fd.write('# Modules\n')
    fd.write('module add FastQC/0.11.2\n')
    fd.write('# Code to run the validation\n')
    fd.write('echo "Running fastq on original files"\n')
    for fname in find_files(original_files_dir(config), '*.fastq.gz'):
        fd.write('fastqc "%s" -o "%s" --extract\n' % (fname, qc_dir))

    fd.write('mv %s %s\n' %
             (running_file(config, 'fastq_original_running'),
              running_file(config, 'fastq_original_done')))
    fd.close()


def had_fastq_original_been_run(config):
    ''' Returns "Yes" if the fastqc has been run on the original fastq files
        by looking for the done file.
    '''
    if os.path.isfile(running_file(config, 'fastq_original_done')):
        return 'Yes'
    return 'No'


def do_run_trim_and_pair_original_fastq_files(config):
    ''' Runs the trim and pair on the original fastq files using the
        project set quality score for trimming.
    '''
    quality = current_project_quality_score(config)
    slrum_file = os.path.join(slrum_dir(config),
                              'trim_and_pair_original_files.slrum')

    # Create the SLRUM file to run the adptor check
    print('Writing job file: %s' % slrum_file)
    fd = open(slrum_file, 'w')

    fd.write('#!/bin/bash --login\n')
    fd.write('#SBATCH --job-name=trim_and_pair_original_files\n')
    fd.write('#SBATCH --output=trim_and_pair_original_files_%j.out\n')
    fd.write('#SBATCH --error=trim_and_pair_original_files_%j.err\n')
    fd.write('#SBATCH --exclusive\n')
    fd.write('#SBATCH --ntasks=16\n')
    fd.write('#SBATCH --time=0-05:00\n')
    fd.write('#SBATCH --mem-per-cpu=8000\n')
    fd.write('touch %s\n' %
             running_file(config, 'trim_and_pair_original_files_running'))
    fd.write('# Modules\n')
    fd.write('module add Java\n')
    fd.write('module add Trimmomatic/0.33\n')
    fd.write('# Code to run the validation\n')
    fd.write('echo "Running trim and pair on original files"\n')

    unique = []
    for fname in find_files(original_files_dir(config), '*.fastq.gz'):
        head, tail = os.path.split(fname)
        unique.append('%s_R' % '_'.join(tail.split('_')[0:-2]))
    idens = list(set(unique))

    paired_dir = os.path.join(trim_paired_files_dir(config), 'paired')
    single_dir = os.path.join(trim_paired_files_dir(config), 'unpaired')

    create_dir(paired_dir)
    create_dir(single_dir)

    for iden in idens:
        fd.write('java -Xms64m -Xmx2000m -jar "%s" PE -threads 8 -phred33 '
                 '%s %s %s %s %s %s ILLUMINACLIP:"%s":2:30:10 '
                 'HEADCROP:3 LEADING:20 SLIDINGWINDOW:4:%s MINLEN:200\n' %
                 (trimmomatic_jar(),
                  os.path.join(original_files_dir(config),
                               '%s1_001.fastq.gz' % iden),
                  os.path.join(original_files_dir(config),
                               '%s2_001.fastq.gz' % iden),
                  os.path.join(paired_dir, '%s1_001P.fastq.gz' % iden),
                  os.path.join(single_dir, '%s1_001U.fastq.gz' % iden),
                  os.path.join(paired_dir, '%s2_001P.fastq.gz' % iden),
                  os.path.join(single_dir, '%s2_001U.fastq.gz' % iden),
                  os.path.join(adaptor_dir(config), 'NexteraPE-PE.fa'),
                  quality))
        fd.write('\n')

    fd.write('mv %s %s\n' % (
             running_file(config, 'trim_and_pair_original_files_running'),
             running_file(config, 'trim_and_pair_original_files_done')))
    fd.close()


def do_run_merge_on_trim_and_paired_files(config):
    ''' Runs the merge on the successful trim and paired fastq files
    '''
    slrum_file = os.path.join(slrum_dir(config),
                              'merge_trim_and_paired_files.slrum')

    # Create the SLRUM file to run the adptor check
    print('Writing job file: %s' % slrum_file)
    fd = open(slrum_file, 'w')

    fd.write('#!/bin/bash --login\n')
    fd.write('#SBATCH --job-name=merge_trim_and_paired_files\n')
    fd.write('#SBATCH --output=merge_trim_and_paired_files_%j.out\n')
    fd.write('#SBATCH --error=merge_trim_and_paired_files_%j.err\n')
    fd.write('#SBATCH --exclusive\n')
    fd.write('#SBATCH --ntasks=16\n')
    fd.write('#SBATCH --time=0-05:00\n')
    fd.write('#SBATCH --mem-per-cpu=8000\n')
    fd.write('touch %s\n' %
             running_file(config,
                          'merge_trim_and_paired_files_running'))
    fd.write('# Modules\n')
    fd.write('module add FLASH/1.2.11\n')
    fd.write('# Code to run the procedure\n')
    fd.write('echo "Running merge on the trim and paired files"\n')

    merged_dir = os.path.join(trim_paired_files_dir(config), 'paired-merged')
    merged_dir_qc = os.path.join(merged_dir, 'qc')
    merged_dir_merged = os.path.join(merged_dir, 'merged')
    merged_dir_not_merged = os.path.join(merged_dir, 'not-merged')
    paired_dir = os.path.join(trim_paired_files_dir(config), 'paired')

    unique = []
    for fname in find_files(paired_dir, '*.fastq.gz'):
        head, tail = os.path.split(fname)
        unique.append('%s_R' % '_'.join(tail.split('_')[0:-2]))
    idens = list(set(unique))

    create_dir(merged_dir)
    create_dir(merged_dir_merged)
    create_dir(merged_dir_qc)
    create_dir(merged_dir_not_merged)

    for iden in idens:
        fd.write('flash -d "%s" --max-overlap 450 --min-overlap 10 '
                 '--max-mismatch-density 0.25 -t 1 -z %s %s\n' %
                 (merged_dir,
                  os.path.join(paired_dir, '%s1_001P.fastq.gz' % iden),
                  os.path.join(paired_dir, '%s2_001P.fastq.gz' % iden)))
        fd.write('mv %s %s\n' %
                 (os.path.join(merged_dir, 'out.notCombined_1.fastq.gz'),
                  os.path.join(merged_dir_not_merged, '%s1_001P.fastq.gz'
                               % iden)))
        fd.write('mv %s %s\n' %
                 (os.path.join(merged_dir, 'out.notCombined_2.fastq.gz'),
                  os.path.join(merged_dir_not_merged,
                               '%s2_001P.fastq.gz' % iden)))
        fd.write('mv %s %s\n' %
                 (os.path.join(merged_dir, 'out.extendedFrags.fastq.gz'),
                  os.path.join(merged_dir_merged,
                               '%s_merged.fastq.gz' % iden)))
        fd.write('mv %s %s\n' %
                 (os.path.join(merged_dir, 'out.hist'),
                  os.path.join(merged_dir_qc, '%s_hist.txt' % iden)))
        fd.write('mv %s %s\n' %
                 (os.path.join(merged_dir, 'out.histogram'),
                  os.path.join(merged_dir_qc, '%s_histogram.txt' % iden)))
        fd.write('\n')

    fd.write('mv %s %s\n' %
             (running_file(config, 'merge_trim_and_paired_files_running'),
              running_file(config, 'merge_trim_and_paired_files_done')))
    fd.close()


def do_remove_merges_under_length(config):
    ''' Removes the merges that are under length. In this case the length is
        hard coded to 450 but in the future this could be changed to allow a
        user set length.
    '''
    slrum_file = os.path.join(slrum_dir(config),
                              'merged_files_length_selection.slrum')

    # Create the SLRUM file to run the adptor check
    print('Writing job file: %s' % slrum_file)
    fd = open(slrum_file, 'w')

    fd.write('#!/bin/bash --login\n')
    fd.write('#SBATCH --job-name=merged_files_length_selected\n')
    fd.write('#SBATCH --output=merged_files_length_selected_%j.out\n')
    fd.write('#SBATCH --error=merged_files_length_selected_%j.err\n')
    fd.write('#SBATCH --exclusive\n')
    fd.write('#SBATCH --ntasks=16\n')
    fd.write('#SBATCH --time=0-05:00\n')
    fd.write('#SBATCH --mem-per-cpu=8000\n')
    fd.write('touch %s\n' %
             running_file(config, 'merged_files_length_selected_running'))
    fd.write('# Modules\n')
    fd.write('module add Java\n')
    fd.write('module add Trimmomatic/0.33\n')
    fd.write('# Code to run the procedure\n')
    fd.write('echo "Running selection on merged files"\n')

    merged_dir = os.path.join(trim_paired_files_dir(config), 'paired-merged')
    merged_merged_dir = os.path.join(merged_dir, 'merged')
    selected_dir = os.path.join(merged_dir, 'merged-length-selected')
    files = find_files(merged_merged_dir, '*.fastq.gz')

    create_dir(selected_dir)

    for fname in files:
        head, tail = os.path.split(fname)
        fd.write('java -Xms64m -Xmx2000m -jar "%s" SE -threads 8 -phred33 %s'
                 '%s MINLEN:450\n' %
                 (trimmomatic_jar(), fname,
                  os.path.join(selected_dir, tail)))
        fd.write('\n')

    fd.write('mv %s %s\n' %
             (running_file(config, 'merged_files_length_selected_running'),
              running_file(config, 'merged_files_length_selected_done')))
    fd.close()


def do_convert_and_calapse_to_fasta(config):
    ''' Convert the length selected merged fastq files to fasta and
        collapse them down i.e. dereplicate them.
    '''
    slrum_file = os.path.join(
        slrum_dir(config),
        'selected_merged_files_convert_and_calapse.slrum')

    # Create the SLRUM file to run the adptor check
    print('Writing job file: %s' % slrum_file)
    fd = open(slrum_file, 'w')

    fd.write('#!/bin/bash --login\n')
    fd.write('#SBATCH --job-name=selected_merged_files_convert_and_calapse\n')
    fd.write(
        '#SBATCH --output=selected_merged_files_convert_and_calapse_%j.out\n')
    fd.write(
        '#SBATCH --error=selected_merged_files_convert_and_calapse_%j.err\n')
    fd.write('#SBATCH --exclusive\n')
    fd.write('#SBATCH --ntasks=16\n')
    fd.write('#SBATCH --time=0-05:00\n')
    fd.write('#SBATCH --mem-per-cpu=8000\n')
    fd.write('touch %s\n' %
             running_file(config,
                          'selected_merged_files_convert_and_calapse_running'))
    fd.write('# Modules\n')
    fd.write('module add fastx_toolkit/0.0.14\n')
    fd.write('# Code to run the procedure\n')
    fd.write('echo "Running convert and calapse to fasta on merged files"\n')

    paried_merged_dir = os.path.join(
        trim_paired_files_dir(config), 'paired-merged')
    merged_selected_dir = os.path.join(
        paried_merged_dir, 'merged-length-selected')
    fasta_dir = os.path.join(
        paried_merged_dir, 'merged-length-selected-fasta')
    fasta_calapsed_dir = os.path.join(
        paried_merged_dir, 'merged-length-selected-calapsed')
    files = find_files(merged_selected_dir, '*.fastq.gz')

    create_dir(fasta_dir)
    create_dir(fasta_calapsed_dir)

    for fname in files:
        head, tail = os.path.split(fname)
        fasta_name = '%s.fasta' % tail.split('.')[0]
        fd.write('zcat %s | fastq_to_fasta -Q33 -o %s\n' %
                 (fname, os.path.join(fasta_dir, fasta_name)))
        fd.write('fastx_collapser -i %s -o %s\n' %
                 (os.path.join(fasta_dir, fasta_name),
                  os.path.join(fasta_calapsed_dir, fasta_name)))
        fd.write('\n')

    fd.write('mv %s %s\n' %
             (running_file(
              config,
              'selected_merged_files_convert_and_calapse_running'),
              running_file(
              config,
              'selected_merged_files_convert_and_calapse_done')))
    fd.close()


def had_trim_and_pair_on_original_been_run(config):
    ''' Has the trim and pair on the original fastq files been run,
        return "Yes" if it has.
    '''
    if os.path.isfile(running_file(
                      config, 'trim_and_pair_original_files_done')):
        return 'Yes'
    return 'No'


def had_merge_on_trim_and_paired(config):
    ''' Has the merge on the trim and pair fastq been run yet, if it has
        return "Yes"
    '''
    if os.path.isfile(running_file(
                      config, 'merge_trim_and_paired_files_done')):
        return 'Yes'
    return 'No'


def had_merged_files_been_length_selected(config):
    ''' Has the merged files been selected for length yet? If it has then
        return "Yes"
    '''
    if os.path.isfile(running_file(
                      config, 'merged_files_length_selected_done')):
        return 'Yes'
    return 'No'


def had_selected_merged_files_been_converted_to_fasta(config):
    ''' Has the merged length selected files been converted to fasta yet?
        If they have then return 'Yes'
    '''
    if os.path.isfile(running_file(
                      config,
                      'selected_merged_files_convert_and_calapse_done')):
        return 'Yes'
    return 'No'


def top_menu(config):
    ''' Display the menu and run the corrosponding command. Each command
        will generate a slrum file that can be run with the sbatch command.
    '''
    key_dispatch = {'e': [do_end, 'Exit'],
                    'E': [do_end, 'Exit'],
                    '1': [do_set_current_project_dir,
                          'Set current project director'],
                    '2': [do_set_adaptor_dir,
                          'Set current adaptor director'],
                    '3': [do_set_quality_score,
                          'Set project quality score'],
                    '4': [do_setup_pipeline_in_current_project_dir,
                          'Set-up pipeline in current project director'],
                    '5': [do_run_adaptor_check_on_original_fastq_files,
                          'Generate adaptor check'],
                    '6': [do_run_fastqc_validator_on_original_fastq_files,
                          'Generate fastq validator'],
                    '7': [do_run_fastqc_on_original_fastq_files,
                          'Generate FastQC'],
                    '8': [do_run_trim_and_pair_original_fastq_files,
                          'Generate trim and pair'],
                    '9': [do_run_merge_on_trim_and_paired_files,
                          'Generate merge on trim and paired'],
                    'a': [do_remove_merges_under_length,
                          'Generate remove under length'],
                    'b': [do_convert_and_calapse_to_fasta,
                          'Generate fasta and collapse'], }

    menu_order = ['1', '2', '3', '4', '5', '6', '7', '8', '9',
                  'a', 'b', 'e']

    project_strings_and_state = [
        ['Current project dir',
         project_dir],
        ['Current adaptor dir',
         adaptor_dir],
        ['Project dir state',
         current_project_state],
        ['Project quality score',
         current_project_quality_score],
        ['Has adaptor check been run?',
         had_adaptor_check_been_run],
        ['Has fastq been validated?',
         had_fastq_been_validated],
        ['Has fastq been run on original?',
         had_fastq_original_been_run],
        ['Has trim and pair been run on original?',
         had_trim_and_pair_on_original_been_run],
        ['Has merge been run on trim and paired files?',
         had_merge_on_trim_and_paired],
        ['Has merged files been length selected?',
         had_merged_files_been_length_selected],
        ['Has selected merged files been collapsed to fasta?',
         had_selected_merged_files_been_converted_to_fasta]]

    while True:
        print('')
        print('Plant Illumina Pipeline')
        print('-----------------------')
        if 'project_dir' in config.section('general'):
            for string, func in project_strings_and_state:
                print('%s: %s', func(config))
            print('')

        for key in menu_order:
            description = key_dispatch[key][1]
            print(f'{key}: {description}')

        print('')
        print('Enter option:'),
        option = input()

        if option in key_dispatch:
            key_dispatch[option](config)


if __name__ == '__main__':
    config = Configuration()
    top_menu(config)
