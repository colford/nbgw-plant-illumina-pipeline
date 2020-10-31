#!/usr/bin/python
# -*- coding: utf-8 -*-
###############################################################################
# Reads in the CSV output of the Vsearch blast process. It combines replicates
# matches with the same IDs and splits the CSV file in to two files. One being
# the match data with BLAST output and the other being the list of IDs that
# go in to making up that row.
###############################################################################

import sys
import logging
import re
import os
import argparse
from collections import defaultdict
import xlsxwriter as xw


blast_spp_score = re.compile('.*\(\d+\)')

def create_dir(path):
    try:
        os.makedirs(path)
    except:
        pass

def extract_sample_id(centroid):
    parts = centroid.split('-')
    if '=' in parts[0]:
        return parts[0].split('=')[1].strip()
    return parts[0].strip()


class Vsearch_record(object):
    ''' A single vsearch record which is a line from the CSV file. There may
        not be an exact one to one match between this record and the line
        read in since we might find a match on ID and auto matched taxa, in
        which case we will compress it to a single record and add in the
        sequences and IDs that go to make up both records.
    '''

    def __init__(self, blast_row_csv):
        super(Vsearch_record, self).__init__()
        property_defaults = {
            'centroid': '',
            'sample_id': '',
            'number_of_sequences': 0,
            'number_combined': 0,
            'highest_bit_score': 0,
            'family_raw': '',
            'family': '',
            'genus': '',
            'species_raw': '',
            'species': '',
            'blast_raw': [],
            'combined_raw': []
        }
        self.__dict__.update(property_defaults)
        self._row = blast_row_csv
        self._parse_row()

    def __repr__(self):
        str = ("\n"
               "Centroid:                  %s\n"
               "Sample_id:                 %s\n"
               "Number sequences:          %s\n"
               "Number combined sequences: %s\n"
               "Auto match family:         %s\n"
               "Auto match genus:          %s\n"
               "Auto match species:        %s\n"
               "Auto match species raw:    %s\n"
               "Auto match family raw:     %s\n"
               "Blast output:              %s\n"
               "Combined with:             %s\n"
               %
               (self.centroid,
                self.sample_id,
                self.number_of_sequences,
                self.number_combined,
                self.family,
                self.genus,
                self.species,
                self.species_raw,
                self.family_raw,
                self.blast_raw,
                self.combined_raw))
        return str

    def __str__(self):
        return self.__repr__()

    def _extract_derep_size(self):
        ''' Parses out the number of depre sequences that make up the
            number of sequences.

            centroid=<id>;seqs=<derep>;size=<total seqs>
        '''
        self.number_combined = int(self.centroid.split(';')[1].split('=')[1])

    def _extract_sample_id(self):
        parts = self.centroid.split('=')[1].split(';')[0].split('-')
        self.sample_id = '-'.join(parts[0:2])

    def _extract_taxa(self):
        ''' Removes things like '%' or 'x' etc
            Then sets the family, genus, species as well as it can
        '''
        self._extract_family()
        self.species_raw = self.species_raw.strip('%').replace(' x', '')

        # The species will have 'Genus Species' or just genus
        parts = self.species_raw.split()
        self.genus = parts[0]
        if len(parts) > 2:
            raise Exception('ERROR: Gen/Spp has more than two!! %s' % parts)
        elif len(parts) == 2:
            self.species = parts[1]

        if self.genus == self.family:
            # We are here because we probably couldn't get to species
            # and hence have a load of 'family genus'
            self._extract_genus()

    def _extract_genus(self):
        ''' We only run this when the genus == family and hence we couldn't
            get to species.
        '''
        genus = []
        for val in set(self.family_raw):
            genus.append(val.split()[1])

        self.genus = sorted(list(set(genus)))
        self.genus = '/'.join(self.genus)

    def _extract_family(self):
        ''' There can be more than one family/genus recognised.
        '''
        family = []
        family_raw = []
        for taxa in self.family_raw:
            family_raw.append(taxa.strip('%').replace(' x', ''))
        self.family_raw = family_raw
        for val in set(self.family_raw):
            family.append(val.split()[0])

        self.family = sorted(list(set(family)))
        self.family = '/'.join(self.family)

    def __eq__(self, other):
        return self.centroid == other.centroid

    def merge(self, other):
        self.number_of_sequences += other.number_of_sequences
        self.number_combined += other.number_combined
        self.combined_raw.extend(other.combined_raw)

    def can_merge(self, other):
        return (self.sample_id == other.sample_id and
                self.family == other.family and
                self.genus == other.genus and
                self.species == other.species)

    def _parse_row(self):
        ''' We expect a row to consist of the following columns:

            0       Centroid code
            1       Number of sequences that make up this clusterd sequence
            2       Top bit score
            3       Genus/Species
            4       Family/Genus (list)
            blank   Start of Blast results
            5       Blast results (list)
            blank   Start of combined sequences that make up this centroid
            6       Combined sequences
            7       Number of sequences

            Note that the total number of sequences that make up this
            centroid should equal the total number of sequences summed for
            each ID.
        '''
        parts = self._row.split(',')
        self.centroid = parts[0]
        self.number_of_sequences = int(parts[1])
        self.highest_bit_score = float(parts[2])
        self.species_raw = parts[3]

        blast_start = -1
        combined_start = -1
        for inx, val in enumerate(parts):
            if blast_start < 0 and val == '':
                blast_start = inx
            elif combined_start < 0 and val == '':
                combined_start = inx

        self.family_raw = parts[4:blast_start]
        self.blast_raw = parts[blast_start + 1:combined_start]
        self.combined_raw = parts[combined_start + 1:]

        self._extract_sample_id()
        self._extract_derep_size()
        self._extract_taxa()


def debug_records(records):
    for sample_id in records:
        logger.debug('-> %s' % sample_id)
        for rec in records[sample_id]:
            try:
                logger.debug('%s' % rec)
            except Exception as x:
                logger.debug('ERROR: %s' % rec.centroid)


def merge(records):
    ''' Go through the records and merge any within a sample_id
        that have the same family/genus/species
    '''
    for sample_id in records:
        start_pos = 0
        logger.info('**** -> %s (%d)' % (sample_id, len(records[sample_id])))
        while True:
            merge_obj = records[sample_id][start_pos]
            logger.debug('Merging obj: %s' % merge_obj)
            used = []
            for rec_obj in records[sample_id][start_pos + 1:]:
                logger.debug('Testing obj: %s' % rec_obj)
                if merge_obj.can_merge(rec_obj):
                    logger.debug('Can merge!!!')
                    merge_obj.merge(rec_obj)
                    used.append(rec_obj)
                else:
                    logger.debug('**** No merge!!!')
            for obj in used:
                logger.debug('Deleting: %s' % obj.centroid)
                records[sample_id].remove(obj)
            logger.debug('<-> %d (%d)' % (start_pos, len(records[sample_id])))
            start_pos += 1
            if start_pos >= len(records[sample_id]):
                logger.debug('Final: %s' % records[sample_id])
                break


def output(outdir, primer, records):
    ''' Outputs the data for the manual process of assigning the correct
        family, genus, species and then after that we can make it in to a
        matrix for each sample vs species.
    '''
    fh = open('%s/%s-centroid-ids.csv' % (outdir, primer), 'w')
    wb = xw.Workbook('%s/%s-now-manual-edit.xlsx' % (outdir, primer))
    ws = wb.add_worksheet()

    static = ['Centroid',
              'Total Sequences',
              'Derep Sequences',
              'Family Decision',
              'Genus Decision',
              'Species Decision',
              'Auto Family',
              'Auto Genus',
              'Auto Species',
              'Blast References']

    row_count = 0
    ws.write_row(row_count, 0, tuple(static))
    row_count += 1

    for sid in records:
        for rec in records[sid]:
            row = [rec.centroid,
                   rec.number_of_sequences,
                   rec.number_combined,
                   '', '', '',
                   rec.family,
                   rec.genus,
                   rec.species]

            for blast in rec.blast_raw:
                row.append(blast)

            ws.write_row(row_count, 0, tuple(row))
            row_count += 1

            fh.write('%s,%s\n' % (rec.centroid, ','.join(rec.combined_raw)))

    fh.close()
    wb.close()


def process_csv(outdir, fname):
    ''' Create a dictionary of research records.
    '''
    records = defaultdict(list)
    primer = os.path.basename(fname).split('.')[0].split('-')[-1]
    with open(fname, 'r') as reader:
        for row in reader:
            rec = Vsearch_record(row.rstrip())
            records[rec.sample_id].append(rec)

    merge(records)
    debug_records(records)
    output(outdir, primer, records)


if __name__ == '__main__':
    log_level = logging.INFO

    argp = argparse.ArgumentParser(
        description='Process VSEARCH CSV to make it easier to use.')
    argp.add_argument(
        '-c',
        '--csv',
        dest='csv_fname',
        help='INPUT: Vsearch pipeline output CSV.',
        required=True)
    argp.add_argument(
            '--project_dir',
            dest='project_dir',
            required=True,
            help='Project directory')
    argp.add_argument(
        '-d',
        '--debug',
        action='store_true',
        help='Debug output')
    args = argp.parse_args()

    if args.debug:
        log_level = logging.DEBUG

    primer = os.path.basename(args.csv_fname).split('.')[0].split('-')[-1]
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
        datefmt='%m-%d %H:%M',
        filename='%s-process-vsearch-blast-output-for-excel.log' % primer, 
        filemode='w')
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    logging.getLogger('').addHandler(console)
    logger = logging.getLogger(__name__)

    logging.basicConfig(stream=sys.stdout, level=log_level)
    logger.debug('Debug is turned on...')

    man_dir = '%s/manual-checking' % args.project_dir
    logger.debug('Creating manual-checking directory: %s' % man_dir)
    create_dir(man_dir)
    process_csv(man_dir, args.csv_fname)
