#!/usr/bin/env python3
"""
A tool to spot check BAM files by sending 10 reads to the NCBI BLAST server.

Input:
    BAM file OR FASTA file OR nucleotide sequence OR RID from previous BLAST job

Requirements:
    Python 3.5 or above
    samtools

Written by: Greg Fedewa
Version: 0.2.0
Copyright: 2017.11.15
"""

import requests
import xml.etree.ElementTree as ET
import subprocess
import argparse
import logging
import sys
from time import sleep as sleep
import os


parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter,
                                 add_help=True
                                 )
input_group = parser.add_mutually_exclusive_group(required=True)
input_group.add_argument('-b', '--bam', type=str, required=False,
                         help='BAM file to spot check, including path.')
input_group.add_argument('-f', '--fastx', type=str, required=False,
                         help='FASTA file to spot check, including path.')
input_group.add_argument('-s', '--seq', type=str, required=False,
                         help='Sequence to spot check.')
input_group.add_argument('-R', '--RID', type=str, required=False,
                         help='RID of BLAST job to look up.')
parser.add_argument('-d', '--details', type=str, required=False,
                    choices=['F', 'A', 'S', 'C'], nargs='+',
                    help='Additional details to print: \n\
    F = Full name of hit \n\
    A = Accession number \n\
    S = Score values \n\
    C = Complete details *not finished*')
### below comes from http://stackoverflow.com/q/14097061/78845
parser.add_argument("-v", "--verbose", help="increase output verbosity",
                    action="store_true")
args = parser.parse_args()
if args.verbose:
    logging.basicConfig(level=logging.DEBUG)

logging.debug('Only shown in debug mode')

###


def bam_to_fasta(bam_file, wanted_read_count):
    """Generates the string that will be sent to NCBI from a BAM file. """
    fraction_size = '{0:.15f}'.format(wanted_read_count / int(subprocess.run(['samtools', 'view', '-c', bam_file],
                                                               stdout=subprocess.PIPE).stdout))
    wanted_bam = subprocess.run(['samtools', 'view', '-s', fraction_size, bam_file],
                                universal_newlines=True,
                                stdout=subprocess.PIPE).stdout
    batch_string = '>\n'+'\n>\n'.join([each_line.split('\t')[9] for each_line in wanted_bam.rstrip().split('\n')])
    return batch_string


def fastx_subsample(fastx_file, wanted_read_count):
    """Parse the FASTx file to get some sequences to send.
    It currently decides if it is a FASTA vs FASTQ entry by entry.
    Please don't mix FASTQ and FASTA entries...
    """
    fastx_data = open(fastx_file, 'r')
    sequence_list = []
    counter = wanted_read_count
    temp_sequence = ''
    # the following doesn't seem very pythonic... but it should work...
    for _ in fastx_data:  # start going through the FASTx file
        if _.startswith('>'):  # Found a FASTA entry!
            while not next(fastx_data).startswith(('>', '@', '+')):  # for the love of everything, please don't use mixed FASTQ/FASTA files...
                temp_line = fastx_data.readline().rstrip()
                temp_sequence = temp_sequence + temp_line
            sequence_list.append(temp_sequence)
            temp_sequence = ''
            counter -= 1
        elif _.startswith('@'):  # Found a FASTQ entry!
            temp_sequence = fastx_data.readline().rstrip()  # After the name comes the sequence
            fastx_data.readline()  # After sequences comes '+', let us skip it
            fastx_data.readline()  # After '+' comes quality, let us skip it
            sequence_list.append(temp_sequence)
            temp_sequence = ''
            counter -= 1
        if counter <= 0:
            break
    batch_string = '>\n'+'\n>\n'.join(sequence_list)
    return batch_string


def send_query(QUERY, EMAIL, PROGRAM='blastn', DATABASE='nr', **kwards):
    """Sends a search query and returns the RID used to get the search results. """
    blast_api_parameters = {'CMD': 'Put',
                            'QUERY': QUERY,
                            'PROGRAM': PROGRAM,
                            'DATABASE': DATABASE,
                            'FILTER': 'L',
                            'FORMAT_TYPE': 'XML',
                            'EMAIL': EMAIL,
                            'TOOL': 'bam_spot_check.py'}
    if 'MEGABLAST' in kwards:
        blast_api_parameters['MEGABLAST'] = 'on'
    query = requests.get('https://blast.ncbi.nlm.nih.gov/Blast.cgi',
                         params=blast_api_parameters)
    RID_value = list(filter(lambda x: 'RID = ' in x, query.text.splitlines()))[0].split(' ')[-1]
    RTOE_value = int(list(filter(lambda x: 'RTOE = ' in x, query.text.splitlines()))[0].split(' ')[-1])
    print('Job {} has been given an estimated completion time of {} sec'.format(RID_value, RTOE_value))
    return RID_value, RTOE_value


def get_hits(RID, RTOE, EMAIL):
    """Generates a list of hits from the RID given."""
    blast_results_parameters = {'CMD': 'Get',
                                'FORMAT_TYPE': 'XML',
                                'NCBI_GI': 'F',
                                'RID': RID,
                                'HITLIST_SIZE': '1',
                                'EMAIL': EMAIL,
                                'TOOL': 'bam_spot_check.py'}
    search_info_parameters = {'CMD': 'Get',
                                'FORMAT_OBJECT': 'SearchInfo',
                                'RID': RID,
                                'EMAIL': EMAIL,
                                'TOOL': 'bam_spot_check.py'}
    delay_timer = RTOE
    sleep(max([delay_timer, 3]))
    while True:
            search_info = requests.get('https://blast.ncbi.nlm.nih.gov/Blast.cgi',
                                         params=search_info_parameters)
            status_value = list(filter(lambda x: 'Status=' in x, search_info.text.splitlines()))[0].split('=')[-1]
            if status_value == 'READY':
                break
            elif status_value == 'WAITING':
                sleep(delay_timer)
                print('Now waiting {} seconds'.format(delay_timer))
                delay_timer = min([delay_timer*1.5, 120])
            elif status_value == 'UNKNOWN':
                sys.exit('Blast status: UKNOWN')
    sleep(3)
    blast_results = requests.get('https://blast.ncbi.nlm.nih.gov/Blast.cgi',
                                         params=blast_results_parameters)
    root = ET.fromstring(blast_results.text)
    hits = []
    for BLAST_hits in root.findall('BlastOutput_iterations/Iteration/Iteration_hits'):
        hit = BLAST_hits.find('Hit')
        if hit is not None:
            if not args.details:
                short_name = [' '.join(hit.find('Hit_def').text.split()[0:2])]
                hits.append(short_name)
            else:
                temp_hit = []
                if 'F' in args.details:
                    # full name of hit
                    temp_hit.append(hit.find('Hit_def').text)
                else:
                    temp_hit.append(' '.join(hit.find('Hit_def').text.split()[0:2]))
                if 'A' in args.details:
                    # Accession of hit
                    temp_hit.append(hit.find('Hit_accession').text)
                if 'S' in args.details:
                    # bit score
                    temp_hit.append(next(hit.iter('Hsp_bit-score')).text)
                    # e-value
                    temp_hit.append(next(hit.iter('Hsp_evalue')).text)
                    # query coverage
                    temp_hit.append(str(round(int(next(hit.iter('Hsp_align-len')).text) /
                                        float(next(hit.iter('Hsp_query-to')).text), 2)*100))
                    # identities
                    temp_hit.append(str(round(int(next(hit.iter('Hsp_positive')).text) /
                                        float(next(hit.iter('Hsp_align-len')).text), 2)*100))
                if 'C' in args.details:
                    print('C is not yet supported')
                hits.append(temp_hit)
        else:
            NA_length = 0
            if not args.details:
                temp_hit = ['No significant similarity found']
                hits.append(temp_hit)
            else:
                if 'A' in args.details:
                    NA_length += 1
                if 'S' in args.details:
                    NA_length += 4
                if 'C' in args.details:
                    NA_length = 5
                temp_hit = ['No significant similarity found']
                temp_hit.extend(['N/A']*NA_length)
                hits.append(temp_hit)
    return hits


def main():
    #  Make sure that the user has a valid email to send to the BLAST server
    try:
        email_file = open('{}/bam_spot_check.email'.format(os.path.split(os.path.realpath(__file__))[0]), 'r')
        current_email = email_file.read()
    except FileNotFoundError:
        email_file = open('{}/bam_spot_check.email'.format(os.path.split(os.path.realpath(__file__))[0]), 'w')
        current_email = ''
    if current_email == '':
        while True:
            print('Please enter a valid email:')
            email1 = input()
            print('Confirm by reentering email:')
            email2 = input()
            if email1 == email2:
                current_email = email1
                email_file.write(email1)
                email_file.close()
                print('Email saved')
                break
            elif email1 != email2:
                print('Emails don\'t match.')

    kwargs = {'MEGABLAST': 'on'}
    if args.bam:
        hit_list = get_hits(*send_query(bam_to_fasta(args.bam, 10), EMAIL=current_email, **kwargs), EMAIL=current_email)
    elif args.fastx:
        # need to check if FASTA or FASTQ
        # need to only get the first ten?
        input_query = fastx_subsample(args.fastx, 10)
        hit_list = get_hits(*send_query(input_query, EMAIL=current_email, **kwargs), EMAIL=current_email)
    elif args.seq:
        hit_list = get_hits(*send_query(args.seq, EMAIL=current_email, **kwargs), EMAIL=current_email)
    elif args.RID:
        hit_list = get_hits(RID=args.RID, RTOE=0, EMAIL=current_email)
    # make print pretty if you have any details to print
    column_names = []
    #if args.details:
        # full name of hit
        #column_names.append('Hit name')
    column_names.append('Hit name')
    if args.details:
        if 'A' in args.details:
            column_names.append('Accession')
        if 'S' in args.details:
            column_names.append('bit score')
            column_names.append('e-value')
            column_names.append('Coverage')
            column_names.append('Identity')
        if 'C' in args.details:
            pass
    temp_cols = [column_names] + hit_list
    #modified from https://stackoverflow.com/questions/9989334/create-nice-column-output-in-python
    column_widths = [max(map(len, column)) for column in zip(*temp_cols)]
    # Print the titles in bold and underlined
    print('\033[1m'+'\033[4m' +
          '  '.join((title.ljust(spacing) for title, spacing in zip(column_names, column_widths))) +
          '\033[0m')
    for each_hit in hit_list:
        print('  '.join((detail.ljust(spacing) for detail, spacing in zip(each_hit, column_widths))))


if __name__ == '__main__':
    main()
