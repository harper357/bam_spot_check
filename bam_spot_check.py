#!/usr/bin/env python3
"""
A tool to spot check BAM files by sending 10 reads to the NCBI BLAST server.

Input:
    BAM file OR FASTA file OR nucleotide sequence

Requirements:
    Python 3.5 or above
    samtools

Written by: Greg Fedewa
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
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 add_help=True
                                 )
input_group = parser.add_mutually_exclusive_group(required=True)
input_group.add_argument('-b', '--bam', type=str, required=False,
                    help='BAM file to spot check, including path.')
input_group.add_argument('-f', '--fasta', type=str, required=False,
                    help='FASTA file to spot check, including path.')
input_group.add_argument('-s', '--seq', type=str, required=False,
                    help='Sequence to spot check.')
#parser.add_argument('-f', '--fasta', type=str, required=True,
#                    choices=['F', 'A', 'D', 'P'], nargs='*',
#                    help='Programs you want the pipeline to run')

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
    fraction_size = str(wanted_read_count / int(subprocess.run(['samtools', 'view', '-c', bam_file],
                                                               stdout=subprocess.PIPE).stdout))
    wanted_bam = subprocess.run(['samtools', 'view', '-s', fraction_size, bam_file],
                                universal_newlines=True,
                                stdout=subprocess.PIPE).stdout
    batch_string = '>\n'+'\n>\n'.join([each_line.split('\t')[9] for each_line in wanted_bam.rstrip().split('\n')])
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
    RTOE_value = int(list(filter(lambda x:'RTOE = ' in x, query.text.splitlines()))[0].split(' ')[-1])
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
    sleep(max([delay_timer,3]))
    while True:            
            search_info = requests.get('https://blast.ncbi.nlm.nih.gov/Blast.cgi',
                                         params=search_info_parameters)
            status_value = list(filter(lambda x:'Status=' in x, search_info.text.splitlines()))[0].split('=')[-1]
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
    # Incase you are trying to make sence of this below, here is how it breaks down:
    # for each_query in root[8]:
    #   for x in each_query.findall(".//Hit_def")[0:5]:
    #     ' '.join(x.text.split()[0:2])
    # List comprehensions are weird and hard...
    hits = [' '.join(x.text.split()[0:2]) for each_query in root[8] for x in each_query.findall(".//Hit_def")[0:1]]
    return hits


def main():
    try:
        email_file = open('{}/bam_spot_check.email'.format(os.path.split(os.path.realpath(__file__))[0]),'r')
        current_email = email_file.read()
    except FileNotFoundError:
        email_file = open('{}/bam_spot_check.email'.format(os.path.split(os.path.realpath(__file__))[0]),'w')
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

    kwargs = {'MEGABLAST':'on'}
    if args.bam:
        hit_list = get_hits(*send_query(bam_to_fasta(args.bam, 10), EMAIL=current_email, **kwargs), EMAIL=current_email)
    elif args.fasta:
        input_query = open(args.fasta, 'r').read().rstrip()
        hit_list = get_hits(*send_query(input_query, EMAIL=current_email, **kwargs), EMAIL=current_email)
    elif args.seq:
        hit_list = get_hits(*send_query(args.seq, EMAIL=current_email, **kwargs), EMAIL=current_email)
    for each_hit in hit_list:
        print(each_hit)

if __name__ == '__main__':
    main()
