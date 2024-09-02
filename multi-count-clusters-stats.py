import vcf 
from Bio import SeqIO
import itertools
from itertools import chain, combinations, product
import numpy as np
import gc
import math
from collections import defaultdict

import re
import gzip
import pickle

from os import system

import argparse
import subprocess
import sys
import logging

# Define ks as a global variable
ks = [21,31]

# Function to extract sample name from a VCF file name
def extract_sample_name(vcf_file_name):
    match = re.search(r'HG00\d+', vcf_file_name)
    if match:
        return match.group(0)
    else:
        raise ValueError(f"Sample name not found in {vcf_file_name}")

# Function to read VCF file names from a text file
def read_vcf_list(file_list_txt):
    with open(file_list_txt, 'r') as file:
        vcf_files = file.read().splitlines()
    return vcf_files

# Define special_match function
def special_match(strg, search=re.compile(r'[a-zA-Z]').search):
    return bool(search(strg))

# Function to cluster variants for each VCF file
def cluster_variants(vcf_file, sample='Sample', k=21, output_path=''):
    logging.info(f'Clustering variants from VCF file {vcf_file}')
    
    # Define output file name based on sample name and k
    output_file = f"{output_path}{sample}_k{k}_clusters_stats.txt"
    skipped_file = f"{output_path}{sample}_k{k}_skipped-vcf-record.txt"
    
    clusters = defaultdict(dict)
    
    with open(skipped_file, 'w') as skipped, open(output_file, 'w') as output:
        ill_vcf_reader_iterator = vcf.Reader(open(vcf_file, 'r'))
        ill_vcf_reader = []
multi-count-clusters-stats.py

