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

        for record in ill_vcf_reader_iterator:
            alleles = []
            if str(record.genotype(sample)['GT']) in ['1/0', '0/1']:
                alleles.append((record.POS-1, int(record.POS-1 + len(record.REF)), record.REF))
                
            for alt in record.ALT:
                check_alt = special_match(str(alt))
                if not check_alt:
                    alleles = []
                    skipped.write(str(record) + '\n')
                    break
                alleles.append((record.POS - 1, int(record.POS - 1 + len(record.REF)), str(alt)))

            if alleles:
                chrom = record.CHROM
                if chrom not in clusters:
                    clusters[chrom] = defaultdict(list)
                clusters[chrom][record.POS].extend(alleles)
                ill_vcf_reader.append({
                    'CHROM': record.CHROM,
                    'AFF_START': record.POS - 1,
                    'END': int(record.POS-1 + len(record.REF)),
                    'REF_LEN': len(record.REF),
                    'GT': record.genotype(sample)['GT'],
                })

        # Write statistics or results to the output file
        output.write(f"Sample: {sample}, k: {k}\n")
        output.write(f"Total variants processed: {len(ill_vcf_reader)}\n")
        
        for variant in ill_vcf_reader:
            output.write(f"CHROM: {variant['CHROM']}, Start: {variant['AFF_START']}, End: {variant['END']}\n")

    # Collect cluster statistics and write to a separate stats file
    stats_file = f"{output_path}{sample}_clusters_statistics.txt"
    with open(stats_file, 'w') as stats:
        for k in ks:
            clusters = cluster_variants(vcf_file, sample, k, output_path)
            k_clusters_freq = defaultdict(int)
            stats.write(f'Stats for sample {sample} with {k}-mers clusters: \n')
            tot_clusters = 0
            tot_k_clusters = 0
            for chrom in clusters.keys():
                for c in clusters[chrom]:
                    tot_clusters += 1
                    l = len(clusters[chrom][c])
                    k_clusters_freq[l] += 1
                    if len(clusters[chrom][c]) > 1:
                        tot_k_clusters += 1
            keys_list = list(k_clusters_freq.keys())
            keys_list.sort()
            sorted_k_clusters_freq = {i: k_clusters_freq[i] for i in keys_list}
            stats.write(f'Tot clusters:\t{tot_clusters}\n')
            stats.write(f'Tot {k} clusters:\t{tot_k_clusters}\t({(tot_k_clusters/tot_clusters)*100:.2f})%\n')
            for freq in sorted_k_clusters_freq:
                stats.write(f'{freq}\t{sorted_k_clusters_freq[freq]}\n')
    
    # Clean up and garbage collection if necessary
    gc.collect()
    return clusters

# Main function to iterate over VCF files and process them
def main(file_list_txt, output_path=''):
    vcf_files = read_vcf_list(file_list_txt)

    for vcf_file in vcf_files:
        sample_name = extract_sample_name(vcf_file)
        for k in ks:
            cluster_variants(vcf_file, sample=sample_name, k=k, output_path=output_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process VCF files.')
    parser.add_argument('file_list', help='Text file containing list of VCF files')
    parser.add_argument('--output_path', default='', help='Output directory path')
    args = parser.parse_args()
    
    main(args.file_list, args.output_path)

