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
import csv
from os import system
import argparse
import subprocess
import sys
import logging

# proviamo a vedere se il branch funziona
# vediamo come funziona il conflitto


# Define ks as a global variable
ks = [21, 31]

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
    
    # Define output file names based on sample name and k
    output_file = f"{output_path}{sample}_k{k}_clusters_stats.txt"
    skipped_file = f"{output_path}{sample}_k{k}_skipped-vcf-record.txt"
    
    clusters = defaultdict(dict)
    
    with open(skipped_file, 'w') as skipped:
        ill_vcf_reader_iterator = vcf.Reader(open(vcf_file, 'r'))
        ill_vcf_reader = []
        for record in ill_vcf_reader_iterator:
            alleles = []

            if str(record.genotype(sample)['GT']) == '1/1':
                continue
            
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
                
    clusters = {}
    c = 0
    for p in range(0, len(ill_vcf_reader)):
        if ill_vcf_reader[p]['CHROM'] not in clusters.keys():
            c = 0
            clusters[ill_vcf_reader[p]['CHROM']] = {}
            clusters[ill_vcf_reader[p]['CHROM']][c] = [ill_vcf_reader[p]]    
            continue
        if ill_vcf_reader[p]['AFF_START'] - ill_vcf_reader[p-1]['AFF_START'] <= k:
            clusters[ill_vcf_reader[p]['CHROM']][c].append(ill_vcf_reader[p])
            continue
        c += 1
        clusters[ill_vcf_reader[p]['CHROM']][c] = [ill_vcf_reader[p]]

    logging.info('Variants clustered successfully')

    # Collect cluster statistics and write to a separate stats file
    stats_file = f"{output_path}{sample}_k{k}_clusters_stats.txt"
    k_clusters_freq = defaultdict(int)
    tot_clusters = 0
        num_clusters_gt_1 = 0
    total_len_gt_1 = 0

    with open(stats_file, 'w') as stats:
        for chrom in clusters.keys():
            for c in clusters[chrom]:
                tot_clusters += 1
                l = len(clusters[chrom][c])
                k_clusters_freq[l] += 1
                if len(clusters[chrom][c]) > 1:
                    num_clusters_gt_1 += 1
                    total_len_gt_1 += l

        avg_len_gt_1 = total_len_gt_1 / num_clusters_gt_1 if num_clusters_gt_1 > 0 else 0
        keys_list = list(k_clusters_freq.keys())
        keys_list.sort()
        sorted_k_clusters_freq = {i: k_clusters_freq[i] for i in keys_list}
        stats.write(f'Tot clusters:\t{tot_clusters}\n')
        stats.write(f'Tot {k} clusters:\t{num_clusters_gt_1}\t({(num_clusters_gt_1/tot_clusters)*100:.2f})%\n')
        for freq in sorted_k_clusters_freq:
            stats.write(f'{freq}\t{sorted_k_clusters_freq[freq]}\n')

    return tot_clusters, num_clusters_gt_1, avg_len_gt_1

# Main function to iterate over VCF files and process them
def main(file_list_txt, output_path=''):
    vcf_files = read_vcf_list(file_list_txt)

    # Open the CSV file for writing
    csv_file = f"{output_path}cluster_statistics.csv"
    with open(csv_file, 'w', newline='') as csvfile:
        fieldnames = ['Sample', 'k', 'Total Clusters', 'Clusters > 1', 'Avg Len of Clusters > 1']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for vcf_file in vcf_files:
            sample_name = extract_sample_name(vcf_file)
            for k in ks:
                tot_clusters, num_clusters_gt_1, avg_len_gt_1 = cluster_variants(vcf_file, sample=sample_name, k=k, output_path=output_path)

                # Write the statistics to the CSV file
                writer.writerow({
                    'Sample': sample_name,
                    'k': k,
                    'Total Clusters': tot_clusters,
                    'Clusters > 1': num_clusters_gt_1,
                    'Avg Len of Clusters > 1': avg_len_gt_1
                })

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process VCF files.')
    parser.add_argument('file_list', help='Text file containing list of VCF files')
    parser.add_argument('--output_path', default='', help='Output directory path')
    args = parser.parse_args()

    main(args.file_list, args.output_path)

