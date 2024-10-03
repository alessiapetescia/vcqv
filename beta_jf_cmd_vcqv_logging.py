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

sys.path.append('/home/alessiap/lib/python3.10/')
import dna_jellyfish

# Setup logging
def setup_logging(log_file):
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )

def load_ref(ref, ext='fasta'):
    logging.info('Started loading the reference')
    reference = SeqIO.parse(open(ref), ext)
    ref = {}
    genome_len = 0
    
    for fasta in reference:
        name, sequence = fasta.id, str(fasta.seq)
        ref[name] = sequence
        genome_len += len(sequence)
    
    logging.info('Reference loaded successfully')
    return ref, genome_len

def run_jf_count(reads_file, k, output, threads, min_count):
    logging.info(f'Running Jellyfish count for {reads_file}')
    cmd = [
        'jellyfish', 'count',
        '-m', str(k),
        '-s', '50M',
        '-t', str(threads),
        '-o', output,
        '-c', str(min_count),
        '-C',
        reads_file
    ]
    subprocess.run(cmd, check=True)
    logging.info(f'Jellyfish count completed and output saved to {output}')

def load_query_table(file):
    logging.info(f'Loading query table from {file}')
    qf = dna_jellyfish.QueryMerFile(file)
    logging.info(f'Query table loaded from {file}')
    return qf

def load_iterable_jf(file):
    logging.info(f'Loading iterable from {file}')
    it = dna_jellyfish.ReadMerFile(file)
    logging.info(f'Iterable loaded from {file}')
    return it

def query_hash_table(h, kmer):
    mer = dna_jellyfish.MerDNA(kmer)
    mer.canonicalize()
    return h[mer]

def cartesian_product_generator(iterables):
    #logging.info('Generating Cartesian product of iterables')
    if len(iterables) == 0:
        yield ()
        return
    
    stack = [(0, ())]
    while stack:
        idx, current_combination = stack.pop()
        if idx == len(iterables):
            yield current_combination
            continue
        
        for item in iterables[idx]:
            new_combination = current_combination + (item,)
            stack.append((idx + 1, new_combination))

def get_phased_variants(ref, alleles, read_mers, alleles_genotype, cartesian_product_generator=cartesian_product_generator, k=21):
    #logging.info('Phasing variants')
    best_haps = []
    
    for i in range(0, len(alleles_genotype)):
        if alleles_genotype[i] == '1/2' or alleles_genotype[i] == '2/1':
            alleles_genotype[i] = '0/1'
        elif alleles_genotype[i] == '1/0':
            alleles_genotype[i] = '0/1'
    
    best_score = float('inf')
    combs = cartesian_product_generator(alleles)

    for comb in combs:
        start_hap_pos = comb[0][0] - k + 1
        hap = [] 
        prev_pos = start_hap_pos  
        for var in comb:
            hap.append(ref[prev_pos:var[0]])
            hap.append(var[2])
            prev_pos = var[1]
        
        end_hap_pos = prev_pos + k - 1
        hap.append(ref[prev_pos:end_hap_pos])
        haplotype = ''.join(hap)
        errors = 0
        
        for p in range(0, len(haplotype) - k + 1):
            kmer_seq = haplotype[p:p+k]
            mer_count = query_hash_table(read_mers, kmer_seq)
            if not mer_count:
                errors += 1
        
        haplotype = haplotype[k-1: len(haplotype)-k+1]

        if errors < best_score:
            best_haps = [(start_hap_pos, end_hap_pos, haplotype, errors, comb)]
            best_score = errors
        elif errors == best_score:
            curr_pair_genotype = []
            for v1, v2 in zip(best_haps[0][4], comb):
                var1, var2 = hash(v1), hash(v2)
                if var1 == var2:
                    curr_pair_genotype.append('1/1')
                else:
                    curr_pair_genotype.append('0/1')
            if curr_pair_genotype == alleles_genotype:
                best_haps = [best_haps[0], (start_hap_pos, end_hap_pos, haplotype, errors, comb)]
    
    if len(best_haps) == 2:
        best_haps.append('1/2')
    else:
        best_haps.append('1/1')
        
    #logging.info('Phased variants identified')
    return best_haps

def phase_chunks(ref, alleles, read_mers, alleles_genotype=None, cartesian_product_generator=cartesian_product_generator, k=21, is_first_chunk=None, is_last_chunk=None, get_genotype=None):
    #logging.info('Phasing chunks')
    best_haps = []
    best_score = float('inf')
    combs = cartesian_product_generator(alleles)
    
    if not get_genotype:
        for i in range(0, len(alleles_genotype)):
            if alleles_genotype[i] == '1/2' or alleles_genotype[i] == '2/1':
                alleles_genotype[i] = '0/1'
            elif alleles_genotype[i] == '1/0':
                alleles_genotype[i] = '0/1'
                
    for comb in combs:
        if is_first_chunk:
            start_hap_pos = comb[0][0] - k + 1
        else:
            start_hap_pos = comb[0][0]
        
        hap = [] 
        prev_pos = start_hap_pos
        for var in comb:
            hap.append(ref[prev_pos:var[0]])
            hap.append(var[2])
            prev_pos = var[1]
        
        if is_last_chunk:
            end_hap_pos = prev_pos + k - 1
            hap.append(ref[prev_pos:end_hap_pos])
        else:
            end_hap_pos = prev_pos
        haplotype = ''.join(hap)
        errors = 0
        
        for p in range(0, len(haplotype) - k + 1):
            kmer_seq = haplotype[p:p+k]
            mer_count = query_hash_table(read_mers, kmer_seq)
            if not mer_count:
                errors += 1
        
        if get_genotype:
            if errors < best_score:
                best_haps = [(start_hap_pos, end_hap_pos, haplotype, errors, comb)]
                best_score = errors
            elif errors == best_score:
                best_haps = [best_haps[0], (start_hap_pos, end_hap_pos, haplotype, errors, comb)]
        else:
            if errors < best_score:
                best_haps = [(start_hap_pos, end_hap_pos, haplotype, errors, comb)]
                best_score = errors
            elif errors == best_score:
                curr_pair_genotype = []
                for v1, v2 in zip(best_haps[0][4], comb):
                    var1, var2 = hash(v1), hash(v2)
                    if var1 == var2:
                        curr_pair_genotype.append('1/1')
                    else:
                        curr_pair_genotype.append('0/1')
                if curr_pair_genotype == alleles_genotype:
                    best_haps = [best_haps[0], (start_hap_pos, end_hap_pos, haplotype, errors, comb)]
    
    if get_genotype:
        if len(best_haps) == 2:
            best_haps = [list(best_haps[0]), list(best_haps[1])]
            best_haps[0][0], best_haps[0][1], best_haps[0][2] = best_haps[0][0]+k-1, best_haps[0][1]-k+1, best_haps[0][2][k-1:len(best_haps[0][2])-k+1]
            best_haps[1][0], best_haps[1][1], best_haps[1][2] = best_haps[1][0]+k-1, best_haps[1][1]-k+1, best_haps[1][2][k-1:len(best_haps[1][2])-k+1]
            best_haps.append('1/2')
        else:
            best_haps = [list(best_haps[0])]
            best_haps[0][0], best_haps[0][1], best_haps[0][2] = best_haps[0][0]+k-1, best_haps[0][1]-k+1, best_haps[0][2][k-1:len(best_haps[0][2])-k+1]
            best_haps.append('1/1')

    #logging.info('Chunks phased successfully')
    return best_haps

def get_phased_haplotypes(ref, alleles, read_mers, alleles_genotype, get_phased_variants=get_phased_variants, cartesian_product_generator=cartesian_product_generator, k=21, max_vars=15):
    #logging.info('Getting phased haplotypes')
    if len(alleles) <= max_vars:
        best_haps = get_phased_variants(ref, alleles, read_mers, alleles_genotype, cartesian_product_generator, k)
    else:
        chunks = [alleles[x:x+max_vars] for x in range(0, len(alleles), max_vars)]
        chunks_genotype = [alleles_genotype[x:x+max_vars] for x in range(0, len(alleles), max_vars)]
        best_partial_haps = []
        
        best_first_chunk = phase_chunks(ref, chunks[0], read_mers, chunks_genotype[0], cartesian_product_generator, k, is_first_chunk='y')
        best_partial_haps.append(best_first_chunk)
        
        for c in range(1, len(chunks) - 1):
            chunk = chunks[c]
            best_inner_chunk = phase_chunks(ref, chunk, read_mers, chunks_genotype[c], cartesian_product_generator, k)
            best_partial_haps.append(best_inner_chunk)
        
        best_last_chunk = phase_chunks(ref, chunks[-1], read_mers, chunks_genotype[-1], cartesian_product_generator, k, is_last_chunk='y')
        best_partial_haps.append(best_last_chunk)
        
        best_haps = phase_chunks(ref, best_partial_haps, read_mers, chunks_genotype, cartesian_product_generator, k, get_genotype='y')
 
    #logging.info('Phased haplotypes obtained')
    return best_haps

def special_match(strg, search=re.compile(r'[a-zA-Z]').search):
    return bool(search(strg))

def cluster_variants(vcf_file, sample='Sample', k=21, output_path=''):
    logging.info(f'Clustering variants from VCF file {vcf_file}')
    with open(output_path + 'skipped-vcf-record.txt','w') as skipped:
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
                ill_vcf_reader.append({
                    'CHROM': record.CHROM,
                    'AFF_START': record.POS - 1,
                    'END': int(record.POS-1 + len(record.REF)),
                    'REF_LEN': len(record.REF),
                    'GT': record.genotype(sample)['GT'],
                    'ALLELES': alleles
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
    return clusters

def create_pseudo_haplotypes(ref, clusters, read_mers, output_haps='pseudo_haps.fasta'):
    logging.info('Creating pseudo-haplotypes')
    with open(output_haps, 'w') as haps_file:
        haps_file.write('')
        best_haps_list = []

    with open(output_haps, 'a') as haps_file:
        for chromosome in clusters.keys():
            logging.info(f'Processing chromosome {chromosome}')
            pseudo_hap_1 = list(ref[chromosome])
            pseudo_hap_2 = list(ref[chromosome])

            for cluster in list(clusters[chromosome].keys()):
                if len(clusters[chromosome][cluster]) == 1:
                    curr_pos = clusters[chromosome][cluster][0]['AFF_START']
                    len_ref_allele = clusters[chromosome][cluster][0]['REF_LEN']

                    if clusters[chromosome][cluster][0]['GT'] == '0/1':
                        pseudo_hap_2[curr_pos:curr_pos+len_ref_allele] = ['']*len_ref_allele
                        alt_allele = clusters[chromosome][cluster][0]['ALLELES'][-1][2]
                        pseudo_hap_2[curr_pos] = alt_allele

                    elif clusters[chromosome][cluster][0]['GT'] == '1/1':
                        pseudo_hap_1[curr_pos:curr_pos+len_ref_allele] = pseudo_hap_2[curr_pos:curr_pos+len_ref_allele] = ['']*len_ref_allele
                        alt_allele = clusters[chromosome][cluster][0]['ALLELES'][-1][2]
                        pseudo_hap_1[curr_pos] = pseudo_hap_2[curr_pos] = alt_allele

                    elif clusters[chromosome][cluster][0]['GT'] == '1/2':
                        pseudo_hap_1[curr_pos:curr_pos+len_ref_allele] = pseudo_hap_2[curr_pos:curr_pos+len_ref_allele] = ['']*len_ref_allele
                        alt_allele_1 = clusters[chromosome][cluster][0]['ALLELES'][-2][2]
                        alt_allele_2 = clusters[chromosome][cluster][0]['ALLELES'][-1][2]
                        pseudo_hap_1[curr_pos] = alt_allele_1
                        pseudo_hap_2[curr_pos] = alt_allele_2

                elif len(clusters[chromosome][cluster]) > 1:
                    alleles_to_combine = [dic['ALLELES'] for dic in clusters[chromosome][cluster]]
                    alleles_genotype = [dic['GT'] for dic in clusters[chromosome][cluster]]
                    
                    best_haps = get_phased_haplotypes(ref[chromosome], alleles_to_combine, read_mers, alleles_genotype)
                    best_haps_list.append(best_haps)
                    curr_pos = best_haps[0][0]
                    len_ref_allele = best_haps[0][1] - curr_pos
                    pseudo_hap_1[curr_pos:curr_pos + len_ref_allele] = pseudo_hap_2[curr_pos:curr_pos + len_ref_allele] = ['']*len_ref_allele

                    if best_haps[-1] == '1/1':
                        alt_allele = best_haps[0][2]
                        pseudo_hap_1[curr_pos] = pseudo_hap_2[curr_pos] = alt_allele

                    elif best_haps[-1] == '1/2':
                        alt_allele_1 = best_haps[0][2]
                        alt_allele_2 = best_haps[1][2]
                        pseudo_hap_1[curr_pos] = alt_allele_1
                        pseudo_hap_2[curr_pos] = alt_allele_2

            pseudo_hap_1 = ''.join(pseudo_hap_1)
            pseudo_hap_2 = ''.join(pseudo_hap_2)

            haps_file.write(f'>{chromosome}_pseudo_hap_1\n{pseudo_hap_1}\n')
            haps_file.write(f'>{chromosome}_pseudo_hap_2\n{pseudo_hap_2}\n')

    logging.info('Pseudo-haplotypes created successfully')
    return best_haps_list

def calculate_coverage(k, genome_len, read_mers_file):
    logging.info('Calculating coverage')
    cmd = f'jellyfish stats {read_mers_file} > stats.txt'
    subprocess.run(cmd, shell=True)
    with open('stats.txt', 'r') as file:
        f = file.readlines()
        tot_k = int(f[2].split()[-1].strip())
    coverage = tot_k / genome_len
    logging.info(f'Coverage calculated: {coverage}')
    return coverage

def QV(hap_iterable, reads_query_table, coverage, k):
    logging.info('Calculating Quality Value (QV)')
    errors = 0 
    tot_hap_counts = 0
    
    hap_mers = load_iterable_jf(hap_iterable)
    
    for hap_mer, hap_mer_count in hap_mers:
        tot_hap_counts += hap_mer_count
        read_mer_count_raw = query_hash_table(reads_query_table, hap_mer) 
        
        if read_mer_count_raw:
            read_mer_count = max(int((read_mer_count_raw / coverage)*2), 1)
        
        if not read_mer_count:
            errors += hap_mer_count
        elif hap_mer_count > read_mer_count:
            errors += (hap_mer_count - read_mer_count)

    errors += 0.0000000000001
    err_probability = 1 - (1 - (errors / tot_hap_counts))**(1 / k)
    qv = (-10) * math.log10(err_probability)

    with open('QV_results.txt', 'w') as qfile:
        qfile.write(f'Erroneous kmers: {errors}\n')
        qfile.write(f'Total hap kmers: {tot_hap_counts}\n')
        qfile.write(f'QV value: {qv}')
    
    logging.info(f'Quality Value (QV) calculated: {qv}')
    return qv

def main(args):
    log_file = args.log_file if args.log_file else 'process.log'
    setup_logging(log_file)
    logging.info('Starting main process')

    k = args.kmer_length
    
    ref, genome_len = load_ref(args.reference)
    
    run_jf_count(args.reads, k, args.output + 'reads_mer_counts.jf', args.threads, args.min_count)
    read_mers = load_query_table(args.output + 'reads_mer_counts.jf')
    read_mers_coverage = calculate_coverage(k, genome_len, args.output + 'reads_mer_counts.jf')
    read_mers_coverage = 3.92
    clusters = cluster_variants(args.vcf, args.sample, k, args.output)
    
    best_haps = create_pseudo_haplotypes(ref, clusters, read_mers, args.output + 'pseudo_haps.fasta')
    run_jf_count(args.output + 'pseudo_haps.fasta', k, args.output + 'hap_mer_counts.jf', args.threads, 1)
    
    qv = QV(args.output + 'hap_mer_counts.jf', read_mers, read_mers_coverage, k)
    
    logging.info(f'Process completed. Quality Value: {qv}')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="VCQC Beta Command Line Tool")
    parser.add_argument('--threads', type=int, required=True, help="Number of threads to use in Jellyfish count")
    parser.add_argument('--min_count', type=int, required=True, help="Minimum count parameter for Jellyfish count")
    parser.add_argument('--output', type=str, required=True, help="Output path prefix for produced files")
    parser.add_argument('--vcf', type=str, required=True, help="VCF file")
    parser.add_argument('--reads', type=str, required=True, help="Reads file")
    parser.add_argument('--reference', type=str, required=True, help="Reference genome file")
    parser.add_argument('--sample', type=str, required=True, help="Sample name")
    parser.add_argument('--kmer_length', type=int, required=True, help="K-mer length")
    parser.add_argument('--log_file', type=str, help="Log file name (optional)")

    args = parser.parse_args()
    main(args)

