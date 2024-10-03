import random
import gzip
from Bio import SeqIO
from os import system 

def simulate_genome(reference_file, output_vcf, output_reads1, output_reads2, output_haplotype1, output_haplotype2, coverage, error_rate, k, linked_variant_frequency, seed=None):
    # Set random seed for reproducibility
    random.seed(seed)
    
    # Load reference genome
    reference = load_reference(reference_file)
    print('Reference is loaded')
    
    # Simulate variations (SNPs, insertions, deletions)
    variants, linked_variants = simulate_variations(reference, linked_variant_frequency, k)
    print('Variants are produced')
    
    # Write VCF file
    write_vcf(output_vcf, variants,reference)
    print('VCF file is generated')
    
    # Extract haplotype sequences
    haplotype1_sequence = extract_haplotype(reference, variants, genotype='1|0')
    haplotype2_sequence = extract_haplotype(reference, variants, genotype='0|1')
    print('The two haplotypes are generated')
    
    
    # Generate reads for haplotype 1 and write to file
    write_reads(output_reads1, generate_reads(haplotype1_sequence, coverage, error_rate))
    print('Reads for hap1 are generated')
    
    # Generate reads for haplotype 2 and write to file
    write_reads(output_reads2, generate_reads(haplotype2_sequence, coverage, error_rate))
    print('Reads for hap2 are generated')
    
    # Write haplotype sequences to compressed FASTA files
    write_haplotype(output_haplotype1, haplotype1_sequence)
    write_haplotype(output_haplotype2, haplotype2_sequence)
    print('Haplotypes are saved')
    
    return variants, linked_variants, haplotype1_sequence, haplotype2_sequence


def load_reference(reference_file):
    """
    Load reference genome from a FASTA file.

    Parameters:
        reference_file (str): Path to the reference genome file in FASTA format.

    Returns:
        dict: A dictionary containing the reference genome with chromosome names as keys and sequences as values.
    """
    reference = {}
    handle = SeqIO.parse(open(reference_file),'fasta')
    
    for fasta in handle:
        reference[fasta.id] = str(fasta.seq)
    return reference

import random

# TODO: allow also DELS
def generate_linked_variants(reference, k, curr_pos):
    # Generate linked variants within a distance k
    linked_variants = []
    chr_name, chr_sequence = list(reference.items())[0]  # Assume there is only one chromosome
    prev_position = curr_pos + k + 1
    last_deleted_position = -1
    num_variants = random.randint(2, 5) # TODO: get it up to 40 again
    for _ in range(num_variants):
        # Choose position within the allowed range
        max_position = min(len(chr_sequence) - 1, prev_position + k - 1)
        position = random.randint(prev_position, max_position)
            
        ref_base = chr_sequence[position]
        
        if ref_base == 'N':
            prev_position +=1
            continue
            
        variant_type = random.choice(['SNP', 'INS'])  # Removed 'DEL' from variant_type choices
        if variant_type == 'SNP':
            alt_base = random.choice(['A', 'C', 'G', 'T'])
            while alt_base == ref_base:
                alt_base = random.choice(['A', 'C', 'G', 'T'])
        elif variant_type == 'INS':
            alt_len = random.randint(2, 5) # TODO: fix to accept ins on len 1
            alt_base = chr_sequence[position] + ''.join(random.choice(['A', 'C', 'G', 'T']) for _ in range(alt_len))
        genotype = random.choice(['1|0', '0|1', '1|1'])  # Randomly assign phased genotype
        linked_variants.append((chr_name, position, ref_base, alt_base, variant_type, genotype))
        prev_position = position + 1
        
    return linked_variants


# TODO: insert the case 1/2
def simulate_variations(reference, linked_variant_frequency, k):
    # Simulate variations (SNPs, insertions, deletions)
    variants = []
    linked_variants = []
    last_linked_position = -1  # Initialize the last position of the last linked variant
    last_deleted_position = -1  # Initialize the ending position of the last deletion
    prev_pos = -1
    for chr_name, chr_sequence in reference.items():
        for position in range(len(chr_sequence)):
            # Check if position is after the last linked variant position
            if position <= last_linked_position + k + 2:
                continue
            
            # Skip variant simulation if reference base is 'N'
            if chr_sequence[position] == 'N':
                continue
            
            # Skip variant simulation if position is within the deleted portion
            if position <= last_deleted_position:
                continue
            
            if position < prev_pos + k + 2:
                continue
            
            # Choose which type of variant to simulate for this position
            variant_type = random.choice(['linked_variants', 'SNP', 'INS', 'DEL'])
  
            if variant_type == 'linked_variants' and random.random() < linked_variant_frequency:
                linked_var = generate_linked_variants({chr_name: chr_sequence}, k, position)
                if linked_var:
                    linked_variants.append(linked_var)
                    last_linked_position = max(variant[1] for variant in linked_var)
                    variants.extend(linked_var)
                
            elif variant_type == 'SNP' and random.random() < 0.01:
                ref_base = chr_sequence[position]
                alt_base = random.choice(['A', 'C', 'G', 'T'])
                while alt_base == ref_base:
                    alt_base = random.choice(['A', 'C', 'G', 'T'])
                genotype = random.choice(['1|0', '0|1', '1|1'])  # Randomly assign phased genotype
                variants.append((chr_name, position, ref_base, alt_base, 'SNP', genotype))
                prev_pos = position
                
            elif variant_type in ['INS', 'DEL'] and random.random() < 0.01: #0.00001
                if variant_type == 'INS':
                    ref_base = chr_sequence[position]
                    alt_len = random.randint(1, 5) # TODO: fix to accept ins on len 1
                    alt_base = chr_sequence[position] + ''.join(random.choice(['A', 'C', 'G', 'T']) for _ in range(alt_len))
                else:  # Deletion
                    alt_len = random.randint(1, 5)
                    alt_base = chr_sequence[position - 1]  # Alternate allele is the base before the deletion
                    ref_base = chr_sequence[position:position + alt_len]  # Use the deleted portion as reference allele
                    position -= 1  # Adjust position to the starting position of the deletion
                    last_deleted_position = position + alt_len - 1  # Update the ending position of the last deletion
                genotype = random.choice(['1|0', '0|1', '1|1'])  # Randomly assign phased genotype
                variants.append((chr_name, position, ref_base, alt_base, variant_type, genotype))
                prev_pos = position
    return variants, linked_variants


# TODO: add REVCOMP reads
def generate_reads(reference, coverage, error_rate, read_length=150):
    # Generate reads randomly choosing starting position on the genome
    # For simplicity, equally sample from two haplotypes
    len_genome = 0
    for chrom in reference.keys():
        len_genome += len(reference[chrom])
    
    n_reads = int((coverage * len_genome) / read_length)
    
    for _ in range(n_reads):
        chr_name = random.choice(list(reference.keys()))
        chr_sequence = reference[chr_name]
        start_pos = random.randint(0, len(chr_sequence) - read_length)
        #start_pos = random.randint(0, len(chr_sequence) - read_length*2)  # Ensure the read fits within the sequence
        end_pos = start_pos + read_length
        read = chr_sequence[start_pos:end_pos]
        # Introduce errors based on error rate
        for i in range(len(read)):
            if random.random() < error_rate:
                if random.random() < 0.5:  # 50% chance of introducing a small insertion or deletion
                    if random.random() < 0.5:  # 50% chance of insertion
                        insertion_len = random.randint(1, 5)
                        insertion = ''.join(random.choice(['A', 'C', 'G', 'T']) for _ in range(insertion_len))
                        read = read[:i] + insertion + read[i:]
                        
                    else:  # Deletion
                        deletion_len = random.randint(1, 5)  
                        read = read[:i] + read[i + deletion_len:] 

                else:  # Base substitution
                    read = read[:i] + random.choice(['A', 'C', 'G', 'T']) + read[i+1:]
        if len(read) < read_length:
            actual_length = len(read)
            read = read + chr_sequence[end_pos: read_length - actual_length] # add missing ref bases to reach the desired read_length
        elif len(read) > read_length:
            read = read[0:read_length]  # trimm bases exceeding the desired read_length
        yield chr_name, start_pos, read


def write_vcf(output_vcf, variants, reference):
    # Open the output VCF file for writing
    with open(output_vcf, 'w') as vcf_file:
        # Write header lines
        vcf_file.write('##fileformat=VCFv4.3\n')
        
        # Write contig definitions
        for chr_name, chr_length in reference.items():
            vcf_file.write(f'##contig=<ID={chr_name},length={len(chr_length)}>\n')
        
        # Write INFO and FORMAT meta-information lines
        vcf_file.write('##INFO=<ID=TYPE,Number=1,Type=String,Description="Type of variation">\n')
        vcf_file.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        
        # Write column headers
        vcf_file.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n')
        
        # Write variant records
        for chr_name, position, ref, alt, variation_type, genotype in variants:
            # Adjust position to 1-based
            position += 1
            
            if genotype == '1|0' or genotype == '0|1':
                genotype = '0/1'
            elif genotype == '1|1':
                genotype = '1/1'
            
            # Format the variant record
            record_fields = [
                chr_name,
                str(position),
                '.',  # ID field
                ref,
                alt,
                '.',  # QUAL field
                '.',  # FILTER field
                f'TYPE={variation_type}',  # INFO field
                'GT',  # FORMAT field
                genotype  # Sample genotype
            ]
            
            # Write the variant record
            vcf_file.write('\t'.join(record_fields) + '\n')

def write_reads(output_reads, reads_generator):
    # Write reads to a file
    with open(output_reads, 'w') as f:
        for chr_name, position, read in reads_generator:
            f.write(f'>{chr_name}:{position}\n')
            f.write(f'{read}\n')

def extract_haplotype(reference, variants, genotype):
    # Extract haplotype sequence based on genotype
    haplotype_sequence = {}
    for chr_name, chr_sequence in reference.items():
        haplotype_sequence[chr_name] = list(chr_sequence)
    for chr_name, position, ref, alt, variation_type, var_genotype in variants:
        if var_genotype == genotype or var_genotype == '1|1':
            len_ref = len(ref)
            haplotype_sequence[chr_name][position:position+len_ref] = ['']*len_ref
            haplotype_sequence[chr_name][position] = alt  # Use the last base of the alt allele
    for chr_name in haplotype_sequence:
        haplotype_sequence[chr_name] = ''.join(haplotype_sequence[chr_name])
    return haplotype_sequence



def write_haplotype(output_haplotype, haplotype_sequence):
    # Write haplotype sequence to a compressed FASTA file
    #with gzip.open(output_haplotype, 'wt') as f:
    with open(output_haplotype, 'w') as f:
        for chr_name, sequence in haplotype_sequence.items():
            f.write(f'>{chr_name}\n')
            f.write(f'{sequence}\n')


# Example usage:
reference_file = 'mini-chr20.fasta'
output_vcf = 'ground_truth/output.vcf'
output_reads1 = 'ground_truth/reads_haplotype1.fasta'
output_reads2 = 'ground_truth/reads_haplotype2.fasta'
output_haplotype1 = 'ground_truth/haplotype1.fasta'
output_haplotype2 = 'ground_truth/haplotype2.fasta'
coverage = 200
error_rate = 0
k = 21  # Maximum distance between linked variants
linked_variant_frequency = 0.01  # Frequency of generating linked variants
seed = 1246  # SEED USATO PER TUTTO
#seed = 3333  # Set seed for reproducibility
#seed = 666
variants, linked_variants, haplotype1_sequence, haplotype2_sequence = simulate_genome(reference_file, output_vcf, output_reads1, output_reads2, output_haplotype1, output_haplotype2, coverage, error_rate, k, linked_variant_frequency, seed)

# Merge the two read sets
merge_reads_cmd = 'cat ' + output_reads1 + ' ' + output_reads2 + ' > ground_truth/merged_reads.fasta'
system(merge_reads_cmd) 
# Merge the two haplotypes
merge_hap_cmd = 'cat ' + output_haplotype1 + ' ' + output_haplotype2 + ' > ground_truth/merged_haps.fasta'
system(merge_hap_cmd) 

print("Generated", len(variants), "variants")
