from typing import List
from typing import Tuple

def count_dna_nucleotides(sequence: str) -> List[int]:
    """
    https://rosalind.info/problems/dna/
    
    Given: A DNA string of length at most 1000 nt.
    Return: Four integers counting the respective number of times that the symbols 'A', 'C', 'G', and 'T' occur.
    """
    bases = {'A': 0, 'C': 0, 'G': 0, 'T':0}
    for s in sequence:
        bases[s] += 1
        
    return bases

def transcribe_dna_into_rna(sequence: str) -> str:
    """
    An RNA string is a string formed from the alphabet containing 'A', 'C', 'G', and 'U'.

    Given a DNA string t corresponding to a coding strand, its transcribed RNA string u
    is formed by replacing all occurrences of 'T' in t with 'U' in u.

    Given: A DNA string t having length at most 1000 nt.
    Return: The transcribed RNA string of t.
    """
    return ''.join([s if s != 'T' else 'U' for s in sequence])

def reverse_complement_dna(sequence: str) -> str:
    """
    In DNA strings, symbols 'A' and 'T' are complements of each other, as are 'C' and 'G'.

    The reverse complement of a DNA string s is the string sc formed by reversing the symbols of s,
    then taking the complement of each symbol (e.g., the reverse complement of "GTCA" is "TGAC").

    Given: A DNA string s of length at most 1000 bp.
    Return: The reverse complement sc of s.
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join([complement[s] for s in sequence][::-1])

def rabbits_and_recurrence_relations(months: int, multiplier: int) -> int:
    """
    A sequence is an ordered collection of objects (usually numbers), which are allowed to repeat. 
    Sequences can be finite or infinite. Two examples are the finite sequence (π,−2‾√,0,π) 
    and the infinite sequence of odd numbers (1,3,5,7,9,…). 
    We use the notation an to represent the n-th term of a sequence.

    A recurrence relation is a way of defining the terms of a sequence 
    with respect to the values of previous terms. 
    In the case of Fibonacci's rabbits from the introduction, 
    any given month will contain the rabbits that were alive the previous month, 
    plus any new offspring. A key observation is that the number of offspring 
    in any month is equal to the number of rabbits that were alive two months prior. 
    As a result, if Fn represents the number of rabbit pairs alive after the n-th month,
    then we obtain the Fibonacci sequence having terms Fn that are defined by 
    the recurrence relation Fn=Fn−1+Fn−2(with F1=F2=1 to initiate the sequence). 
    Although the sequence bears Fibonacci's name, it was known to Indian mathematicians over two millennia ago.

    When finding the n-th term of a sequence defined by a recurrence relation, 
    we can simply use the recurrence relation to generate terms for progressively larger values of n. 
    This problem introduces us to the computational technique of dynamic programming, 
    which successively builds up solutions by using the answers to smaller cases.

    Given: Positive integers n≤40 and k≤5.
    Return: The total number of rabbit pairs that will be present after n months, 
    if we begin with 1 pair and in each generation, every pair of reproduction-age rabbits 
    produces a litter of k rabbit pairs (instead of only 1 pair).
    """
    if months <= 2:
        return 1
    else:
        return rabbits_and_recurrence_relations(months - 1, multiplier) + (rabbits_and_recurrence_relations(months - 2, multiplier) * multiplier)


def fasta_to_dict(path_to_fasta: str) -> dict:
    parsed_fasta = {}
    current_key = ""
    try:
        with open(path_to_fasta, 'r') as fasta:
            for line in fasta:
                line = line.replace("\n", "") 
                if line[0] == ">":
                    current_key = f"{line[1:]}" 
                    # parsed_fasta[current_key] = []
                    parsed_fasta[current_key] = ""
                else:
                    # If each line has multiple pieces of info
                    # parsed_fasta[current_key].append(line)

                    # If each line is a continuation of the sequence
                    parsed_fasta[current_key] = parsed_fasta[current_key] + line
    except FileNotFoundError:
        print("File not found.")
    
    return parsed_fasta
    


def computing_gc_content(path_to_fasta: str) -> Tuple[str, float]: 
    """
    The GC-content of a DNA string is given by the percentage 
    of symbols in the string that are 'C' or 'G'. 
    For example, the GC-content of "AGCTATAG" is 37.5%. 
    Note that the reverse complement of any DNA string has the same GC-content.

    DNA strings must be labeled when they are consolidated into a database. 
    A commonly used method of string labeling is called FASTA format. 
    In this format, the string is introduced by a line that begins with '>', 
    followed by some labeling information. Subsequent lines contain the string itself; 
    the first line to begin with '>' indicates the label of the next string.

    In Rosalind's implementation, a string in FASTA format will be labeled by the ID "Rosalind_xxxx", 
    where "xxxx" denotes a four-digit code between 0000 and 9999.

    Given: At most 10 DNA strings in FASTA format (of length at most 1 kbp each).
    Return: The ID of the string having the highest GC-content, 
    followed by the GC-content of that string. 
    Rosalind allows for a default error of 0.001 in all decimal answers unless otherwise stated; 
    please see the note on absolute error below.
    """
    fasta_dict = fasta_to_dict(path_to_fasta)
    total_count, gc_count = 0, 0
    highest_gc_name = ""
    highest_gc_percentage = 0
    
    for name, sequence in fasta_dict.items():
        for base in sequence:
            if base == 'G' or base == 'C':
                gc_count += 1
            total_count += 1
        gc_percentage = gc_count / total_count
        if gc_percentage > highest_gc_percentage:
            highest_gc_percentage = gc_percentage
            highest_gc_name = name
        gc_count, total_count = 0, 0
    
    return highest_gc_name, highest_gc_percentage * 100

def count_point_mutations(sequence_1: str, sequence_2: str) -> int:
    """
    Given two strings s and t of equal length, 
    the Hamming distance between s and t, denoted dH(s,t), 
    is the number of corresponding symbols that differ in s and t.

    Given: Two DNA strings s and t of equal length (not exceeding 1 kbp).
    Return: The Hamming distance dH(s,t).
    """
    hamming_distance = 0
    for index in range(0, len(sequence_1)):
        if sequence_1[index] != sequence_2[index]:
            hamming_distance += 1
            
    return hamming_distance
