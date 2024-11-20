from typing import List

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

