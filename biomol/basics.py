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