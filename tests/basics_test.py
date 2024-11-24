import pytest
from biomol import basics

# Count DNA nucleotide
@pytest.mark.parametrize("input, expected_output", [
    ("AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC", {'A': 20, 'C': 12, 'G': 17, 'T':21}),
    ("GTACTTAACCGCAAAGGGATTAACAACACTGCGAATGATCAGCGCTCCCACAGGCCGACCTTCCGTGCACCCAGTCATTTGCTCCGCACTACAAACATGTCAGGACCACAAAAACGGGTAGTTTCCCAACGAATCGACAGCGAAGCGGGCTCGACGGCCTCCCGCACACAAAGCTACTCCGTATCTACTGCGCCCTGGGAGAACTATTCCGCAGAGTTAGCGCAGGCCAGTGTAGTTGCCCCATAATGTGCGTGTAGGGACGCAATCAACCCTGACAAGGCAGACATAGGCACGACTCGCGAGCCTGTCATTCAATAACGATGTTTTAAATTCGACTGGTTGCCTTCATTCCTTCACAGGTCTACTCTGCGATTACCACCCTGGTGACAAGTAGCCCGCCCATAGTATCGTCCCTATCTGCCTCCGACTGGATAAGCATACCTGGCGGCCCCTACCAGCGATGACCGAACCTCACTTTGTATAGATCTCACTCCGCGTCAACTAGATGCCTCCTGGTTAGCCCAGTACCCCAGTCTGAAGCCCATTTAATGGCCTCGTTGAGTCAGGCTCATCACTCCCACATTGGACGCATTGGATAGCAGTACGTTCTAACGCGCCCAGTTTATTAAAAATTGGCAGATGTACCAACACAGCTGCAAAAGTACCCCCTTATTCAGATTGAGATTTCTGGTTACACAGTAGGTGACAGGTAAAACGAAATCTCTCTTCGCGTTATCCGGACGTCGGAGCGTGACGCCATCTATAGAGGCTGCGAACTGCGCCATCAAGCGTTGGCAGTCCCTCTTTGGAGCCCTCAGGCATACTCGATCCATAAATCTATTATGGATCAGGACCAGGATACGAGCCTCCGGCCAGTGCCGGCATGAATGTGGTGAAGGTATGCGTCTGGTCATCTCACAACTCTAATGGTCGGGGGTGGTAGAAAGAGCGCCTG", {'A': 239, 'C': 276, 'G': 224, 'T':216})
])
def test_count_dna_nucleotides(input, expected_output):
    assert basics.count_dna_nucleotides(input) == expected_output
    

# Transcribe DNA into RNA
@pytest.mark.parametrize("input, expected_output", [
    ("GATGGAACTTGACTACGTAAATT", "GAUGGAACUUGACUACGUAAAUU"),
    ("TAGACCGCAAATGTAGCATGATACTACTTTGGCGCATAACACTCCGCTCCTAGCTGGCTTAGTCCTTATTCTTGTCAGGCTAGTGTGGATTATCTGCTCTGGCTGCAAAGCCGTCATCCAAGTTGCGCACGTTACCTACGGGTTTTCATAGCCACCTCACCATGAGTGACCCGTGGTGGCTGTTTTAACCGAAACAATCAGATGGTAGTTTCTCCAATGGTATGTACCTTAGCCTTGCAGATCACAAAGGCACTCGCTATCACATACGGACTTAAGAAACCTTACCCATCTCGGAGGAGTGATTCTTGGTCTGGGATATCGCAGAACATTCATAGGATTCCACCTCCCCCATAAATACAAATAAGGCCGATAGACAAGGAGATCTGTCGTAGGTGATATCACCCGAAGAGATGCTCGCATACATAGCGCTAATCTCAAGGATGACGAGGTCCCTGGAACCATATCCCACATGCCTGTCGAAGCTACGAGAATACCTTACAGCAAATATTAGACACCCAGCCGCCGACTAGTACGAGCATAGGGGTTTCGGCAACTCTCTGTTTACGAGATATTTAAGTTTCTTATTCTGCGTCGGCAGATCCGGGGCGAGACGAGTGCTTGCTCAGGGTCAAGAGGCATGAGATTATTCACCTAGGATTAGACGATACACAGAATGATTGACCATTACCCAACGCATTCTAATTAGAACCTGTGATACAGGTGCCAACAAGGGCATTAGGTGTGAGCCGGCCTCGCAGGCTCTACTTTTACCGGTATGCAGGCGAGAACCTAGCAGGCGATTAGCACGTATCGTGTAACCCTGGAGAGCCCCTAGTCCTGGTGGAAAAAACTTACACGGGGCCACTAGTAGAACTGTATTGATAGGTGACAGATGAACTTAGCTTGCAATGGCAGGCCGAAGTAGGTACCCTTCCCGTG",
     "UAGACCGCAAAUGUAGCAUGAUACUACUUUGGCGCAUAACACUCCGCUCCUAGCUGGCUUAGUCCUUAUUCUUGUCAGGCUAGUGUGGAUUAUCUGCUCUGGCUGCAAAGCCGUCAUCCAAGUUGCGCACGUUACCUACGGGUUUUCAUAGCCACCUCACCAUGAGUGACCCGUGGUGGCUGUUUUAACCGAAACAAUCAGAUGGUAGUUUCUCCAAUGGUAUGUACCUUAGCCUUGCAGAUCACAAAGGCACUCGCUAUCACAUACGGACUUAAGAAACCUUACCCAUCUCGGAGGAGUGAUUCUUGGUCUGGGAUAUCGCAGAACAUUCAUAGGAUUCCACCUCCCCCAUAAAUACAAAUAAGGCCGAUAGACAAGGAGAUCUGUCGUAGGUGAUAUCACCCGAAGAGAUGCUCGCAUACAUAGCGCUAAUCUCAAGGAUGACGAGGUCCCUGGAACCAUAUCCCACAUGCCUGUCGAAGCUACGAGAAUACCUUACAGCAAAUAUUAGACACCCAGCCGCCGACUAGUACGAGCAUAGGGGUUUCGGCAACUCUCUGUUUACGAGAUAUUUAAGUUUCUUAUUCUGCGUCGGCAGAUCCGGGGCGAGACGAGUGCUUGCUCAGGGUCAAGAGGCAUGAGAUUAUUCACCUAGGAUUAGACGAUACACAGAAUGAUUGACCAUUACCCAACGCAUUCUAAUUAGAACCUGUGAUACAGGUGCCAACAAGGGCAUUAGGUGUGAGCCGGCCUCGCAGGCUCUACUUUUACCGGUAUGCAGGCGAGAACCUAGCAGGCGAUUAGCACGUAUCGUGUAACCCUGGAGAGCCCCUAGUCCUGGUGGAAAAAACUUACACGGGGCCACUAGUAGAACUGUAUUGAUAGGUGACAGAUGAACUUAGCUUGCAAUGGCAGGCCGAAGUAGGUACCCUUCCCGUG")
])
def test_transcribe_dna_into_rna(input, expected_output):
    assert basics.transcribe_dna_into_rna(input) == expected_output


# Reverse complement of DNA
@pytest.mark.parametrize("input, expected_output", [
    ("AAAACCCGGT", "ACCGGGTTTT"),
    ("TTTGGAAGCTTTCCTACGCGTGCACTGGTGGACGCCATTCTCCCTACAAGGTAGTGCCGGAAACCCTGGCGCACTAGCAAAATCGCAATGCTAACAAGGTGACTTGTGAGATCCGCGGTTAATTACCTAGGATCATCCTCGTGTGCCGTATTTTAGGGCCACCAGGTAGCTTGGACCCCTTTTCTTGTCCCATTTGGGCATCGAGGAATGCCCGATTCATAGGTGTAAACCGGCTGGTCTGGCCAGTAATGTCTTAATCTATCCTACCCCTTTTACTCCACTCTTATGGCCGGAGAGGTACCCGTGCGCCCGCCCCGTTCTGGGCTTCAGTTGGCTCTAATCGCACACGGAAGCCCCCCGAAGATGTGCCGCTCTCGGTTGGCTACGAGGGGAGGCCCCAGAATGGTCATATCTTCGCTCGCGTGAATTCGCCTGGAATAGGGCGGGTCCCCGAGAGGTTGCCTTCATGGAAAGTCGCGCATGCCATCCAGGCAGGGGTCCCGTGGATCGCAGAACAGAACGTGTTCCCAGTCTGATCCTTCCGAATAAGGAAGATCAGGTCCACATACACTGTTGGTCGCTTTTACGGAGCATGTAGTCCCAGCTAGCCCCACTGAGCTCGCTCAAGTGCAACTACGTGAGGTACGGCCGATTCCATTGCACTACGAATACACCGGCGGTCTTCCATACATATTGGCCGAAATATTAATATGAGCTAAGTCTACGGAATCATAGCCCGTGGTGTTTCCGGTGAAGCAGGGGAGTGGCGCGCTGTTTGTTATGTCGGAGGGTAGACCATGTGGAAGAGACGCATATGTGGGTTGAACCCATACGCGCCTTAGTCCATGTGGCTACTCCCGGGGACAGGTCCTGG",
     "CCAGGACCTGTCCCCGGGAGTAGCCACATGGACTAAGGCGCGTATGGGTTCAACCCACATATGCGTCTCTTCCACATGGTCTACCCTCCGACATAACAAACAGCGCGCCACTCCCCTGCTTCACCGGAAACACCACGGGCTATGATTCCGTAGACTTAGCTCATATTAATATTTCGGCCAATATGTATGGAAGACCGCCGGTGTATTCGTAGTGCAATGGAATCGGCCGTACCTCACGTAGTTGCACTTGAGCGAGCTCAGTGGGGCTAGCTGGGACTACATGCTCCGTAAAAGCGACCAACAGTGTATGTGGACCTGATCTTCCTTATTCGGAAGGATCAGACTGGGAACACGTTCTGTTCTGCGATCCACGGGACCCCTGCCTGGATGGCATGCGCGACTTTCCATGAAGGCAACCTCTCGGGGACCCGCCCTATTCCAGGCGAATTCACGCGAGCGAAGATATGACCATTCTGGGGCCTCCCCTCGTAGCCAACCGAGAGCGGCACATCTTCGGGGGGCTTCCGTGTGCGATTAGAGCCAACTGAAGCCCAGAACGGGGCGGGCGCACGGGTACCTCTCCGGCCATAAGAGTGGAGTAAAAGGGGTAGGATAGATTAAGACATTACTGGCCAGACCAGCCGGTTTACACCTATGAATCGGGCATTCCTCGATGCCCAAATGGGACAAGAAAAGGGGTCCAAGCTACCTGGTGGCCCTAAAATACGGCACACGAGGATGATCCTAGGTAATTAACCGCGGATCTCACAAGTCACCTTGTTAGCATTGCGATTTTGCTAGTGCGCCAGGGTTTCCGGCACTACCTTGTAGGGAGAATGGCGTCCACCAGTGCACGCGTAGGAAAGCTTCCAAA")
])
def test_reverse_complement_dna(input, expected_output):
    assert basics.reverse_complement_dna(input) == expected_output
    

# Rabbits and recurrence relations
@pytest.mark.parametrize("input, expected_output", [
    ((5, 3), 19),
    ((28, 3), 3855438727)
])
def test_rabbits_and_recurrence_relations(input, expected_output):
    assert basics.rabbits_and_recurrence_relations(*input) == expected_output
    
    
# Computing GC content
@pytest.mark.parametrize("input, expected_output", [
    ("tests/test_cases/test.txt", ('Rosalind_0808', 60.91954022988506)),
    ("tests/test_cases/rosalind_gc.txt", ('Rosalind_3456', 52.430196483971045))
])
def test_computing_gc_content(input, expected_output):
    assert basics.computing_gc_content(input) == expected_output


# Count point mutations
@pytest.mark.parametrize("input, expected_output", [
    (("GAGCCTACTAACGGGAT", "CATCGTAATGACGGCCT"), 7),
    (("TACACTACTCCTAAAGCGATCCGCAGCGCTAGCGGTGGATTCATCTCGGTGCCAGTCAATGACGGGATTACGCCGCCATTGCTGGGATAATCCTCGCTCCACCACCGTAACGGGTGCCCCAACCGCGTGATCTATGACGCTCACTTGAAACTAGACGGCTCACTTTCACTCAAACGCTTCGTATCTCATTGTAGACCACATTTTCTCGACGCAGATATCTCGAGCTCCAGATAGTAGCGCTCGATCTCCACTCAATTCTTAGATCAGGCGAATAGCTTGTGGTAGCTGTCTTAACATCTCCTTACCCTAAAGGTATATCAGGCTTCCCTCGCAGCAGTGCTGGAGGACTTTGATACTGAGCAAAGGCTAGTCTAGAGCGGTACTCTTCGCACTGGTTCCGGATCTGTCACTTTCAATTAACAAGGACGGTATTTGACGCCTCTCTATCAATATGCAATAATGCCGCAAAAAGGGCTGCCGATGAGGTGAGCACTGGTTTGGTTGAAGTCTAACCACCGTCGCGCACCCCAGGTTAGGCTAGCTAATAATGGATGCGAGAGCAGGCTCCACGGCCATTCCGCCGTATGTTGGGATGTTGTTACATCAAATAGTGAGATTATCCTCACGCCCTGAAAAATCACGACCAAAACCTATTTGGGGCACCCGGTCTCTTTGCGGACTAGTGGACTTTGGTGGATTTACGATTAAATTTACACACCTGCATTAAGCGCATAAACTAATGTGGCGGCCAGTGTGTTTGGACACGATGATCCAGGTCTATGCCAAGCAGGTGGGAAGGACGGATGGTCGAGATCCCACTTCGGTCGATGAAGCACGATTTTCTATAAGCCCGAATACTATGGTCGGGGGTATTTCCGCAGGACCTTGGCTGACTACATGCTCCTAGTGTC", "TTCACCCGAACCAAGGCGAGCCGCAGTCTGTCGCGCCTGTCCGCCTCGGATCGCTTCATTTAGGGAAATTTTCAGCGTGTCCTGGAACTTTCGTTTTGCCGCCATCGTTGAGGCAGCCTAAAACGCCGAATCAGTGCGGCCCCAACGACGCTAGACAGTTCTTTATCGCGAGCTTGTAGTGCATTTGATTCTCGACCACCTGTACGTGCCGGACCGTGCTACTCTCCCAGCCTGTGATGTTCTTGGTGCGGGCTGTAAAGAGACTACACGACTATCTTATAGGCGAAGTATGATCCTCTGTACACCTCAAATCTTTCGCTGTTTATGTTCGCGGCAGGAAAGCCTGATTGTCGCAAAGAGGGAAATCTAGTCAAGACCTATTTTCACGGCAACGGGTAGGGCCTGGTCACTTTCCATTTCGCGGTTAAAGTTTGGGCTCCAAGAAAGGTGTGCCGAATTATCCCGCAAGAACATTATTCCACTTCGTAGAGTTTCACATAGTTAAAGGTCTCGCAGTTCCGTGCCTAGATGCGTAAGGAACCTGATGGGCTAACTTCGTGCAGGCAGAACCCGCAACGCGCCCAATTGCGAAATCCGATTTCCTCTACTACTGAGACTTAACTTACGCGCTGATAAGACTCCATCAATAGCTATGATAGCCGCCCGTTCGCATACCAAATTTGTTCAACTCCACGTATGTACGAATACAGTAACCGACCATACGTAAGGAAATCCACCAATGTGGAAGCCCGTGTGGTTGATGACCGCTATGCGCGAGTATGCAAATCCTGATTATAGGGCACCTGTTGCAGATCCCGCGTCTATCGTGCAAGTTAATTTGTGGGCAAGGCTCCATATTGCGTTGCGTGGTTTTCTCTCGCATACATTAGAGTCAACGTATTGTTAGTTCA"), 462)
])
def test_count_point_mutations(input, expected_output):
    assert basics.count_point_mutations(*input) == expected_output

    
# Mendel's first law
@pytest.mark.parametrize("input, expected_output", [
    ((2, 2, 2), 0.7833333333333333),
    ((18, 24, 22), 0.7202380952380952)
])
def test_(input, expected_output):
    assert basics.mendels_first_law(*input) == expected_output