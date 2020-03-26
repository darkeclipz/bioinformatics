from importlib import reload  
import numpy as np


def count_dna_nucleotides(dna_string):
    """
    Counts how many times the nucleotides A, C, G, and T occur in the DNA string.
    """
    D = {k:0 for k in ['A', 'C', 'G', 'T']}
    for nucleotide in dna_string:
        D[nucleotide] += 1
    return D


def transcribe(dna_string):
    """
    An RNA string is formed from the alphabet containing A, C, G, and U.
    Given a DNA string s corresponding to a coding strand, its transcribed 
    RNA string u is formed by replacing all occurances of T in s with U in U.
    """
    return dna_string.replace('T', 'U')


def complement(dna_string):
    """
    Returns the complement of a DNA base. Each base bonds to a base in the
    opposite strand. Adenine always bonds to thymine, and cytosine always
    bonds with guanine.
    """
    bonds = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join([bonds[nucleotide] for nucleotide in dna_string])


def gc_content(dna_string):
    """
    The GC-content of a DNA string is the percentage of symbols in the
    string that are C or G. If this is >50%, this suggests it might
    be prokaryotes.
    """
    D = count_dna_nucleotides(dna_string)
    return (D['G'] + D['C']) / len(dna_string)


def gc_content_fasta(filename):
    """
    Calculate the GC content for a list of DNA strings stored in FASTA format
    in a file.
    """
    f = open(filename, 'r')
    txt = f.read().split('>')[1:]
    S = [s.split('\n', 1) for s in txt]
    S = [('>'+name, dna_string.replace('\n', '')) for name, dna_string in S]
    return {name: gc_content(dna_string)*100 for name, dna_string in S}
    

def print_dictionary(D):
    """
    Prints a dictionary, for use in the console.
    """
    for k, v in D.items():
        print(k, v)
    print('Dictionary contains {} elements.'.format(len(D)))


def hamming_distance(s, t):
    """ 
    The Hamming distance between s and t, is the number of
    corresponding symbols that differ in s and t.
    """
    if len(s) != len(t):
        raise ValueError('strings must have equal length.')
    return sum([1 for i in range(len(s)) if s[i] != t[i]])


def pattern_count(text, pattern):
    """
    Count how many times a k-mer (string of length k) occurs
    in a text.
    """
    count = 0
    for i in range(len(text) - len(pattern)):
        if text[i:i+len(pattern)] == pattern:
            count += 1
    return count