"""
Transcribing DNA into RNA.

An RNA string is formed from the alphabet containing A, C, G, and U.
Given a DNA string s corresponding to a coding strand, its transcribed 
RNA string u is formed by replacing all occurances of T in s with U in U.
"""
def transcribe(dna_string):
    return dna_string.replace('T', 'U')

print(transcribe('GATGGAACTTGACTACGTAAATT'))