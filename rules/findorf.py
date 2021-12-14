# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 11:57:20 2021

@author: cow082
"""

# standard library modules
import re

# 3rd party modules
from Bio import SeqIO
from Bio.Seq import Seq

__version__ = 1.0


def findORF(record):
    '''returns the longest ORF in a seqrecord
    '''
    
    def fixPartialCodonError(S):
        '''fixes this error'''
        while True:
            if len(S) % 3 != 0:
                S = S[:-1] # trim partial codon
            else:
                break
        return S        
    
    def removeGaps(S):
        '''avoid translation errors
        e.g. if seq contains codons such as --- or xxx '''
        return Seq(str(S).upper().replace('-', '').replace('X', 'N'))

    sequence = str(removeGaps(record.seq))
    record.seq
    
    if not 'ATG' in sequence:
        return '' # no ATG --> no ORF
    
    startP = re.compile('ATG')
    longest = (0,)
    for m in startP.finditer(sequence):
        candidate_orf = fixPartialCodonError(Seq(sequence)[m.start():])
        protein = candidate_orf.translate(to_stop=True, cds=False)
        if len(protein) > longest[0]:
            longest = (len(protein),
                       m.start(),
                       str(protein),
                       sequence[m.start():m.start()+len(protein)*3+3])
    
    cds_start = longest[1]
    cds_end = longest[1]+len(longest[2])*3
    predicted_ORF = record[cds_start:cds_end].seq
    
    return predicted_ORF


def test1():
    """
    """
    HOME = "C:/Users/cow082/aivpipe"
    fasta = '%s/NC_004910.1.fasta' % HOME
    record = SeqIO.read(fasta, 'fasta')
    gene = record.seq
    orf = findORF(record)
    pep = orf.translate()
    print("\n\n".join([str(s) for s in [gene, orf, pep]]))


def test2():
    """
    """
    HOME = "C:/Users/cow082/aivpipe"
    fasta = '%s/21-02023-01_H5N3_complete_genome.fasta' % HOME
    records = SeqIO.parse(fasta, 'fasta')
    for i, record in enumerate(records):
        gene = record.seq
        orf = findORF(record)
        pep = orf.translate()
        print(i, "\n\n".join([str(s) for s in [gene, orf, pep]])+"\n")
    

def test3():
    """
    """
    HOME = "C:/Users/cow082/aivpipe"
    fasta = '%s/21-02023-03_H7N8_partial_genome.fasta' % HOME
    records = SeqIO.parse(fasta, 'fasta')
    for i, record in enumerate(records):
        gene = record.seq
        orf = findORF(record)
        pep = orf.translate()
        print(i, "\n\n".join([str(s) for s in [gene, orf, pep]])+"\n")


def main():
#    test1()
#    test2()
    test3()


if __name__ == "__main__":
    main()
