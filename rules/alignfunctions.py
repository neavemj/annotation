# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 13:46:39 2021

@author: cow082
"""
# 3rd party modules
from Bio import pairwise2
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import substitution_matrices, MultipleSeqAlignment

# custom modules
import findorf

blosum62 = substitution_matrices.load('BLOSUM62')

def globalAlignment(refseq, testseq, nucleotide=False):
    '''Perform a global pairwise alignment of 2 seqRecords. 
    Should be equivalent to EMBOSS Needle
    '''
    try:
        seq1 = refseq
        seq2 = testseq
        orf1 = findorf.findORF(seq1)
        orf2 = findorf.findORF(seq2)
        pep1 = orf1.translate(to_stop=True)
        pep2 = orf2.translate(to_stop=True)
        
        if str(pep1) == str(pep2):
            print('testseq (%s) is identical to refseq (%s)' % (testseq.id, refseq.id))
#            return

        # perform global pairwise sequence alignment, equivalent to Needle
        try:
            alignments = pairwise2.align.globalds(
                    pep1, # refseq
                    pep2, # testseq
                    blosum62,  # substitution matrix: same as Needle default
                    -10,       # gapopen penalty: same as Needle default
                    -0.5,      # gapextend penalty: same as Needle default
                    penalize_end_gaps = (len(pep1) > 2000), # if seq1 < 2000, don't penalise end gaps!
                    )
        except SystemError:
            print('%s: Caused a SystemError' % testseq.id)
            return
        
        # format alignment as a string of 3 lines
        formatted = pairwise2.format_alignment(*alignments[0]).split('\n')
        
        # convert the aligned sequences back into seqrecords
        pro1 = SeqRecord(
            Seq(formatted[0]),
            id=seq1.id,
            name=seq1.name,
            description=seq1.description,
            )
        pro2 = SeqRecord(
            Seq(formatted[2]), 
            id=seq2.id,
            name=seq2.name,
            description=seq2.description,
            )
        
        # assemble into a Biopython MSA object
        align1 = MultipleSeqAlignment([pro1, pro2])

        # create gapped nucleotide sequences from aligned proteins
        gapped_orf1 = gapsFromPeptide(str(align1[0].seq), str(orf1))
        gapped_orf2 = gapsFromPeptide(str(align1[1].seq), str(orf2))
        
        # create SeqRecord objects from gapped nucleotide sequences
        record1 = SeqRecord(
            Seq(
                    gapped_orf1, 
#                    alphabet=Gapped(IUPAC.IUPACUnambiguousDNA()),
                    ),
            id=seq1.id,
            name=seq1.name,
            description=seq1.description,
            )
        record2 = SeqRecord(
            Seq(
                    gapped_orf2, 
#                    alphabet=Gapped(IUPAC.IUPACUnambiguousDNA()),
                    ),
            id=seq2.id,
            name=seq2.name,
            description=seq2.description,
            )
    
        # assemble gapped nucleotide sequences into a Biopython MSA object
        align2 = MultipleSeqAlignment([record1, record2])
    
        '''
        NOTE: codonalign confirms that my gapsFromPeptide function works 
        correctly. In future we could use codonalign to replace several steps, 
        since it returns an alignment of seqrecord objects. 
        
        NOTE: fails if sequences contain ambiguous bases eg Y
        '''
#        msa = codonalign.build(align1, [SeqRecord(s) for s in (orf1, orf2)])
#        assert msa[0].seq==record1.seq
#        assert msa[1].seq==record2.seq
#        print(msa.format("clustal"))

    except AssertionError:
        print('%s: AssertionError' % testseq.id)
        
    if nucleotide:
        return align2
    else: # peptide alignment
        return align1


def gapsFromPeptide(peptide_seq, nucleotide_seq):
    '''Transfers gaps from aligned peptide seq into codon partitioned 
    nucleotide seq (codon alignment)

    peptide_seq is an aligned peptide sequence with gaps that need to be 
    transferred to nucleotide_seq
    
    nucleotide_seq is an un-aligned dna sequence whose codons translate to 
    peptide_seq

    https://www.biostars.org/p/89741/
    
    # NOTE: codonalign does the same thing.
    '''

    def chunks(l, n):
        '''Yield successive n-sized chunks from l'''
        for i in range(0, len(l), n):
            yield l[i:i+n]

    codons = [codon for codon in chunks(nucleotide_seq, 3)]  #splits nucleotides into codons (triplets)
    gappedCodons = []
    codonCount = 0
    for aa in peptide_seq:  #adds '---' gaps to nucleotide seq corresponding to peptide
        if aa != '-':
            gappedCodons.append(codons[codonCount])
            codonCount += 1
        else:
            gappedCodons.append('---')
    return ''.join(gappedCodons)

