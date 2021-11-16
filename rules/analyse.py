# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 12:19:46 2021

@author: cow082
"""

# standard library modules
import subprocess
import io
import argparse
import os
import re

# 3rd party modules
from Bio import SeqIO, Entrez
import pandas as pd

# custom modules
import findorf
#import alignfunctions

__version__ = 1.0


def recordToFastaString(record):
    """convert a seqrecord to a fasta string
    """
    out_handle = io.StringIO()
    SeqIO.write(record, out_handle, "fasta")
    return out_handle.getvalue()


def blast(args, fasta):
    """
    """
    # http://www.metagenomics.wiki/tools/blast/blastn-output-format-6
    fields = 'qseqid sseqid qlen slen pident length mismatch gapopen qstart qend sstart send evalue bitscore'# qseq sseq'
    
    record = SeqIO.read(fasta, 'fasta')
    
    p = subprocess.run([
        args.program, 
        '-db', args.db, 
        '-max_hsps', '1', 
        '-strand', 'both', 
        '-outfmt', '6 %s' % fields,
            ], 
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE,
        input=recordToFastaString(record), 
        encoding='ascii',
            )

    try:
        assert p.returncode == 0
    except AssertionError:
        return 'BLAST error'
    
    data = pd.read_csv(
        io.StringIO(p.stdout),
        header=None,
        names=fields.split(),
        sep="\t",
        )
        
    try:
        top_hit = data.iloc[0:]
    except IndexError: # no hits at all
        return 'no BLAST hits'
    
    return data.iloc[:10] # top 10 hits


def findCleavageSite(s):
    """using regex, return the probable cleavage site sequence and its coordinates
    
    s <str> is an amino acid sequence
    """
    pattern = re.compile('N[V,I][H,R,P,L,Q].{1,15}G[L,I]F')
    match = pattern.search(s)
    start, end = match.span()
    if match:
        return match.group()[2:], start+2, end, len(s)


def analyse(args):
    """locate fasta files and run blastn and cleavage site analysis
    """
    # register with NCBI
    Entrez.email = args.email
    
    fasta_files = [os.path.join(root, f) for (root, dirs, files) in
            os.walk(os.path.abspath(args.input_dir)) 
            for f in files if f.endswith(".fasta")]
    
    for fasta in fasta_files:
        basename = os.path.basename(fasta).split(".")[0]
        segment = basename.split("_")[1]
        record = SeqIO.read(fasta, 'fasta')
        orf = findorf.findORF(record)
        peptide = str(orf.translate())
        
        # blast analysis
        blastresults = blast(args, fasta)
        blastresults.to_csv("%s/%s.blast.csv" % (args.output_dir, basename))
        
#        with Entrez.efetch(
#                db="nucleotide", rettype="gb", retmode="text", 
#                id=blastresults["sseqid"].iloc[0]) as handle:
#            top_hit = SeqIO.read(handle, "gb") 
#        
#        aln = alignfunctions(record, top_hit)
#        
#        with open("%s/%s.aln.fasta" % (args.output_dir, basename), 'w') as F:
#            SeqIO.write(aln, handle, "fasta")
        
        
        # cleavage site prediction
        if segment == "HA": # haemagglutinin
            try:
                cs , start, end, length = findCleavageSite(peptide)
                with open("%s/%s.cs.txt" % (args.output_dir, basename), 'w') as F:
                    F.write(" ".join([str(s) for s in [cs , start, end, length]]))
            except ValueError:
                print("cs not found")
        
        
        


def argParser():
    """
    """
    parser = argparse.ArgumentParser(description = "Performs a BLAST search")
    parser.add_argument('-input_dir', required=True, help = 'directory containing fasta files')
    parser.add_argument('-output_dir', required=True, help = 'where to save output?')
    parser.add_argument('-email', required=True, help = 'user email for NCBI queries')
    parser.add_argument('-program', default='tblastx', help = 'which blast program to run?')
    parser.add_argument('-db', default="influenza_A_genbank", help = 'which database to search?')
    return parser

            
def main():
    """
    """
    args = argParser().parse_args([
#            "-input_dir", "C:/Users/cow082/aivpipe/irma_assembly/datasets/21-02023-0001/FLU-avian-acdp/irma_output",
            "-input_dir", "C:/Users/cow082/aivpipe/irma_assembly/datasets/21-02023-0003/FLU-avian-acdp/irma_output",
            "-output_dir", "C:/Users/cow082/aivpipe/blastresults",
            "-email", "cow082@csiro.au",
            ])
    
    analyse(args)
    

if __name__ == "__main__":
    main()
