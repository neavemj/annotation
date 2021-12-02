# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 13:40:54 2021

@author: cow082
"""

# standard library modules
import subprocess
import io
import argparse
import os

# 3rd party modules
from Bio import SeqIO, Entrez
import pandas as pd

# custom modules
import findorf

__version__ = 1.0


def getRecord(email, accession):
    """retrieve a single GBK record in gbk format
    """
    Entrez.email = email
    
    return Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text").read()
    
    
def getLatestGenbankRecords(email, outdir, target=1000000, request_limit=10000):
    """routine to obtain all Influenza A records in GenBank (in batches of 10K)
    """
    Entrez.email = email

    # retrieve list of accession numbers
#    term = '"Influenza A virus"[porgn:__txid11320] AND "Australia"'
    term = '"Influenza A virus"[porgn:__txid11320]'
    handle = Entrez.esearch(db="nucleotide", retmax=int(1e6), term=term, idtype="acc") # default retmax is 20
    record = Entrez.read(handle)
    IdList = record["IdList"]
    count = len(IdList)
    handle.close()
    
    # retrive the desired number of records in batches
    target = min(target, count)
    query = IdList[:target]
    i = 0
    print("downloading %d of %d available records" % (target, count))
    
    while i < target:
        print([i,i+request_limit])
        handle = Entrez.efetch(db="nucleotide", id=query[i:i+request_limit], rettype="fasta", retmode="text")
        SeqIO.write(SeqIO.parse(handle, 'fasta'), 
                    '%s/genbank_influenza_A_11.10.2021_%d_%d.fasta' % (outdir, i+1, i+request_limit), "fasta")       
        handle.close()
        i += request_limit
#        time.sleep(0.35) # max 3 requests per second (not needed for batches as they take longer to process)
    
    return len(IdList)
    # as of 11.10.21 there are 18560  records for '"Influenza A virus"[porgn:__txid11320] AND "Australia"'
    # as of 11.10.21 there are 792486 records for '"Influenza A virus"[porgn:__txid11320]'


def countRecords():
    """verify the number of fasta records containe within the downloaded files
    """
    pwd = "C:/Users/cow082/aivpipe"
    starts = [i for i in range(0, 800000, 10000)]
    fasta_files = ["%s/genbank_influenza_A_11.10.2021_%d_%d.fasta" % (pwd, start+1, start+10000) for start in starts]
    count = 0
    for fasta in fasta_files:
        records = open(fasta, 'r').read().count('>')
        count += open(fasta, 'r').read().count('>')
        print(fasta, records)
        
    print(count)
    # 160000


def makeblastdb(fasta, title, dbtype="nucl", verbose=False):
    """
    to create:
        makeblastdb -in influenza_A_genbank.fasta -out influenza_A_genbank -title influenza_A_genbank -dbtype nucl -parse_seqids
    
    to validate:
        blastdbcmd -db influenza_A_genbank -info
    """   
    p = subprocess.run(
        [
            "makeblastdb",
            "-in", fasta,
            "-out", title,
            "-title", title,
            "-dbtype", dbtype,
            "-parse_seqids",
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, 
        encoding="ascii",
        )
    
    if verbose:
        print(io.StringIO(p.stdout).read())
        print(io.StringIO(p.stderr).read())
        
    return p.returncode


def blastRemote(record, program="blastp"):
    """seems quite slow to do an external query
    """
    args = [
        program, 
        "-db", "nr", 
        "-remote", 
        "-entrez_query", "Influenza A virus [Organism]",
        ]
        
    p = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    out, _ = p.communicate(record.encode())
    return out.decode()


def recordToFastaString(record):
    """convert a seqrecord to a fasta string
    """
    out_handle = io.StringIO()
    SeqIO.write(record, out_handle, "fasta")
    return out_handle.getvalue()


def runBlast(args):
    """perform a blast search
    """
    
    # http://www.metagenomics.wiki/tools/blast/blastn-output-format-6
    fields = 'qseqid sseqid qlen slen pident length mismatch gapopen qstart qend sstart send evalue bitscore'# qseq sseq'
    
    record = SeqIO.read(args.fasta, 'fasta')
    
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
    
#    print(io.StringIO(p.stdout).read())
#    print(io.StringIO(p.stderr).read())
    
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
    
#    NOTE: data['qseqid'] == record.id
    
    try:
        top_hit = data.iloc[0:]
    except IndexError: # no hits at all
        return 'no BLAST hits'
    
    return data.iloc[:10] # top 10 hits


def processName(s):
    """extract parts of the seq name
    
    returns:
        1. sample (usually the accession number)
        2. segment_number (for influenza this will be 1-8)
        3. segment_name (e.g. HA)
        4. segment_subtype (e.g. H7)
    """
    parts = s.split("_")
    
    label = "_".join(parts[:-1])
    seg_info = parts[-1].split("-")
    
    seg = seg_info[0]
    tail = seg.lstrip('0123456789')
    head = seg[:len(seg)-len(tail)]

    if len(seg_info) == 1:
        return label, head, tail, None
    elif len(seg_info) == 2:
        return label, head, tail, seg_info[1]
    else:
        print("processName encountered an unexpected format")


def test1():
    """retrieve one record from GenBank
    """
    record = getRecord(
        email="chris.cowled@csiro.au", 
        accession="gb|CY022693.1|",
        )
    print(record)


def test2():
    """retrieve all IAV records from GenBank
    """
    email="chris.cowled@csiro.au"
    pwd = "C:/Users/cow082/aivpipe"
    count, records = getLatestGenbankRecords(
        email=email, 
        outdir=pwd, 
        target=10, 
        request_limit=3,
        )
    for record in records:
        print(record)
    print(count)
    
   
def test3():
    """blast one record against nr using remote blast
    """
    pwd = "C:/Users/cow082/aivpipe"
    fasta = "%s/21-02023-03_H7N8_partial_genome.fasta" % pwd
    records = SeqIO.parse(fasta, 'fasta')
    for i, record in enumerate(records):
        blastresults = blastRemote(record)
        print(blastresults)
        
        
def test4():
    """blast multiple records against nr using remote blast
    """
    pwd = "C:/Users/cow082/aivpipe"
    fasta = "%s/21-02023-03_H7N8_partial_genome.fasta" % pwd
    records = SeqIO.parse(fasta, 'fasta')
    for i, record in enumerate(records):
        if i > 0:
            break
        orf = findorf.findORF(record)
        pep = orf.translate()
        blastresults = blastRemote(pep, program='blastp')
        print(blastresults)


def test5():
    """build a local database
    
    fasta is result of ncbi Nucleotide search: 
        "Influenza A virus"[porgn:__txid11320] AND "Australia" 
    we could automate this search using biopython Entrez"""
    pwd = "C:/Users/cow082/aivpipe"
    fasta = "%s/influenza_A_genbank.fasta" % pwd
    makeblastdb(fasta, 'influenza_A_genbank')


def test6():
    """blast multiple records against a local database
    """
    pwd = "C:/Users/cow082/aivpipe"
    fasta = "%s/21-02023-01_H5N3_complete_genome.fasta" % pwd
    records = SeqIO.parse(fasta, 'fasta')
    for i, record in enumerate(records):
        args = argParser().parse_args([
            "-record", record,
            "-program", "blastn",
            "-db", "influenza_A_genbank",
            ])
        sample, segment_number, segment_name, segment_subtype = processName(record.name)          
        blastresults = runBlast(args)
        blastresults.to_csv("%s/blastresults.%s.csv" % (pwd, segment_name))

           
def test7():
    """use os.path.walk to locate files in a folder
    """
    input_dir = "C:/Users/cow082/aivpipe/irma_assembly/datasets/21-02023-0001/FLU-avian-acdp/irma_output"
    output_dir = "C:/Users/cow082/aivpipe/blastresults"
    fasta_files = [os.path.join(root, f) for (root, dirs, files) in
            os.walk(os.path.abspath(input_dir)) 
            for f in files if f.endswith(".fasta")]
    for fasta in fasta_files:        
        blastresults = runBlast(argParser().parse_args([
            "-fasta", fasta,
            "-program", "blastn",
            "-db", "influenza_A_genbank",
            ]))
        blastresults.to_csv("%s/blastresults.%s.csv" % (
            output_dir, os.path.basename(fasta).split(".")[0]))


def argParser():
    """
    """
    parser = argparse.ArgumentParser(description = "Perform a BLAST search")
    parser.add_argument('-fasta', required=True, help = 'fasta file containing a query sequence')
    parser.add_argument('-program', default='blastp', help = 'which blast program to run?')
    parser.add_argument('-db', required=True, help = 'which database to search?')
    return parser

            
def main():
#    test1()
#    test2()
#    test3()
#    test4()
#    test5()
#    test6()
    test7()


if __name__ == "__main__":
    main()
