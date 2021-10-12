# -*- coding: utf-8 -*-
"""
Created on Fri Jul 30 13:07:49 2021

@author: cow082

see also spread_genomics/BLASTGenotypes.py for additional ideas
"""

# standard library modules
import subprocess
import io
import re
from collections import Counter
import time

# 3rd party modules
from Bio import SeqIO, Entrez
import pandas as pd

# custom modules
import findorf

__version__ = 1.0


def getLatestGenbankRecords(target=1000000, request_limit=10000):
    """
    """
    Entrez.email = "chris.cowled@csiro.au"
    
    HOME = "C:/Users/cow082/aivpipe"
    
    # retrieve list of accession numbers
    term = '"Influenza A virus"[porgn:__txid11320] AND "Australia"'
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
        SeqIO.write(SeqIO.parse(handle, 'fasta'), '%s/genbank_influenza_A_11.10.2021_%d_%d.fasta' % (HOME, i+1, i+request_limit), "fasta")       
        handle.close()
        i += request_limit
#        time.sleep(0.35) # max 3 requests per second (not needed for batches as they take longer to process)
    
    return len(IdList)
    # as of 11.10.21 there are 18560  records for '"Influenza A virus"[porgn:__txid11320] AND "Australia"'
    # as of 11.10.21 there are 792486 records for '"Influenza A virus"[porgn:__txid11320]'
    


def makeblastdb(fasta, title, dbtype="nucl"):
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
    
    stdout = io.StringIO(p.stdout).read()
    stderr = io.StringIO(p.stderr).read()
    print(stdout, stderr)
    
    return p.returncode

    
def blast(record, program, 
          db="nr", organism="Influenza A virus [Organism]"):
    """seems quite slow to do an external query
    """
    args = [program, "-db", db]
    
    if db == "nr":
        args.extend(["-remote", "-entrez_query", organism])
        
    p = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    out, _ = p.communicate(record.encode())
    return out.decode()

    
def test1():
    """ blast nr using remote blast
    """
    HOME = "C:/Users/cow082/aivpipe"
    fasta = '%s/21-02023-03_H7N8_partial_genome.fasta' % HOME
    records = SeqIO.parse(fasta, 'fasta')
    for i, record in enumerate(records):
        if i > 0:
            break
        orf = findorf.findORF(record)
        pep = orf.translate()
        blastresults = blast(pep, program='blastp')
        print(blastresults)


def test2():
    """ build a local database
    
    fasta is result of ncbi Nucleotide search: 
            "Influenza A virus"[porgn:__txid11320] AND "Australia" 
    we could automate this search using biopython Entrez"""
    HOME = "C:/Users/cow082/aivpipe"
    fasta = '%s/influenza_A_genbank.fasta' % HOME
    makeblastdb(fasta, 'influenza_A_genbank')


def recordToFastaString(record):
    """convert a seqrecord to a fasta string
    """
    out_handle = io.StringIO()
    SeqIO.write(record, out_handle, "fasta")
    return out_handle.getvalue()


def runBlast(record, program, db):
    """perform a blast search
    """
    
#    validateBlastDb(genotype_index)
        
    # http://www.metagenomics.wiki/tools/blast/blastn-output-format-6
    fields = 'qseqid sseqid qlen slen pident length mismatch gapopen qstart qend sstart send evalue bitscore'# qseq sseq'
    
    p = subprocess.run([
        program, 
        '-db', db, 
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
    except AssertionError as e:
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
        return 'no hits'
    
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
        print("processName() encountered an unexpected format")
        return None


def test3():
    """blast against local database
    """
    HOME = "C:/Users/cow082/aivpipe"
    fasta = '%s/21-02023-03_H7N8_partial_genome.fasta' % HOME
    fasta = '%s/21-02023-01_H5N3_complete_genome.fasta' % HOME
    records = SeqIO.parse(fasta, 'fasta')
    seg_index = {}
    cleavage_site = None
    for i, record in enumerate(records):
        sample, segment_number, segment_name, segment_subtype = processName(record.name)
#        print(segment_number, segment_name, segment_subtype)
        """
        1 PB2 None
        2 PB1 None
        3 PA None
        4 HA H5
        5 NP None
        6 NA N3
        7 MP None
        8 NS None
        """
        seg_index[segment_name] = segment_subtype
        
        orf = findorf.findORF(record)
#        nuc = str(orf)
        pep = orf.translate()
        
        if segment_subtype == 'H5':
            pq = pep.find('PQ')
            glf = pep.find('GLF')
            cleavage_site = pep[pq:glf+3]
            cleavage_site = str(cleavage_site)
            print(str(pep).replace(cleavage_site, cleavage_site.lower()))
        elif segment_subtype == 'H7':
            pe = pep.find('PE')
            glf = pep.find('GLF')
            cleavage_site = pep[pe:glf+3]
            cleavage_site = str(cleavage_site)
            print(str(pep).replace(cleavage_site, cleavage_site.lower()))
            
#        blastresults = runBlast(record, program='blastn', db='influenza_A_genbank')
#        print(blastresults)
        
    try:
        subtype = seg_index['HA']+seg_index['NA']
    except IndexError:
        subtype = None
        
    print(subtype)
    # H5N3
    
    print(cleavage_site)
    # 


def test4():
    """pull down all IAV records from genbank
    """
    count, records = getLatestGenbankRecords(target=10, request_limit=3)
    for record in records:
        print(record)
    print(count)


def downloadGBK():
    count = getLatestGenbankRecords()
    print(count)
    
    
def test5():
    Entrez.email = "chris.cowled@csiro.au"
    handle = Entrez.efetch(db="nucleotide", id="gb|CY022693.1|", rettype="gb", retmode="text")
    print(handle.read())


def findCleavageSite(record):
    """
    """
    pep = findorf.findORF(record).translate()
    pattern = re.compile('P.{1,15}GLF') # P, followed by 1-15 of any residue, followed by GLF
    match = pattern.search(str(pep))
    start, end = match.span()
    if match:
        return match.group(), start, end, len(pep)


def test6():
    HOME = "C:/Users/cow082/aivpipe"
#    fasta = '%s/influenza_A_genbank.fasta' % HOME
    fasta = '%s/genbank_influenza_A_11.10.2021.fasta' % HOME
    records = SeqIO.parse(fasta, 'fasta')
    yes, no = 0, 0
    with open('%s/genbank_cleavage_sites_11.10.2021.csv' % HOME, 'w') as F:
        F.write(','.join(['accession', 'cleavage_site', 'start', 'end', 'segment_length', 'description'])+'\n')
        for i, record in enumerate(records):
#            if i > 100: return
            try:
                cleavage_site, start, end, length = findCleavageSite(record)
                print(record.name, cleavage_site, start, end, length)
                F.write(','.join([record.name, cleavage_site, str(start), str(end), str(length), record.description])+'\n')
                yes += 1
            except AttributeError:
                no += 1
    print(yes, no)


def test6B():
    HOME = "C:/Users/cow082/aivpipe"
    starts = [i for i in range(0, 800000, 10000)]
    fasta_files = ["%s/genbank_influenza_A_11.10.2021_%d_%d.fasta" % (HOME, start+1, start+10000) for start in starts]
    yes, no = 0, 0
    with open('%s/genbank_cleavage_sites_ALL_11.10.2021.csv' % HOME, 'w') as F:
        F.write(','.join(['accession', 'cleavage_site', 'start', 'end', 'segment_length', 'description'])+'\n')
        for fasta in fasta_files:
            records = SeqIO.parse(fasta, 'fasta')
            for i, record in enumerate(records):
                try:
                    cleavage_site, start, end, length = findCleavageSite(record)
#                    print(record.name, cleavage_site, start, end, length)
                    F.write(','.join([record.name, cleavage_site, str(start), str(end), str(length), record.description])+'\n')
                    yes += 1
                except Exception:
                    no += 1
            print(fasta)
    print(yes, no)


def test6C():
    HOME = "C:/Users/cow082/aivpipe"
    starts = [i for i in range(0, 800000, 10000)]
    fasta_files = ["%s/genbank_influenza_A_11.10.2021_%d_%d.fasta" % (HOME, start+1, start+10000) for start in starts]
    count = 0
    for fasta in fasta_files:
        records = open(fasta, 'r').read().count('>')
        count += open(fasta, 'r').read().count('>')
        print(fasta, records)
    print(count)
    # 160000
        
    

def test7():
    HOME = "C:/Users/cow082/aivpipe"
    
    # count GBK cleavage sites
#    sites = '%s/genbank_cleavage_sites.csv' % HOME
#    sites = '%s/genbank_cleavage_sites_11.10.2021.csv' % HOME
    sites = '%s/genbank_cleavage_sites_ALL_11.10.2021.csv' % HOME
    df = pd.read_csv(sites)
    counts = Counter(df['cleavage_site'])
    df2 = pd.DataFrame(zip(counts.items()))
    df3 = df2.apply(lambda row: pd.Series(row[0]), axis=1)
    df3.columns = ['cleavage_site', 'count']
    df3.sort_values('count', ascending=False, inplace=True)
    df3.reset_index(inplace=True, drop=True)
#    print(df3)
#    print(sum(df3['count']))
    # 871
    
    # count OIE cleavage sites
    sites = '%s/oie_cleavage_sites.csv' % HOME
    df4 = pd.read_csv(sites)
    df4['Cleavage site consensus'] = df4['Cleavage site consensus'].apply(lambda s: s.replace('/', ''))
    counts = Counter(df4['Cleavage site consensus'])
    df5 = pd.DataFrame(zip(counts.items()))
    df6 = df5.apply(lambda row: pd.Series(row[0]), axis=1)
    df6.columns = ['cleavage_site', 'count']
    df6.sort_values('count', ascending=False, inplace=True)
    df6.reset_index(inplace=True, drop=True)
#    print(df6)
    
    # count combined GBK + OIE cleavage sites
    df7 = pd.concat([df3, df6], ignore_index=True)
    df7.reset_index(inplace=True, drop=True)
    counts = Counter(df7['cleavage_site'])
    df8 = pd.DataFrame(zip(counts.items()))
    df9 = df8.apply(lambda row: pd.Series(row[0]), axis=1)
    df9.columns = ['cleavage_site', 'count']
    df9.sort_values('count', ascending=False, inplace=True)
    df9.reset_index(inplace=True, drop=True)
#    print(df9)
    
    sites = list(df9['cleavage_site'])
    for x, X in enumerate(sites):
        for y, Y in enumerate(sites[x+1:]):
            if X in Y:
                print(x, y, X, Y)
            elif Y in X:
                print(x, y, X, Y)
            # 82 26 PRNVPQIESRGLF PQIESRGLF
    
    gbk = set(df3['cleavage_site'])
    oie = set(df6['cleavage_site'])
    
    print()
    print(len(gbk.union(oie)))
    print(len(gbk.intersection(oie)))
    print(len(gbk.difference(oie)))
    print(len(oie.difference(gbk)))

    print()
    print('\n'.join(sorted(list(gbk.intersection(oie)))))


def main():
#    test1()
#    test2()
#    test3()
#    test4()
#    test5()
#    test6()
    
#    downloadGBK()
#    test6C()
    test6B()
#    test7()


if __name__ == "__main__":
    main()
