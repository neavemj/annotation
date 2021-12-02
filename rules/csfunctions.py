# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 14:12:32 2021

@author: cow082
"""

# standard library modules
import re
from collections import Counter

# 3rd party modules
from Bio import SeqIO
import numpy as np
import pandas as pd
import matplotlib_venn as venn
from matplotlib import pyplot as plt

# custom modules
import findorf
from blastfunctions import processName

__version__ = 1.0


def findCleavageSite(record):
    """using regex, return the probable cleavage site sequence and its coordinates
    """
    pep = findorf.findORF(record).translate()
    pattern = re.compile('N[V,I][H,R,P,L,Q].{1,15}G[L,I]F')
    match = pattern.search(str(pep))
    start, end = match.span()
    if match:
        return match.group()[2:], start+2, end, len(pep)


def findSubType(s):
    """using regex, identify the subtype, e.g. H5N1
    """
    pattern = re.compile('H[0-9]{1,2}N[0-9]{1,2}')
    match = pattern.search(s)
    if match:
        return match.group()
    else:
        return ""
    

def isIn(s, tags):
    """check to see if any of a list of tags is a substring of s
    """
    for t in tags:
        if t.lower() in s.lower():
            return True
    return False


def pprintCounter(c):
    """pretty print a Counter (dict) object
    """
    d = dict(c)
    for k in sorted(d.keys()):
        print("%s\t%s" % (str(k), str(d[k])))
        

def test1():
    """basic preliminary effort to identify cleavage sites
    """
    pwd = "C:/Users/cow082/aivpipe"
#    fasta = '%s/21-02023-03_H7N8_partial_genome.fasta' % pwd
    fasta = '%s/21-02023-01_H5N3_complete_genome.fasta' % pwd
    
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
#            print(str(pep).replace(cleavage_site, cleavage_site.lower()))
        elif segment_subtype == 'H7':
            pe = pep.find('PE')
            glf = pep.find('GLF')
            cleavage_site = pep[pe:glf+3]
            cleavage_site = str(cleavage_site)
#            print(str(pep).replace(cleavage_site, cleavage_site.lower()))
        
    try:
        subtype = seg_index['HA']+seg_index['NA']
    except IndexError:
        subtype = None
        
    print(subtype)
    # H5N3
    
    print(cleavage_site)
    # PQKATRGLF
   

def test2():
    """read a fasta file and identify the cleavage sites using regex
    """
    
    # input
    pwd = "C:/Users/cow082/aivpipe"
    fasta = '%s/genbank_influenza_A_11.10.2021.fasta' % pwd
    
    # output
    sites = '%s/genbank_cleavage_sites_ALL_11.10.2021.csv' % pwd
    
    fields = ['accession', 'cleavage_site', 'start', 'end', 'segment_length', 'description']
    records = SeqIO.parse(fasta, 'fasta')
    yes, no = 0, 0
    
    with open(sites, 'w') as F:
        F.write(','.join(fields)+'\n')
        for i, record in enumerate(records):
#            if i > 100: return
            try:
                cleavage_site, start, end, length = findCleavageSite(record)
#                print(record.name, cleavage_site, start, end, length)
                F.write(','.join([record.name, cleavage_site, str(start), str(end), str(length), record.description])+'\n')
                yes += 1
            except AttributeError:
                no += 1
                
#    print(yes, no)


def test3():
    """same as test2 (above) but checks to see whether HA in label
    """
    # input
    pwd = "C:/Users/cow082/aivpipe"
    starts = [i for i in range(0, 800000, 10000)]
    fasta_files = ["%s/genbank_influenza_A_11.10.2021_%d_%d.fasta" % (pwd, start+1, start+10000) for start in starts]
    
    # output
    sites = '%s/genbank_cleavage_sites_ALL_15.10.2021.csv' % pwd

    tags = ['hemagglutinin', 'haemagglutinin', 'segment 4', '(HA)', 'HA gene', 'H gene']
    fields = ['accession', 'cleavage_site', 'start', 'end', 'segment_length', 'description']
    yes, no, notHA = 0, 0, 0
    
    with open(sites, 'w') as F:
        F.write(','.join(fields)+'\n')
        for fasta in fasta_files:
            records = SeqIO.parse(fasta, 'fasta')
            for i, record in enumerate(records):
                if isIn(record.description, tags):
                    try:
                        cleavage_site, start, end, length = findCleavageSite(record)
                        F.write(','.join([record.name, cleavage_site, str(start), str(end), str(length), record.description])+'\n')
                        yes += 1
                    except Exception:
                        no += 1
                else:
                    notHA += 1
            
    print('\n'.join([
        "yes: %d" % yes, 
        "no: %d" % no, 
        "notHA: %d" % notHA,
        ]))


def test4():
    """review the original set of cleavage sites file and try to determine which of the 
    sequences are actually HA, vs which are not, and any indeterminate or 
    unusual ones
    """
    # input
    pwd = "C:/Users/cow082/aivpipe"
    sites = '%s/genbank_cleavage_sites_ALL_11.10.2021.csv' % pwd
    
    # output
    unusual_sites = '%s/genbank_cleavage_sites_unusual_15102021.csv' % pwd
    
    tags = ['hemagglutinin', 'haemagglutinin', 'segment 4', '(HA)', 'HA gene', 'H gene']
    
    df = pd.read_csv(sites, dtype=str)
    
    df['full_description'] = df.apply(lambda row: "%s %s %s %s" % (
        str(row['description']), str(row['x']), str(row['y']), str(row['z'])), axis=1)
        
    HA = df[df.apply(lambda row: isIn(row['full_description'], tags), axis=1)]
#    print(HA)
#    [88472 rows x 9 columns]
    
    NOT_HA = df[df.apply(lambda row: not(isIn(row['full_description'], tags)), axis=1)]
#    print(NOT_HA)
#    [2398 rows x 9 columns]
    
    x = set(NOT_HA['cleavage_site']).difference(set(HA['cleavage_site']))
#    print(len(x))
#    17
    
    unusual = df[df.apply(lambda row: isIn(row['cleavage_site'], list(x)), axis=1)]
#    print(unusual[['start', 'end', 'segment_length', 'cleavage_site']])
#    [34 rows x 9 columns]

#    print('\n'.join(set(list(unusual['full_description']))))
    """ None of thse are haemagglutinin. 
    They are all other genes or weird patent sequences"""
    
    unusual.to_csv(unusual_sites)

    
def test5():
    """set analysis of cleavage sites, plus counting incidences
    """
    # input
    pwd = "C:/Users/cow082/aivpipe"
    sites1 = '%s/genbank_cleavage_sites_ALL_11.10.2021.csv' % pwd
    sites2 = '%s/genbank_cleavage_sites_ALL_15.10.2021.csv' % pwd
    sites3 = '%s/oie_cleavage_sites.csv' % pwd
    
    # count GBK cleavage sites (original regex)
    df1A = pd.read_csv(sites1)
    df2A = pd.DataFrame(zip(Counter(df1A['cleavage_site']).items()))
    df3A = df2A.apply(lambda row: pd.Series(row[0]), axis=1)
    df3A.columns = ['cleavage_site', 'count']
    df3A.sort_values('count', ascending=False, inplace=True)
    df3A.reset_index(inplace=True, drop=True)
#    print(df3A)
#    [655 rows x 2 columns]
#    print(sum(df3A['count']))
#    90870    
    
#     count GBK cleavage sites (new regex plus description filter)
    df1B = pd.read_csv(sites2)
    df2B = pd.DataFrame(zip(Counter(df1B['cleavage_site']).items()))
    df3B = df2B.apply(lambda row: pd.Series(row[0]), axis=1)
    df3B.columns = ['cleavage_site', 'count']
    df3B.sort_values('count', ascending=False, inplace=True)
    df3B.reset_index(inplace=True, drop=True)
#    print(df3B)
#    [566 rows x 2 columns]
#    print(sum(df3B['count']))
#    120208
    
    # count OIE cleavage sites
    df4 = pd.read_csv(sites3)
    df4['Cleavage site consensus'] = df4['Cleavage site consensus'].apply(lambda s: s.replace('/', ''))
    counts = Counter(df4['Cleavage site consensus'])
    df5 = pd.DataFrame(zip(counts.items()))
    df6 = df5.apply(lambda row: pd.Series(row[0]), axis=1)
    df6.columns = ['cleavage_site', 'count']
    df6.sort_values('count', ascending=False, inplace=True)
    df6.reset_index(inplace=True, drop=True)
#    print(df6)
#    [101 rows x 2 columns]
    
    # count combined GBK + OIE cleavage sites
#    df7 = pd.concat([df3A, df3B, df6], ignore_index=True)
    df7 = pd.concat([df3B, df6], ignore_index=True)
    df7.reset_index(inplace=True, drop=True)
    counts = Counter(df7['cleavage_site'])
    df8 = pd.DataFrame(zip(counts.items()))
    df9 = df8.apply(lambda row: pd.Series(row[0]), axis=1)
    df9.columns = ['cleavage_site', 'count']
    df9.sort_values('count', ascending=False, inplace=True)
    df9.reset_index(inplace=True, drop=True)
#    print(df9)
#    [809 rows x 2 columns]
    
    sites = list(df9['cleavage_site'])
    overlap = {}
    for x, X in enumerate(sites):
        for y, Y in enumerate(sites[x+1:]):
            if X in Y:
#                print(y, x, Y, X)
                overlap[Y] = X
            elif Y in X:
#                print(x, y, X, Y)
                overlap[X] = Y
            # 82 26 PRNVPQIESRGLF PQIESRGLF
    
    df9['cleavage_site'] = df9['cleavage_site'].apply(lambda s: overlap.get(s, s))

    sites = list(set(set(df9['cleavage_site'])))
    
    # verify that longer versions of sites have all been replaced by shorter ones
    for x, X in enumerate(sites):
        for y, Y in enumerate(sites[x+1:]):
            if X in Y:
                print(y, x, Y, X)
                # None
            elif Y in X:
                print(x, y, X, Y)
                # None
    
    df3A['cleavage_site'] = df3A['cleavage_site'].apply(lambda s: overlap.get(s, s))
    df3B['cleavage_site'] = df3B['cleavage_site'].apply(lambda s: overlap.get(s, s))
    df6['cleavage_site'] = df6['cleavage_site'].apply(lambda s: overlap.get(s, s))
    
    gbk = set(df3A['cleavage_site']).union(set(df3B['cleavage_site']))
    oie = set(df6['cleavage_site'])
    
#    print(len(gbk.union(oie)))
#    print(len(gbk.intersection(oie)))
#    print(len(gbk.difference(oie)))
#    print(len(oie.difference(gbk)))
#    print('\n'.join(sorted(list(gbk.intersection(oie)))))

    venn.venn2_unweighted([gbk, oie], ["gbk", "oie"])
    plt.show()
    venn.venn2([gbk, oie], ["gbk", "oie"])
    plt.show()
    venn.venn2_unweighted([set(df3A['cleavage_site']), set(df3B['cleavage_site'])], ["gbk1", "gbk2"])
    plt.show()
    
#    print('\n'.join(sorted(list(gbk.union(oie)))))
    
    
def test6():
    """examine all influenza records from genbank (~800,000), and get 
    descriptions (fasta headers)
    """
    # input
    pwd = "C:/Users/cow082/aivpipe"
    starts = [i for i in range(0, 800000, 10000)]
    fasta_files = ["%s/genbank_influenza_A_11.10.2021_%d_%d.fasta" % \
                   (pwd, start+1, start+10000) for start in starts]
    
    # output
    headers = "%s/genbank_influenza_A_11.10.2021_headers.txt" % pwd

    descriptions = []
    for fasta in fasta_files:
        for line in open(fasta, 'r').readlines():
            if line[0] == ">":
                descriptions.append(line.strip())
    with open(headers, 'w') as F:
        F.write("\n".join(descriptions))
        
    print(len(descriptions))
    # 792503
    

def test7():
    """examine all fasta headers, identify HA seqs, compare to the list of 
    accessions I was able to identify using regex to find possible 
    non-canonical cleavage sites
    """
    # input
    pwd = "C:/Users/cow082/aivpipe"
    sites = '%s/genbank_cleavage_sites_ALL_11.10.2021.csv' % pwd
    headers = "%s/genbank_influenza_A_11.10.2021_headers.txt" % pwd
    
    # output
    interesting_records = "%s/interesting.csv" % pwd

    tags = ['hemagglutinin', 'haemagglutinin', 'segment 4', '(HA)', 'HA gene', 'H gene']
    
    df = pd.read_csv(sites)
    accessions = df['accession']
#    print(accessions)
    # Name: accession, Length: 90870, dtype: object
    
    df = pd.read_csv(headers, sep='$', header=None) # $ does not exist in this file, so entire line is placed in column 1
#    print(df)
#    [792503 rows x 1 columns]
    
    HA = df[df.apply(lambda row: isIn(row[0], tags), axis=1)]
    HA.index = HA[0].apply(lambda s: s.split()[0].lstrip(">"))
#    print(HA)
    # [158144 rows x 1 columns]
    
    HA.drop(index=accessions, inplace=True, errors='ignore')
#    print(HA)
    # [68938 rows x 1 columns]
    
    HA.to_csv(interesting_records)


def test8():
    """print some of the interesting records to see if I can find why they were not picked up?
    """
    # inputs
    pwd = "C:/Users/cow082/aivpipe"
    fasta = "%s/genbank_influenza_A_11.10.2021_1_10000.fasta" % pwd
#    fasta = "%s/genbank_influenza_A_11.10.2021_10001_20000.fasta" % pwd
    
    accessions = []
    for line in open(fasta).readlines():
        if line[0] == ">":
            accessions.append(line.strip(">").split()[0])
    
    thisfile = set(accessions)
    
    df = pd.read_csv("%s/interesting.csv" % pwd, index_col=0)
    unknown = set(df["0.1"].index)
    
    overlap = thisfile.intersection(unknown)
#    print(len(overlap))
    # 447
    
    sites = []
    records = SeqIO.parse(fasta, 'fasta')
    for i, record in enumerate(records):
        if record.id in overlap:
            try:
                cleavage_site, start, end, length = findCleavageSite(record)
                if cleavage_site:
#                    print(cleavage_site)
                    sites.append(cleavage_site)
            except AttributeError: # 'NoneType' object has no attribute 'span'
#                pass
                print(record.description, findorf.findORF(record).translate())
#                pep = findorf.findORF(record).translate()
#                print(record.name, 'AttributeError', pep[300:380])
            except TypeError: # translate() takes exactly one argument (0 given)
                pass
#                print(record.name, 'TypeError', record[1:].translate()[:]) 
                                                       
    pprintCounter(Counter(sites))


def test9():
    """analyse headers and construct a csv file with info
    """
    # input
    pwd = "C:/Users/cow082/aivpipe"
    headers = "%s/genbank_influenza_A_11.10.2021_headers.txt" % pwd
    
    # output
    info = "%s/genbank_influenza_info.csv" % pwd
    
    tags = ['hemagglutinin', 'haemagglutinin', 'segment 4', '(HA)', 'HA gene', 'H gene']
    
    df = pd.read_csv(headers, sep='$', header=None)
    df.index = df[0].apply(lambda s: s.split()[0].lstrip(">"))
    
    df["subtype"] = df.apply(lambda row: findSubType(row[0]), axis=1)
    df["is_HA"] = df.apply(lambda row: isIn(row[0], tags), axis=1)
    df["HA_type"] = df.apply(lambda row: row["subtype"].split("N")[0], axis=1)
    
    df.to_csv(info)


def test10():
    """count the number of records of each type in the info file
    """
    # input
    pwd = "C:/Users/cow082/aivpipe"
    info = "%s/genbank_influenza_info.csv" % pwd
    
    df = pd.read_csv(info)
    
    subtypes = Counter(df[df["is_HA"]]['subtype'])
    print(subtypes)
    
    ha_types = Counter(df[df["is_HA"]]['HA_type'])
    print(ha_types)
    
    summary = df.groupby("HA_type").apply(lambda df: Counter(df["is_HA"]))
    print(summary)
    

def test11():
    """set analysis of cleavage sites, plus counting incidences (extention of test5)
    """
    # input
    pwd = "C:/Users/cow082/aivpipe"
    gbk_sites = '%s/genbank_cleavage_sites_ALL_15.10.2021.csv' % pwd
    oie_sites = '%s/oie_cleavage_sites.csv' % pwd
    
    # output
    extended_info = '%s/extended_info.csv' % pwd
    
#     count GBK cleavage sites (new regex plus description filter)
    df1 = pd.read_csv(gbk_sites)
    df2 = pd.DataFrame(zip(Counter(df1['cleavage_site']).items()))
    df3 = df2.apply(lambda row: pd.Series(row[0]), axis=1)
    df3.columns = ['cleavage_site', 'count']
    df3.sort_values('count', ascending=False, inplace=True)
    df3.reset_index(inplace=True, drop=True)
#    print(df3)
#    [566 rows x 2 columns]
#    print(sum(df3['count']))
#    120208
    
    # count OIE cleavage sites
    df4 = pd.read_csv(oie_sites)
    df4['cleavage_site'] = df4['Cleavage site consensus'].apply(lambda s: s.replace('/', ''))
    counts = Counter(df4['cleavage_site'])
    df5 = pd.DataFrame(zip(counts.items()))
    df6 = df5.apply(lambda row: pd.Series(row[0]), axis=1)
    df6.columns = ['cleavage_site', 'count']
    df6.sort_values('count', ascending=False, inplace=True)
    df6.reset_index(inplace=True, drop=True)
#    print(df6)
#    [101 rows x 2 columns]
    
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
#    [613 rows x 2 columns]
    
    # replace longer versions with the shortest version
    sites = list(df9['cleavage_site'])
    overlap = {}
    for x, X in enumerate(sites):
        for y, Y in enumerate(sites[x+1:]):
            if X in Y:
                overlap[Y] = X
            elif Y in X:
                overlap[X] = Y
            
    df1['cleavage_site'] = df1['cleavage_site'].apply(lambda s: overlap.get(s, s))
    df3['cleavage_site'] = df3['cleavage_site'].apply(lambda s: overlap.get(s, s))
    df4['cleavage_site'] = df4['cleavage_site'].apply(lambda s: overlap.get(s, s))
    df6['cleavage_site'] = df6['cleavage_site'].apply(lambda s: overlap.get(s, s))
    df9['cleavage_site'] = df9['cleavage_site'].apply(lambda s: overlap.get(s, s))
    
    # verify that longer versions have all been replaced (should print nothing)
    sites = list(set(set(df9['cleavage_site'])))
    count = 0
    for x, X in enumerate(sites):
        for y, Y in enumerate(sites[x+1:]):
            if X in Y:
                count += 1
            elif Y in X:
                count += 1
    
    assert count == 0
    
    gbk = set(df3['cleavage_site'])
    oie = set(df6['cleavage_site'])
#    all_sites = gbk.union(oie)
#    print('\n'.join(sorted(list(all_sites))))

    info = "%s/genbank_influenza_info.csv" % pwd
    df = pd.read_csv(info)
    HA = df[df["is_HA"]].reset_index(drop=True)
    HA.columns = ["accession", "description", "subtype", "is_HA", "HA_type"]
    df1 = df1[["accession", "cleavage_site", "segment_length", "start", "end"]]
    df10 = HA.merge(df1, how="outer", on="accession")
#    print(df10.T)
    
    def crossCheck(s):
        if s:
            return isIn(str(s), oie)
        else:
            return np.nan
    
    df10['oie_match'] = df10.apply(lambda row: crossCheck(row["cleavage_site"]), axis=1)
    df10.to_csv(extended_info)
    
    print(df10[df10["oie_match"]].reset)


def test12():
    """analyse extended info file and count the number of different cleavge
    sites for each HA sybtype
    """
    # input
    pwd = "C:/Users/cow082/aivpipe"
    extended_info = '%s/extended_info.csv' % pwd
    
    df = pd.read_csv(extended_info)
    
    print(df)
    
    count = df.groupby("HA_type").apply(lambda frame: len(set(frame["cleavage_site"])))
    
    print(count)

def main():
    """
    """
#    test1() # ok
#    test2() # write files
#    test3() # write files
#    test4() # ok
#    test5() # ok
#    test6() # ok
#    test7() # write files
#    test8() # ok
#    test9() # write files
#    test10() # ok
#    test11() # ok
#    test12() # ok
    

if __name__ == "__main__":
    main()
