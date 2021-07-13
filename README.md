# AIV Pipeline Annotation Module
---

Module to take assembled influenza genomes and annotate important features such as cleavage site, pathotyping and subtyping.

Input:
* Fasta file of assembled genome from IRMA assembly module.

Output:
* BLAST analysis of each segment against the NAIWB database (plus NCBI?).
* translated amino acid sequences
* Cleavage site location in the HA segment (e.g., 758 - 823).
* Is this a novel cleavage sequence or have we seen it previously?
* Influenza pathotype (e.g., high path or low path).
* Influenza subtype (e.g., H7N7)

