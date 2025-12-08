#!/bin/env python
# conda install python3 biopython requests
# ./gbff.bnumber.location.py > ncbi.bnumber.location.tsv

#==========================================#
# SCRIPT1-ter
#==========================================#
# Script pour générer le mapping bnumber <-> location depuis un fichier gbff NCBI
# scripts/gbff.bnumber.location.py


import sys
import re
from Bio import SeqIO, Entrez
import gzip

gbff = sys.argv[1] if len(sys.argv)>1 else 'downloads/GCF_000005845.2_ASM584v2_genomic.gbff.gz'

with gzip.open(gbff, "rt") as handle:
    record = SeqIO.read(handle, "genbank")
    #print(record)
    print('bnumber\trank\tstrand\tbegin\tend')
    rank=0
    for f in record.features:
        if f.type=='CDS':
            rank += 1
            #print(f)
            f_loc = str(f.location)
            m = re.search(r"\[(\d+):(\d+)\]\((.)\)", f_loc)
            f_xref = f.qualifiers['db_xref']
            f_id = f.qualifiers['locus_tag']
            f_desc = f.qualifiers['product']
            for i in f_id:
                print(f'{i}\t{rank}\t{m[3]}\t{m[1]}\t{m[2]}')
                
#==========================================#
# END SCRIPT1-ter
#==========================================#