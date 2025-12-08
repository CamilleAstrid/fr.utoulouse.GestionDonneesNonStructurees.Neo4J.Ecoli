#!/bin/env python
# conda install python3 biopython requests
# ./gbff.to.bnumber.GeneID.py > mapping.bnumber_ncbi.tsv

#==========================================#
# SCRIPT1-bis
#==========================================#
# Script pour générer le mapping bnumber <-> GeneID depuis un fichier gbff NCBI
# scripts/gbff.to.bnumber.GeneID.py


import sys
import gzip
from Bio import SeqIO, Entrez

gbff = sys.argv[1] if len(sys.argv)>1 else 'downloads/GCF_000005845.2_ASM584v2_genomic.gbff.gz'

with gzip.open(gbff, "rt") as handle:
    record = SeqIO.read(handle, "genbank")
    #print(record)
    print('bnumber\tdbname\tdbid')
    for f in record.features:
        if f.type=='CDS':
            #print(f)
            f_loc = f.location
            f_xref = f.qualifiers['db_xref']
            f_id = f.qualifiers['locus_tag']
            f_desc = f.qualifiers['product']
            for i in f_id:
                for j in f_xref:
                    (dbsource, dbid) = j.split(':')
                    dbsource = 'UniProtKB' if dbsource.startswith('UniProtKB') else dbsource
                    print(f'{i}\t{dbsource}\t{dbid}')

#==========================================#
# END SCRIPT1-bis
#==========================================#