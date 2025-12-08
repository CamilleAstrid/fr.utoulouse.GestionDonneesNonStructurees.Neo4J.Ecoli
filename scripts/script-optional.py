#!/bin/env python

#==========================================#
# SCRIPT7
#==========================================#
# Script optionel
# Les données peuvent être récupérées depuis les résultats sur gitlab
# https://gitlab.com/rbarriot/data/-/raw/main/M2.GDNS-AP/bnumber.PMID.tsv

import sys
import datetime
import time
import requests
import os.path
import json

sleep_for = 1
n=0
last_sleep = time.time()
GeneIDs = set()
mapping = {}
with open('mapping.bnumber_ncbi.tsv', 'r') as f:
    for line in f.readlines():
        (bnumber, dbname, dbid) = line.strip().split('\t')
        if dbname == 'GeneID':
            GeneIDs.add(dbid)
            mapping[dbid] = bnumber
            resfile = 'ncbi/pubmed/GeneID.'+dbid+'.PMIDlinks.json'
            if os.path.isfile(resfile):
                sys.stderr.write('%s file %s exists for %s skipped\n' % (datetime.datetime.now(), resfile, dbid))
            else:
                url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?retmode=json&dbfrom=gene&db=pubmed&id='+dbid
                sys.stderr.write(f'fetching pmid links of {dbid} from {url}\n')
                r = requests.get(url)
                fh = open(resfile, mode='w')
                fh.write(r.text)
                fh.close()
                # timing
                n += 1
                elapsed = (time.time() - last_sleep)
                rps = n / elapsed
                sys.stderr.write("  %3d requests in %2.2f seconds → %2.1f requests per second" % (n, elapsed, rps))
                if rps > 3 : # max 3 requests/second at NCBI = 9 req in 3 sec
                    sys.stderr.write("sleeping for %s seconds\n" % (sleep_for))
                    time.sleep(sleep_for)
                    last_sleep = time.time()
                    n=0
                elif n>=10: # reset counters after 10 requests
                    last_sleep = time.time()
                    n=0

print('bnumber\tPMID')
for dbid in GeneIDs:
    sys.stderr.write("extracing PMIDs for %s\n" % (dbid))
    # extract PMIDs from XML file
    resfile = 'ncbi/pubmed/GeneID.'+dbid+'.PMIDlinks.json'
    f = open(resfile, 'r')
    d = json.load(f)
    for ls in d['linksets']:    
        for lsd in ls['linksetdbs']:
            if lsd['linkname'] == "gene_pubmed":
                for i in lsd['links']:
                    print(f'{mapping[dbid]}\t{i}')

#==========================================#
# END SCRIPT7
#==========================================#