#!/usr/bin/env python

#==========================================#
# SCRIPT7
#==========================================#
# Script d'enrichissement utilisant une connexion au serveur Neo4J
# scripts/blastsets.neo4j.py

# IMPORTS
# chargement des librairies
import argparse
from os.path import isfile
from scipy.stats import binom, hypergeom
from neo4j import GraphDatabase
import sys


#------------------------------------------------#
# Paramètres de ligne de commande
#------------------------------------------------#
parser = argparse.ArgumentParser(description='Search enriched terms/categories in the provided (gene) set using Neo4j')
parser.add_argument('-q', '--query', required=True, help='Query set (filename or space-separated list of gene IDs).')
parser.add_argument('-t', '--target', dest='target_type', required=True, help='Target node label in Neo4j (e.g. Pathway, TU, PubMed, InterPro, Keyword).')
parser.add_argument('-x', '--taxid', required=False, type=int, default=511145, help='NCBI taxonomy ID of the organism (default: 511145 for E. coli).')
parser.add_argument('-a', '--alpha', required=False, type=float, default=0.05, help='Significance threshold.')
parser.add_argument('-c', '--adjust', required=False, action="store_true", help='Adjust for multiple testing (FDR, Benjamini-Hochberg).')
parser.add_argument('-m', '--measure', required=False, default='binomial', help='Dissimilarity index: binomial (default) or hypergeometric (not yet implemented here).')
parser.add_argument('-l', '--limit', required=False, type=int, default=0, help='Maximum number of results to report.')
parser.add_argument('-v', '--verbose', required=False, action="store_true", help='Talk a lot.')
parser.add_argument('-s', '--species', required=False, type=int, help='Alias for taxid (same as -x/--taxid).')

# paramètres de connexion Neo4j
parser.add_argument('--uri', required=False, default='bolt://localhost:7687', help='Neo4j URI (default: bolt://localhost:7687).')
parser.add_argument('--user', required=False, default='neo4j', help='Neo4j username (default: neo4j).')
parser.add_argument('--password', required=False, default='omybioinfo', help='Neo4j password (default: omybioinfo).')
param = parser.parse_args()

if param.species is not None:
    param.taxid = param.species

#------------------------------------------------#
# Chargement de l’ensemble requête
#------------------------------------------------#
text = param.query
query = set()
if isfile(text):
    with open(text) as f:
        content = ' '.join(f.read().split('\n')).split()
        query |= set(content)
else:  # parse string
    query |= set(text.split())

if param.verbose:
    print(f'[INFO] Query set ({len(query)} IDs): {sorted(query)}', file=sys.stderr)


#------------------------------------------------#
# Récupération des ensembles cibles depuis Neo4j
#------------------------------------------------#

# Mapping type de sommet -> type de relation vers Gene
RELATION_MAP = {
    'Pathway': 'requires', # (p:Pathway)-[:requires]->(g:Gene)
    'TU': 'harbors', # (t:TU)-[:harbors]->(g:Gene)
    'InterPro': 'harbored_by', # (i:InterPro)-[:harbored_by]->(g:Gene)
    'PubMed': 'cites', # (p:PubMed)-[:cites]->(g:Gene)
    'Keyword': 'describes' # (k:Keyword)-[:describes]->(g:Gene)
}

target_label = param.target_type
if target_label != 'GOTerm' and target_label not in RELATION_MAP:
    print(f"[ERREUR] Type de sommet cible inconnu : {target_label}. "
          f"Types supportés : Pathway, TU, PubMed, InterPro, Keyword, GOTerm",
          file=sys.stderr)
    sys.exit(1)

rel_type = RELATION_MAP.get(target_label, None)  # None pour GOTerm


def fetch_target_sets(tx, label, rel, taxid, verbose=False):
    """
    Récupère depuis Neo4j les ensembles cibles :
    - cas général : (s:Label)-[:rel]->(g:Gene)
    - cas GOTerm : (s:GOTerm)-[:annotates]->(g:Gene)<-[:encoded_by]-(:GeneProduct)
      (en n’oubliant pas les associations GOterm → GeneProduct implicites)

    Retourne une liste de dicts :
        { 'id': ..., 'desc': ..., 'elements': [gene1, gene2, ...] }
    """

    if label == 'GOTerm':
        cypher_query = """
        MATCH (s:GOTerm)-[:annotates]->(g:Gene)<-[:encoded_by]-(:GeneProduct)
        WHERE g.organism = $taxid
        RETURN s.id AS id,
               coalesce(s.desc, s.id) AS desc,
               collect(DISTINCT g.id) AS elements
        """
    else:
        cypher_query = f"""
        MATCH (s:{label})-[:{rel}]->(g:Gene)
        WHERE g.organism = $taxid
        RETURN s.id AS id,
               coalesce(s.desc, s.id) AS desc,
               collect(DISTINCT g.id) AS elements
        """

    if verbose:
        print(f"[INFO] Cypher query:\n{cypher_query}", file=sys.stderr)

    result = tx.run(cypher_query, taxid=taxid)
    sets = []
    for record in result:
        sets.append({
            'id': record['id'],
            'desc': record['desc'],
            'elements': list(record['elements'])
        })
    return sets


# Connexion à Neo4j et chargement des ensembles
driver = GraphDatabase.driver(param.uri, auth=(param.user, param.password))

with driver.session() as session:
    sets = session.execute_read(
        fetch_target_sets,
        target_label,
        rel_type,
        param.taxid,
        verbose=param.verbose
    )

driver.close()

if param.verbose:
    print(f"[INFO] Nombre d'ensembles cibles de type {target_label} : {len(sets)}",
          file=sys.stderr)
    if len(sets) > 0:
        print(f"[INFO] Exemples d'ensembles cibles : {sets[0:2]}", file=sys.stderr)


#================================================#
# Tests statistiques
#================================================#

# COMPUTE POPULATION SIZE
# calcul de la taille de la population
population = set()
for s in sets:
    elements = set(s['elements'])
    population |= elements

if param.verbose:
    print(f"[INFO] Taille de la population (nb gènes distincts dans les ensembles cibles) : "
          f"{len(population)}", file=sys.stderr)

# EVALUATE SETS
results = []
query_size = len(query)

if query_size == 0:
    print("[ERREUR] L'ensemble requête est vide.", file=sys.stderr)
    sys.exit(1)

for s in sets:
    elements = set(s['elements'])
    common_elements = elements.intersection(query)

    if param.measure == 'binomial':
        target_size = len(elements)
        population_size = len(population)
        k = len(common_elements)

        if population_size == 0:
            pvalue = 1.0
        else:
            p_success = target_size / float(population_size)
            # P(X >= k) pour X ~ Binom(n = |Q|, p = T/G)
            pvalue = binom.sf(k - 1, query_size, p_success)
    else:
        print(f"Sorry, measure '{param.measure}' not (yet) implemented.", file=sys.stderr)
        sys.exit(1)

    r = {
        'id': s['id'],
        'desc': s['desc'],
        'common.n': len(common_elements),
        'target.n': len(elements),
        'p-value': pvalue,
        'elements.target': elements,
        'elements.common': common_elements
    }
    results.append(r)

if param.verbose:
    print(f"[INFO] Nombre d'ensembles avec p-valeur calculée : {len(results)}",
          file=sys.stderr)

# PRINT SIGNIFICANT RESULTS
# tri par p-valeur croissante
results.sort(key=lambda an_item: an_item['p-value'])

# ajustement pour tests multiples (FDR)
if param.adjust:
    m = len(results)
    if m > 0:
        k_max = 0
        # résultats déjà triés par p-valeur croissante
        for i, r in enumerate(results, start=1):
            # condition BH : P_i <= (i/m) * alpha
            if r['p-value'] <= (i / float(m)) * param.alpha:
                k_max = i
        # si aucun test ne satisfait la condition, aucun résultat ne sera significatif
        if k_max == 0:
            param.alpha = 0.0
        else:
            # seuil effectif = p-valeur du dernier test retenu (rang k_max)
            param.alpha = results[k_max - 1]['p-value']
            
for r in results:
    if r['p-value'] > param.alpha: 
      break
    # OUTPUT
    print("{}\t{}\t{}\t{}/{}\t{}\t{}".format(r['id'],r['p-value'],r['common.n'],r['common.n'],r['target.n'],r['desc'],', '.join(r['elements.common'])))

#==========================================#
# END SCRIPT7
#==========================================#
