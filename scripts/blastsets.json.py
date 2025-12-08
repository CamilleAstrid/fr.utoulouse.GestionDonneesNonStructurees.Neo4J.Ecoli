#!/usr/bin/env python

#==========================================#
# SCRIPT3
#==========================================#
# Script d'enrichissement utilisant les fichiers de référence JSON
# scripts/blastsets.json.py


#================================================#
# Chargement des ensembles
#================================================#

# IMPORTS
# chargement des librairies
import argparse
from os.path import isfile
import json
from scipy.stats import binom, hypergeom

# SCRIPT PARAMETERS
# prend en paramètres : la liste ou le nom de fichier contenant les identifiants composant l’ensemble requête query et le nom de fichier contenant les ensembles cibles target
parser = argparse.ArgumentParser(description='Search enriched terms/categories in the provided (gene) set')
parser.add_argument('-q', '--query', required=True, help='Query set.')
parser.add_argument('-t', '--sets', required=True, help='Target sets filename.')
parser.add_argument('-a', '--alpha', required=False, type=float, default=0.05, help='Significance threshold.')
parser.add_argument('-c', '--adjust', required=False, action="store_true", help='Adjust for multiple testing (FDR).')
parser.add_argument('-m', '--measure', required=False, default='binomial', help='Dissimilarity index: binomial (default), hypergeometric, chi2 or coverage. chi2 and coverage are NOT YET IMPLEMENTED')
parser.add_argument('-l', '--limit', required=False, type=int, default=0, help='Maximum number of results to report.')
parser.add_argument('-v', '--verbose', required=False, action="store_true", help='Talk a lot.')
param = parser.parse_args()

# LOAD QUERY
# charger les identifiants de l’ensemble requête
text = param.query
query = set()
if isfile(text):
    with open(text) as f:
        content = ' '.join(f.read().split('\n')).split()
        query |= set(content)
else: # parse string
    query |= set(text.split())

if param.verbose:
  print(f'query set: {query}')
  
# LOAD REFERENCE SETS
# ensembles cibles chargés avec la librairie jsonlite
sets = json.loads(open(param.sets).read())
if param.verbose:
    print('first target sets: ', sets[0:2])


#================================================#
# Tests statistiques
#================================================#

# COMPUTE POPULATION SIZE
# calcul de la taille de la population
population = set()
for s in sets:  
    elements = set(s['elements' ])
    population |= elements

# EVALUATE SETS
# comparaison des ensembles cibles à l’ensemble requête
results = []
query_size = len(query)
for s in sets:
    elements = set(s['elements' ])
    common_elements = elements.intersection( query )
    if param.measure=='binomial': # par défaut binom.cdf(<=success, attempts, proba). il nous faut calculer p(X>=x)
        # taille de l'ensemble cible
        target_size = len(elements)
        # taille de la population totale
        population_size = len(population)
        # nombre de succès observés
        k = len(common_elements)
        if population_size == 0:
            # pas de population définie
            pvalue = 1.0
        else:
            # Probabilité de succès attendue : T / G
            p_success = target_size / float(population_size)
            # P(X >= k) pour X ~ Binom(n = |Q|, p = T/G)
            pvalue = binom.sf(k - 1, query_size, p_success)
    else:
        print(f'sorry, {param.measure} not (yet) implemented')
        exit(1)
    r = { 'id': s['id'], 'desc': s['desc'], 'common.n':len(common_elements), 'target.n': len(elements), 'p-value': pvalue, 'elements.target': elements, 'elements.common': common_elements }
    results.append( r )
if param.verbose:
  print(results)

# PRINT SIGNIFICANT RESULTS
# afficher uniquement les résultats statistiquement significatifs
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
# END SCRIPT3
#==========================================#