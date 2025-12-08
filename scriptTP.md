# E. coli - Gestion de données non structurées - Applications post-génomiques - Mise en oeuvre

```R
#knitr::opts_chunk$set(echo = TRUE)
options(width = 95) # nombre de colonnes pour l'affichage sur la sorties HTML
# knitr::opts_chunk$set(cache=F)
# knitr::opts_chunk$set(fig.width=6, fig.height=6) # for rstudio
# knitr::opts_chunk$set(fig.width=14, fig.height=14) # for html
```

Environnement conda/mamba
```bash
mamba create -n gdnsap r-tidyverse r-reticulate r-codetools r-dt r-kableextra bioconductor-stringdb r-rmdformats r-factoextra scipy   igraph py2neo monotonic packaging numpy pandas scikit-learn biopython   r-markdown
```

Puis installer neo2R depuis le CRAN
```R
install.packages("neo2R")
```

Librairie utilisée pour lancer quelques commandes bash ou python sous Rstudio pour générer ce RNotebook.
```R
library(reticulate)
```

Répertoires qui seront utilisés pour les données, les scripts, etc. :
```bash
mkdir downloads      # données d'origine téléchargées
mkdir reference.sets # ensembles de référence générés pour la recherche d'enrichissement
mkdir query.sets     # ensembles requête
mkdir neo4j.import   # répertoire pour les fichiers de données à importer dans neo4j
mkdir neo4j.data     # répertoire où neo4j stocke ses données
mkdir scripts        # scripts de recherche d'enrichissement et de prioritization
```

Illustration sur un organisme. Constitution d’une base dédiée à un organisme et son exploitation.

## 1 Approche par ensembles : Enrichissement

Il a été identifié un ou quelques ensembles de gènes d’intérêt chez E. coli. A quels processus biologiques, fonctions moléculaires ou localisations sub-cellulaires peut-on les relier ?

Ensembles : set 1, set 2, set 3

Utiliser l’interface Web http://silico.biotoul.fr/enrichment/ afin de caractériser ces ensembles.

backup wep app: http://enrichment.inn.ovh:8080

### 1.1 Ensembles cibles et script de recherche d’enrichissement

Un premier exercice, relativement simple, va être de reproduire ce type d’analyse. Nous allons utiliser les mots-clés associés aux protéines dans UniProt.

Les étapes vont être les suivantes :
* récupération du protéome : toutes les protéines avec les annotations qui nous intéressent, notamment les keywords pour commencer
* reformatage : liste des protéines annotées par chaque keyword
* écriture d’un script python, qui prend en paramètre une liste de protéines et qui cherchent les keywords enrichis en comparant cette liste aux listes associées à chacun des keywords

#### 1.1.1 Génération des ensembles cibles (UniProt keywords, Interpro domains, EcoCyc pathways, EcoCyc transcription units)

Aller sur UniProt et télécharger le protéome correspondant à E. coli K12 MG1655 ; étapes :

→ sélectionner Species Proteomes puis rechercher Escherichia coli K12 MG1655. Identifier le bon résultat (protéome de référence) et cliquer sur son Proteome ID   
→ Cliquer ensuite sur le lien pour télécharger toutes les protéines   
→ Dans les paramètres pour le téléchargement, spécifier de télécharger toutes les protéines,   
→ au format TSV,   
→ sélectionner les colonnes/champs que nous allons utiliser par la suite :   
* Entry name
* Gene names (ordered locus)
* Gene Ontology IDs
* Interpro
* Keywords

→ Puis télécharger le fichier (uniprotkb_proteome_UP000000625_2023_10_11tsv.gz) dans un répertoire destiné à garder les fichiers d’origine récupérés (par exemple downloads)

```bash
zcat downloads/uniprotkb_proteome_UP000000625_*.tsv.gz | head | column -ts $'\t' 
```

Pour générer un fichier texte au format simple qui donne pour chaque mot-clé, la liste des protéines associées, nous allons utiliser R et les librairies tidyverse et jsonlite.

##### 1.1.1.1 Choix des identifiants de référence

Un des problèmes souvent rencontrés dans ce type d’analyse est d’arriver à identifier : les génomes et les gènes, protéines, etc. En effet, chaque source de données peut utiliser des méthodes de référencement et d’identification qui lui est propre. Pour le protéome d’E. coli, nous allons utiliser un type d’identifiant très souvent utilisé : les bnumbers.

A partir du fichier téléchargé, vous allez donc le charger sous R et le reformater pour avoir un mapping des identifiants UniProt vers les bnumbers.

```R
library(tidyverse)

#Chargement du fichier et son contenu
uniprot <- read_tsv("downloads/uniprotkb_proteome_UP000000625_2024_09_19.tsv.gz", show_col_types = F)
uniprot %>% head
```

A partir de ce tibble, il s’agit d’extraire les bnumbers de la colonne Gene names (ordered locus). La fonction str_extract permet d’extraire un pattern à partir d’une expression régulière (ici, ce sera 'b\\d+') d’une chaîne de caractères.

```R
mapping <- uniprot %>% 
  select(Entry, names = `Gene Names (ordered locus)`) %>%
  mutate(bnumber=str_extract(names, 'b\\d+')) %>%
  select(bnumber, uniprotID=Entry) %>%
  filter(!is.na(bnumber)) %>% # 2023 → P0DQD7 and P0A6D5 are lost (no bnumber)
  arrange(bnumber)
bnumber.uniprot <- mapping
mapping %>% head
```

##### 1.1.1.2 Reformatage des données

A présent, il s’agit de générer la liste des protéines associées à chaque mot-clé (colonne Keywords). Pour cela, on ne gardera que les colonnes Entry et Keywords, une jointure avec mapping va permettre de faire le mapping avec les bnumbers. La fonction separate_rows va nous servir à découper la colonne Keywords et répartir les mots-clés sur des lignes du tibble :

```R
keywords <- uniprot %>% 
  select(uniprotID=Entry, keyword=Keywords) %>%
  right_join(mapping, by="uniprotID") %>% # right join to remove those without bnumber
  separate_rows(keyword, sep=';') %>%
  select(bnumber, keyword) %>%
  arrange(bnumber)
keywords %>% head
```

Maintenant, il ne s’agit plus que de regrouper les bnumbers par mot-clé :

```R
ref_sets <- keywords %>% 
  group_by(keyword) %>%
  summarise(count=n(), elements = list(bnumber)) %>%
  ungroup %>%
  filter(count>1)  %>%
  select(id=keyword, desc=count, elements)
ref_sets
```

Pour utiliser “facilement” ces données avec python, nous allons sauvegarder le tibble au format JSON :

```R
library(jsonlite)

ref_sets %>% 
  toJSON %>% 
  write("reference.sets/uniprot.keywords.sets.json")
```

Apperçu du contenu :
```bash
jq . "reference.sets/uniprot.keywords.sets.json" | head -40
```

#### 1.1.2 Script de recherche d’enrichissement

La recherche d’enrichissement se résume à comparer chacun des ensembles cibles (pour ce premier exemple, constitués des protéines associées à un même mot-clé) au moyen d’un test statistique. Nous utiliserons pour commencer une loi binomiale.

##### 1.1.2.1 Paramètres de la ligne de commande

Ecrire un script python qui prendra en paramètres :
* la liste ou le nom de fichier contenant les identifiants composant l’ensemble requête query
* le nom de fichier contenant les ensembles cibles target

On vous fournit le début du script :
```python
#!/usr/bin/env python

import argparse
from os.path import isfile
import json
from scipy.stats import binom, hypergeom

# SCRIPT PARAMETERS
parser = argparse.ArgumentParser(description='Search enriched terms/categories in the provided (gene) set')
parser.add_argument('-q', '--query', required=True, help='Query set.')
parser.add_argument('-t', '--sets', required=True, help='Target sets filename.')
parser.add_argument('-a', '--alpha', required=False, type=float, default=0.05, help='Significance threshold.')
parser.add_argument('-c', '--adjust', required=False, action="store_true", help='Adjust for multiple testing (FDR).')
parser.add_argument('-m', '--measure', required=False, default='binomial', help='Dissimilarity index: binomial (default), hypergeometric, chi2 or coverage. chi2 and coverage are NOT YET IMPLEMENTED')
parser.add_argument('-l', '--limit', required=False, type=int, default=0, help='Maximum number of results to report.')
parser.add_argument('-v', '--verbose', required=False, action="store_true", help='Talk a lot.')
param = parser.parse_args()
```

##### 1.1.2.2 Ensemble requête

Il s’agit ensuite de charger les identifiants de l’ensemble requête :
```python
# LOAD QUERY
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
```

##### 1.1.2.3 Ensembles cibles

Les ensembles cibles peuvent être chargés avec la librairie jsonlite (importée au début):

```python
# LOAD REFERENCE SETS
sets = json.loads(open(param.sets).read())
if param.verbose:
    print('first target sets: ', sets[0:2])
```

##### 1.1.2.4 Test statistique

Afin d’appliquer un test avec la loi binomiale, nous avons besoin de connaître la taille de la population pour calculer la probabilité de succès ou d’échec à chaque tentative (nombre d’éléments pris au hasard dans la population, avec remise).

Fonction disponible dans scipi.stats → https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.binom.html

Rappel sur la loi binomiale :   
k succès obtenus sur n tentatives, chaque tentative ayant une probabilité p de succès. Pour chaque ensemble cible T, ayant C éléments en commun avec l’ensemble requête Q, il s’agit donc de considérer la probabilité d’avoir au moins C succès en effectuant Q tirages aléatoires avec remise, ayant chacun une probabilité de T/G de succès, G étant le nombre de gènes/protéines considérés pour le tirage.

###### 1.1.2.4.1 Taille de la population

Pour appliquer le test, il nous faut déterminer la probabilité de succès et donc le nombre total d’identifiants possibles. Ajouter une partie dans le script qui calcule le nombre d’identifiants “disponibles” parmi les ensembles cibles :

```python
# COMPUTE POPULATION SIZE
population = set()
for s in sets:
  # TO DO
```

###### 1.1.2.4.2 Comparaison des ensembles cibles à l’ensemble requête

Ensuite, il faut compléter la partie TO DO ci-dessous pour stocker la p-valeur obtenue pour chaque test.

```python
# EVALUATE SETS
results = []
query_size = len(query)
for s in sets:
    elements = set(s['elements' ])
    common_elements = elements.intersection( query )
    if param.measure=='binomial': # par défaut binom.cdf(<=success, attempts, proba). il nous faut calculer p(X>=x)
        pvalue = 100000 # TO DO 
    else:
        print(f'sorry, {param.measure} not (yet) implemented')
        exit(1)
    r = { 'id': s['id'], 'desc': s['desc'], 'common.n':len(common_elements), 'target.n': len(elements), 'p-value': pvalue, 'elements.target': elements, 'elements.common': common_elements }
    results.append( r )
if param.verbose:
  print(results)
```

##### 1.1.2.5 Affichage des résultats significatifs

Une fois l’ensemble des ensembles cibles comparés et évalués, il s’agit d’afficher uniquement les résultats statistiquement significatifs :

```python
# PRINT SIGNIFICANT RESULTS
results.sort(key=lambda an_item: an_item['p-value'])
for r in results:
    if r['p-value'] > param.alpha: 
      break
    # OUTPUT
    print("{}\t{}\t{}/{}\t{}\t{}".format( r['id'], r['p-value'], r['common.n'], r['target.n'], r['desc'], ', '.join(r['elements.common'])))
```

Exemple d’utilisation avec set.01.txt :
```bash
./scripts/blastsets.json.py -q query.sets/set.01.txt -t reference.sets/uniprot.keywords.sets.json  | column -ts $'\t'
```

Nombre de résultats significatifs obtenus :
```bash
./scripts/blastsets.json.py -q query.sets/set.01.txt -t reference.sets/uniprot.keywords.sets.json | wc -l
```

###### 1.1.2.5.1 Ajustement dans le cadre de tests multiples → FDR

Il s’agit maintenant d’ajouter une correction du seuil α dans le cadre de ctests multiples. La FDR consiste à ajuster ce seuil par rapport aux p-valeurs obtenues : Les p-valeurs Pi sont triées de manière croissante (P1 à Pm, m étant le nombre de tests) et on garde les k plus petites p-valeurs vérifiant Pk≤kmα

Nombre de résultats significatifs obtenus :
```bash
./scripts/blastsets.json.py --adjust -q query.sets/set.01.txt -t reference.sets/uniprot.keywords.sets.json | wc -l
```

A partir de toutes ces informations et du code fourni, écrire le script python permettant de faire la recherche d’ensembles cibles similaires.

#### 1.1.3 Autres ensembles d’ensembles cibles

Même chose à faire pour

* les domaines Interpro à partir du protéome téléchargé UniProt Proteome (colonne InterPro). Pour la description des domaines associée leur identifiant, vous pourrez utiliser le format court (InterPro entry list) proposé sur la page téléchargement du site InterPro https://www.ebi.ac.uk/interpro/download/
* les GOTerms référencés dans le même fichier d’UniProt (colonne Gene ontology IDs)
* les voies métaboliques (pathways) et les unités de transcription (TUs) que vous trouverez sur EcoCyc (Faire un Export → to Spreadsheet File... format frame IDs) :
    * gènes (pour le mapping) : https://www.ecocyc.org/group?id=biocyc13-55140-3842501533
    * pathways : https://www.ecocyc.org/group?id=biocyc17-55140-3842483872
    * unités de transcription : https://www.ecocyc.org/group?id=biocyc17-55140-3842483291

Adaptez les traitements effectués pour générer les ensembles cibles de référence pour les mots-clés UniProt, afin de générer les ensembles cibles de référence pour les domaines InterPro, GOTerms, les pathways et les TUs.

Une autre source qui peut s’avérer intéressante est la littérature biomédicale relative à chacun·e des gènes/protéines. Pour cela Entrez du NCBI va nous permettre de constituer, pour chaque article référençant au moins 2 gènes d’E. coli, un ensemble cible dans la section suivante.

##### 1.1.3.1 Ensembles de gènes cités dans une même publication (PubMed)

###### 1.1.3.1.1 Ficher d’annotation du génome

A partir du site du NCBI, avec Entrez, rechercher E. coli K-12 MG1655 : https://www.ncbi.nlm.nih.gov/genome/?term=escherichia+coli+k-12+MG1655

Dans les NCBI Datasets obtenus, le problème, encore une fois est d’identifier la bonne version des données. Pour cet exemple, nous utiliserons le génome de référence et son annotation RefSeq ASM584v2 → RefSeq:GCF_000005845.2

Télécharger le ficher au format GBFF. Nous aurons besoin ensuite du module biopython pour analyser le contenu et extraire les gènes codants (pour des protéines).
```bash
unzip -l downloads/GCF_000005845.2.zip 
```

Extract gbff to compressed file
```bash
unzip -p downloads/GCF_000005845.2.zip ncbi_dataset/data/GCF_000005845.2/genomic.gbff | gzip > downloads/GCF_000005845.2.gbff.gz
```

Nous allons ici ne conserver que les GeneID (pour chercher les références bibliographiques) et la localisation des gènes sur le chromosome :

Attention : pour ne pas surcharger le serveur du NCBI (et se faire bloquer), il faut passer par leur API et limiter le nombre de requêtes effectuées par seconde.

###### 1.1.3.1.2 Mapping des identifiants

Encore une fois, il faut choisir les “bons” identifiants parmi ceux utiliser. Ici, nous utiliserons les bnumbers.

Pour extraire le mapping à partir du fichier GenBank, créer le script suivant (gbff.to.bnumber.GeneID.py) dans le répertoire data et l’exécuter :

```python
#!/bin/env python
# conda install python3 biopython requests
# ./gbff.to.bnumber.GeneID.py > mapping.bnumber_ncbi.tsv

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
```

On obtient à la sortie le mapping entre identifiants :
```bash
./scripts/gbff.to.bnumber.GeneID.py downloads/GCF_000005845.2.gbff.gz > generated.data/mapping.bnumber_ncbi.tsv
head generated.data/mapping.bnumber_ncbi.tsv | column -ts $'\t'
```

###### 1.1.3.1.3 Génération des ensembles cibles (PubMed ids)

Pour chaque gène/identifiant/GeneID, l’utilitaire elink d’Entrez a été utilisé pour récupérer les publications relatives à chacun des gènes. Ces utilitaires peuvent être utilisés directement en accédant à une certaine URL (comme ici), ou bien sont disponibles dans des librairies R, python ou encore comme commande linux dans un shell.

Le script suivant a été utilisé pour récupérer, à partir des GeneID, les identifiants PubMed relatifs à une séquence :

>[!WARNING]   
> Ceci n’est PAS à réaliser → récupérez plutôt l’ensemble des résultats sur gitlab.

```python
#!/bin/env python

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
```

Nous obtenons donc le fichier associant les identifiants de gène aux identifiants PubMed :
```bash
curl -s https://gitlab.com/rbarriot/data/-/raw/main/M2.GDNS-AP/bnumber.PMID.tsv > generated.data/bnumber.PMID.tsv
head generated.data/bnumber.PMID.tsv
```

Utiliser la même technique que pour les mots-clés UniProt afin de générer le fichier au format JSON correspondant qui associe à chaque PMID la liste des gènes (bnumbers) cités dans la publication.

### 1.2 Données supplémentaires sur la localisation des gènes sur le chromosome

Par la suite, nous pourrons aussi être intéressés par la localisation des gènes sur le chromosome.

Effectuez un travail similaire à partir du fichier d’annotation GenBank afin d’extraire la localisation des gènes sur le chromosome et obtenir le tibble suivant :
| bnumber | rank  | strand | begin | end  |
| <chr>	  | <dbl> |	<chr>  | <dbl> | <dbl>|
| b0001	  | 1     | +	   | 189   | 255  |

### 1.3 Intégration dans Neo4J

<img src="M2BBS_GDNS-APG_schema_base_neo4j.png">
Schéma de la base de données

#### 1.3.1 Neo4j installation

Sites et documentations :
* https://neo4j.com/
* https://neo4j.com/developer/get-started/

Récupérer une image docker :
```bash
podman pull neo4j:5.23.0-community
```

Remarques :
* Autre téléchargements possibles sur https://neo4j.com/download-center/#releases
* Documentation sur l’utilisation de l’image : https://hub.docker.com/_/neo4j/

Démarrage du conteneur à partir du shell :
```bash
podman run --rm -it --name neo4j \
   -p 7474:7474 \
   -p 7687:7687 \
   -v ./neo4j.data:/data \
   -v ./neo4j.import:/import \
   -e NEO4J_AUTH=neo4j/omybioinfo \
   neo4j:5.23.0-community
```

Le processus est au premier plan donc pour arrêter le serveur il faut faire Ctrl + C dans le terminal.

Remarque : au lancement, le container vérifie à qui appartiennent les fichiers et répertoires (notamment /data et /var/lib/neo4j/import) dont il se sert, et les réattribue si besoin. Pour que vous puissiez créer/ajouter des fichiers dans le répertoire import, il faudra avoir les permissions.

Utilisation depuis le navigateur (vérifier le port renseigné lors de la précédente commande) : http://localhost:7474/

Au premier lancement du container, le mot de passe de l’utilisateur neo4j est redéfini omybioinfo (mais pas les suivantes une fois que la base de données est initialisée).

Aller sur http://localhost:7474 une première fois.

Remarque : le driver pour R neo4r ne semble plus maintenu et ne fonctionne pas avec les dernières versions de neo4j. Cela ne remets pas en question la constitution de la base de données son interrogation mais il faudra adapter les commandes pour un autre langage (python ?) ou une autre interface pour alimenter la base de données.

Ouverture d’un shell cq sur le conteneur en cours (depuis le shell) :
```bash
podman exec -it neo4j ./bin/cypher-shell
```

#### 1.3.2 Alimentation de la base de données

Librairie neo2R https://cran.r-project.org/web/packages/neo2R/index.html et https://github.com/patzaw/neo2r

```R
if (!require('neo2R')) { # client neo4j 
  install.packages('neo2R')
}

## Loading required package: neo2R

library(neo2R)
neodb <- startGraph(
  "http://localhost:7474",
  check = FALSE,
  username = "neo4j", 
  password = "omybioinfo",
  importPath = paste0(getwd(), "/neo4j.import"),
  .opts = list(ssl_verifypeer=0)
)
cq = function(query, neo4j=neodb, ...) cypher(neo4j, query, arraysAsStrings = F, ...)
```

##### 1.3.2.1 Nettoyage

Attention, cette partie nettoyage permet de supprimer tous les sommets et relations qu’il y a dans la base de données pour recommencer les commandes. Elle ne supprime pas les étiquettes ou types déjà rencontrés (ex: Gene)

Pour supprimer tous les liens et sommets
```R
'MATCH ()-[r]-() DELETE(r)' %>% cq
'MATCH ()-[r]-() RETURN count(r)' %>% cq %>% unlist
'MATCH (n) DELETE(n)' %>% cq
'MATCH (n) RETURN count(n)' %>% cq %>% unlist
```

##### 1.3.2.2 Création du graphe

###### 1.3.2.2.1 Genes

Copier le fichier dans le répertoire d’import de neo4j s’il n’est pas au bon endroit.
```bash
# cp generated.data/ncbi.bnumber.location.tsv neo4j.import/
head neo4j.import/ncbi.bnumber.location.tsv
```

Load CSV
```R
# IMPORT
"
LOAD CSV WITH HEADERS FROM 'file:///ncbi.bnumber.location.tsv' AS row FIELDTERMINATOR '\t' 
CREATE (n:Gene)
SET n = row,
 n.id = row.bnumber,
 n.organism = toInteger('511145'),
 n.rank = toInteger(row.rank),
 n.strand = row.strand,
 n.begin = toInteger(row.begin),
 n.end = toInteger(row.end)
"  %>% cq
```

Vérification
```R
'MATCH (n:Gene) RETURN count(n)' %>% cq %>% unlist
```

La sortie n’est pas très pratique. Nous allons créer un raccourci pour la reformater sous forme de tibble :
```R
cq2tb <- function(res) res$n %>% simplify2array %>% t 
"MATCH (n:Gene) RETURN n LIMIT 2" %>% cq %>% cq2tb
```

Index sinon l’ajout de liens et les requêtes peuvent être très longs
```R
'CREATE INDEX node_gene_id FOR (n:Gene) ON (n.id)' %>% cq
```

###### 1.3.2.2.2 Keywords et liens Keyword → Gene

Load CSV
```R
keywords %>% 
  select(keyword) %>%
  unique %>%
  write_csv("neo4j.import/uniprot.keywords.csv")
'LOAD CSV WITH HEADERS FROM "file:///uniprot.keywords.csv" AS row 
CREATE (n:Keyword)
SET n = row,
 n.id = row.keyword
'  %>% cq
'CREATE INDEX node_keyword_id FOR (n:Keyword) ON (n.id)' %>% cq
```

Vérification
```R
'MATCH (n:Keyword) RETURN count(n)' %>% cq %>% unlist
```

Liens Keyword → Gene
```R
keywords %>% write_csv("neo4j.import/uniprot.keywords.genes.csv")
```

Import du CSV
```R
"MATCH (:Keyword)-[r:describes]->(:Gene) DELETE r" %>% cq
"
LOAD CSV WITH HEADERS FROM 'file:///uniprot.keywords.genes.csv' AS line 
MATCH (k:Keyword),(g:Gene) 
WHERE k.id=line.keyword AND g.id=line.bnumber
WITH k,g 
MERGE (k)-[:describes]->(g)
" %>% cq
```

Vérification
```R
"MATCH (:Keyword)-[r:describes]->(:Gene) RETURN count(r)" %>% cq %>% unlist
```

###### 1.3.2.2.3 PubMed

```R
pmids <- read_tsv("generated.data/bnumber.PMID.tsv", col_types = "cc")
pmids %>% select(PMID) %>% unique %>% write_csv("neo4j.import/ncbi.pmid.csv")
pmids %>% rename(gene_id = bnumber) %>% write_csv("neo4j.import/ncbi.pmid.genes.csv")
```

Import du fichier CSV
```R
"MATCH (:PubMed)-[r:cites]->(:Gene) DELETE r" %>% cq
'MATCH (n:PubMed) DELETE n' %>% cq
'LOAD CSV WITH HEADERS FROM "file:///ncbi.pmid.csv" AS row 
CREATE (n:PubMed)
SET n = row,
 n.id = row.PMID
'  %>% cq
'CREATE INDEX node_pubmed_id FOR (n:PubMed) ON (n.id)' %>% cq
```

Vérification
```R
'MATCH (n:PubMed) RETURN count(n)' %>% cq %>% unlist
```

Import du fichier CSV
```R
"MATCH (:PubMed)-[r:cites]->(:Gene) DELETE r" %>% cq
"
LOAD CSV WITH HEADERS FROM 'file:///ncbi.pmid.genes.csv' AS line 
MATCH (p:PubMed),(g:Gene) 
WHERE p.id=line.PMID AND g.id=line.gene_id
WITH p,g 
MERGE (p)-[:cites]->(g)
" %>% cq
```

Vérification
```R
"MATCH (:PubMed)-[r:cites]->(:Gene) RETURN count(r)" %>% cq %>% unlist
```

Complétez cette partie afin de charger dans la base de données :
* les domaines InterPro,
* les pathways,
* les unités de transcription (TUs),

ainsi que leurs liens avec les gènes.

Résultats attendus :
* Vérification pour les domaines InterPro.
```R
"MATCH (n:InterPro) RETURN count(n)" %>% cq %>% unlist
```
* Liens InterPro → Genes
```R
"MATCH (n:InterPro)-[r:harbored_by]->(:Gene) RETURN count(r)" %>% cq %>% unlist
```
* Vérification des TUs
```R
"MATCH (n:TU) RETURN count(n)" %>% cq %>% unlist
```
* Vérification des liens TU → Gene
```R
"MATCH (n:TU)-[r:harbors]->(:Gene) RETURN count(r)" %>% cq %>% unlist
```
* Vérification des pathways
```R
"MATCH (n:Pathway) RETURN count(n)" %>% cq %>% unlist
```
* Vérification des liens Pathway → Gene
```R
"MATCH (n:Pathway)-[r:requires]->(:Gene) RETURN count(r)" %>% cq %>% unlist
```

#### 1.3.3 Adaptation du script pour l’utilisation de neo4j

Modules python nécessaire :
```bash
mamba install py2neo monotonic packaging
```

Créer un nouveau script en adaptant le précédent qui utilisera, non plus les fichiers au format JSON, mais une connexion au serveur Neo4J pour effectuer les recherches d’enrichissement.

Pour pouvoir construire la requête Neo4J qui identifiera les ensembles cibles, nous avons besoin du type de sommets (PubMed, TU, Pathways, …) ainsi que de l’identifiant taxonomique (en prévision de l’ajout d’autres organismes).

Exemple d’utilisation avec set.01.txt :
```bash
./scripts/blastsets.neo4j.py -q query.sets/set.01.txt -t Pathway --species 511145 -c
```

### 1.4 Ajout de Gene Ontology

Récupération de la dernière version au format OBO :
```bash
curl http://current.geneontology.org/ontology/go-basic.obo -o downloads/go-basic.obo
```

Réutilisation du code de la librairie faite en graphes pour générer les fichiers correspondants aux sommets et arcs à charger dans neo4j à partir du script GeneOntology.py
```bash
curl -s https://gitlab.com/rbarriot/data/-/raw/main/M2.GDNS-AP/GeneOntology.py > scripts/GeneOntology.py
python scripts/GeneOntology.py downloads/go-basic.obo nodes > neo4j.import/go.nodes.tsv
python scripts/GeneOntology.py downloads/go-basic.obo edges is_a > neo4j.import/go.is_a.edges.tsv
python scripts/GeneOntology.py downloads/go-basic.obo edges part_of > neo4j.import/go.part_of.edges.tsv
```

Sommets
```R
"MATCH (:GOTerm)-[r]-() DELETE r" %>% cq
"MATCH (n:GOTerm) DELETE n" %>% cq
"
LOAD CSV WITH HEADERS FROM 'file:///go.nodes.tsv' AS row FIELDTERMINATOR '\t' 
CREATE (n:GOTerm)
SET n.id = row.id,
n.name = row.desc,
n.desc  = row.def,
n.namespace = row.namespace
"  %>% cq
```

Vérification
```R
"MATCH (n:GOTerm) RETURN count(n)" %>% cq %>% unlist
```

Index
```R
"CREATE INDEX node_go_id FOR (n:GOTerm) ON (n.id)" %>% cq
```

Extait du contenu :
```R
"MATCH (n:GOTerm) RETURN n LIMIT 6" %>% cq %>% cq2tb
```

Relations is a
```R
"
LOAD CSV WITH HEADERS FROM 'file:///go.is_a.edges.tsv' AS line FIELDTERMINATOR '\t' 
MATCH (t1:GOTerm),(t2:GOTerm) 
WHERE t1.id=line.term1 AND t2.id=line.term2
WITH t1,t2 
MERGE (t2)-[r:is_a]->(t1)
" %>% cq
```

Vérification
```R
"MATCH ()-[r:is_a]->() RETURN count(r)" %>% cq %>% unlist
```

Relation part of
```R
"
LOAD CSV WITH HEADERS FROM 'file:///go.part_of.edges.tsv' AS line FIELDTERMINATOR '\t' 
MATCH (t1:GOTerm),(t2:GOTerm) 
WHERE t1.id=line.term1 AND t2.id=line.term2
WITH t1,t2 
MERGE (t2)-[r:part_of]->(t1)
" %>% cq
```

Vérification
```R
"MATCH ()-[r:part_of]->() RETURN count(r)" %>% cq %>% unlist
```

Gene associations
```R
GOTerms <- uniprot %>% 
  select(uniprotID=Entry, GOTerm=`Gene Ontology IDs`) %>%
  right_join(mapping, by = join_by(uniprotID)) %>% # right join to remove those without bnumber
  separate_rows(GOTerm, sep='; ') %>%
  select(bnumber, GOTerm) %>%
  arrange(bnumber)
GOTerms
```

Save CSV file for Neo4J
```R
GOTerms %>% write_csv("neo4j.import/uniprot.GOTerm.bnumber.csv")
```

Merge into graph
```R
"
LOAD CSV WITH HEADERS FROM 'file:///uniprot.GOTerm.bnumber.csv' AS line
MATCH (t:GOTerm),(g:Gene) 
WHERE t.id=line.GOTerm AND g.id=line.bnumber
WITH t,g
MERGE (t)-[:annotates]->(g)
" %>% cq
```

Vérification
```R
"MATCH ()-[r:annotates]->() RETURN count(r)" %>% cq %>% unlist
```

#### 1.4.1 Adaptation du script en conséquence

Adapter le script précédent afin qu’il compare les ensembles cibles correspondants aux GOTerm (en n’oubliant pas les associations GOterm → GeneProduct implicites).

Exemple d’utilisation :
```bash
./scripts/blastsets.neo4j.py -q query.sets/set.01.txt -t GOTerm --species 511145 -c
```

#### 1.4.2 Distance sémantique entre GOTerms sur-représentés

Analysez les résultats du script (GOTerm et p-valeurs) avec REVIGO : http://revigo.irb.hr/

## 2 Approches par graphe ou matrice de distances/dissemblances/similarités/probabilités

### 2.1 Génération de matrices de similarité

A partir des liens Keyword, Pathway, PubMed et GOTerm, il s’agit de générer des matrices pivots individus(bnumber)–variables(Keyword ou pahtways ou …) pour ensuite calculer la distance/similarité entre paires de gènes.

Ajout d’une colonne avec des 1 pour les association présentes dans les données de départ :
```R
keyword.mat <- keywords %>% 
  unique %>% 
  mutate(asso=1)
keyword.mat
```

Génération d’une table pivot (les valeurs de la colonne keyword deviennent des colonnes) :
```R
keyword.tab <- keyword.mat %>% 
  pivot_wider(names_from = keyword, values_from = asso, values_fill = 0) %>% 
  as.data.frame
keyword.tab %>% head
```

Transfert des identifiants de gènes de la colonne bnumber aux noms de lignes (nécessaires pour qu’ils soient reportés sur les noms des lignes/colonnes) de la matrice de distance :
```R
rownames(keyword.tab) <- keyword.tab$bnumber
keyword.tab %>% head
```

Calcul de la matrice de distance (entre paires de lignes/gènes) avec la méthode binary (cf. ?dist) :
```R
keyword.tab <- keyword.tab %>% dplyr::select(-bnumber)
keyword.dist <- keyword.tab %>% dist(method='binary')
keyword.dist[1:10]
```

Sous forme de scores à la StringDB
```R
keyword.scores <- round( (1-keyword.dist) * 1000) %>% as.matrix
keyword.scores[1:10,1:10]
```

### 2.2 Client R à STRINGdb

Client R Bioconductor :
```R
library(STRINGdb)
```

Autres librairies utiles pour la gestion et l’affichage de tables :
```R
# mise en forme des tableaux:
library(kableExtra) # nicer tables

kabex <- function(o) o %>% kable(format="html", escape=F) %>% kable_styling(bootstrap_options = c("striped"))
library(DT) 
dted <- function(o) o %>% datatable(rownames=F, filter="top", options=list(pageLength = 10), escape=F)
```

Connexion au serveur de STRING (NCBI taxon id de Ecoli est 511145 : https://string-db.org/cgi/input.pl?input_page_active_form=organism )
```R
confidence_threshold <- 333
stringdb <- STRINGdb$new(version='12', species=511145, score_threshold=confidence_threshold, input_directory='downloads')
```

Téléchargement du génome/protéome :
```R
proteins <- stringdb$get_proteins() %>% as_tibble
proteins %>% dted
```

Fichier récupéré
```bash
zcat downloads/511145.protein.info.v12.0.txt.gz | head | column -ts $'\t'
```

Chargement de set.01.txt pour analyse.
```R
s1 <- scan('query.sets/set.01.txt', character()) 
s1
```

Mapping des identifiants vers ceux de StringDB
```R
s1.mapped <- stringdb$mp(s1)
s1.mapped
```

Plot du graphe correspondant
```R
stringdb$plot_network(s1.mapped)
```

Enrichment
```R
enrichment <- stringdb$get_enrichment(s1.mapped)
enrichment %>% dted
```

Sommets voisins d’un ou plusieurs sommets donnés
```R
stringdb$get_neighbors(s1.mapped[1:3])
```

Interactions entre un ensemble de sommets donnés
```R
stringdb$get_interactions(s1.mapped)
```

### 2.3 Données de coexpression à partir des données complètes

Nous allons utiliser les données d’expression déjà intégrées dans STRINGdb. Il est bien sûr possible de calculer un score de co-expression entre les paires de gènes à partir de données de microarray et/ou de RNA-Seq.

#### 2.3.1 STRINGdb detailed links → neo4j

Il faut passer par la page de téléchargement sur https://string-db.org/cgi/download.pl en restreignant l’organisme d’intérêt.

Le fichier complet faisant ~150Go (au 15/10/2021), il est conseillé de d’abord sélectionner un organisme avant le téléchargement (~10Mo pour E. coli). Une copie est disponible sur silico à récupérer sur http://silico.biotoul.fr/enseignement/m2bbs/gdns-ap/

Qu’y a-t-il dedans ?
```bash
zcat downloads/511145.protein.links.detailed.v12.0.txt.gz | head | column -t
```

Chargement sous forme de tibble :
```R
links.detailed <- read_delim("downloads/511145.protein.links.detailed.v12.0.txt.gz", delim=" ", col_types = "ccnnnnnnnn")
links.detailed
```

C’est la colonne coexpression qui nous intéresse pour commencer.

Génération des fichiers CSV pour l’import
```R
links.detailed %>%
  filter(coexpression>0 & protein1 < protein2) %>%
  select(protein1, protein2, coexpression) %>%
  mutate(organism=str_extract(protein1, '^\\d+'), id1=str_extract(protein1, 'b\\d+'), id2=str_extract(protein2, 'b\\d+')) %>%
  select(organism:id2,coexpression) %>%
  write_csv("neo4j.import/string.coexpression.csv")
```

Import dans neo4j
```R
'
LOAD CSV WITH HEADERS FROM "file:///string.coexpression.csv" AS line
MATCH (g1:Gene),(g2:Gene)
 WHERE g1.id=line.id1 AND g2.id=line.id2
WITH g1,g2, toInteger(line.coexpression) AS value 
MERGE (g1)-[r:STRINGdb]-(g2)
 SET r.coexpression=value
'  %>% cq
```

Vérification
```R
'MATCH ()-[r:STRINGdb]-() RETURN count(r)' %>% cq %>% unlist
```

Extraction des liens de coexpression d’au moins 0.4
```R
coex <- "MATCH (g1:Gene)-[r:STRINGdb]-(g2:Gene) WHERE r.coexpression>=995 RETURN g1.id, g2.id, r.coexpression" %>% cq 
tibble(id1=coex$g1.id, id2=coex$g2.id, coexpression=coex$r.coexpression)
```

Extraction d’un sous-graphe : on utilise ici neo4r et la fonction fournie (call_neo4j) avec le paramètre type="graph" afin de récupérer un sous-graphe.
```R
g.coexpr <- 'MATCH p = ()-[r:STRINGdb]-() WHERE r.coexpression>=995 RETURN p'  %>% cq(result="graph")
```

Manipulation pour igraph pour récupérer les propriétés des sommets sous forme de tibble :
```R
g.coexpr$nodes <-  g.coexpr$nodes |> 
  map(\(node) c(neo_id=node$id, node$properties)) %>% 
  bind_rows()
g.coexpr$nodes %>% 
  dted
```

De même pour les arcs/arêtes du graphe à l’aide de la fonction unnest_relationships pour le passage dans igraph : besoin de réordonner les liens :
```R
g.coexpr$relationships <- g.coexpr$relationships |> 
  map(\(edge) c(startNode=edge$startNode, endNode=edge$endNode, type=edge$type, edge$properties)) %>% 
  bind_rows

g.coexpr$relationships %>% 
  dted
```

Chargement de la librairies et positionnement des valeurs par défaut pour certains paramètres :
```R
library(igraph)

igraph.options(vertex.color=NA)

## Warning: `igraph.options()` was deprecated in igraph 2.0.0.
## ℹ Please use `igraph_options()` instead.

igraph.options(vertex.label.cex=.6) # font size
igraph.options(vertex.label.family='sans')
igraph.options(vertex.size=2)
igraph.options(edge.label.cex=.6)
igraph.options(edge.label.family='sans')
```

Création du graphe dans igraph à partir des sommets et des relations (remarque : la première colonne du df passé comme vertices est considérée comme l’identifiant des sommes, donc neo_id ici) :
```R
g.coexpr = graph_from_data_frame(d=g.coexpr$relationships, directed=FALSE, vertices = g.coexpr$nodes)
g.coexpr
```

Plot
```R
plot(g.coexpr, vertex.label=NA, main='coexpression')
```

### 2.4 Autres scores StringDB à intégrer

Faire de même avec les liens correspondant aux
* interactions protéines–protéines → colonne experimental,
* les associations basées sur la conservation des paires de gènes dans le voisinage sur le chromosome dans d’autres génomes → neighborhood,
* les associations basées sur la littérature bio-médicale → textmining, (que l’on pourra comparer aux identifiants PubMed→Gene )
* les associations basées sur les annotations → database, (que l’on pourra comparer aux liens Keyword→Gene, Pathway→Gene et GOTerm→Gene)
* le score combiné → combined_score

### 2.5 Combinaison de matrices de similarité

Calcul du score combiné par STRINGdb (source: https://string-db.org/help/faq/#how-are-the-scores-computed).

En utilisant la librairie igraph et après avoir chargé le graphe g

avec les liens de coexpression, neighborhood et experiment à partir de neo4j :
```R
prior <- 0.041
no_prior <- function(x, prior = 0.041) (ifelse(is.na(x), 0, x) / 1000 - prior) / (1-prior)
s_coexp_nop <- no_prior(E(g)$coexpression)
s_ppi_nop <- no_prior(E(g)$ppi)
s_neighborhood_nop <- no_prior(E(g)$neighborhood)
s_tot_nop <- 1 - (1 - s_coexp_nop) * (1 - s_ppi_nop) * (1 - s_neighborhood_nop)
E(g)$combined_score <- round(1000 * (s_tot_nop + prior *(1 - s_tot_nop)))
```

### 2.6 Prioritisation de gènes

Gènes “training”. Choix du cycle des citrates (pathway TCA). Récupération des gènes impliqués dans le pathway
```R
training.genes <- "MATCH (Pathway {id: 'TCA'})-[:requires]->(g: Gene {organism: 511145}) return g.bnumber " %>% cq  %>% .$g %>% unlist
training.genes
```

Gènes “candidats” (ensemble du génome) :
```R
candidates <- "MATCH (g:Gene {organism: 511145}) RETURN g.bnumber" %>% cq %>% .$g 
candidates
```

Distance d’un gène à l’ensemble de référence
```R
score <- function(gene.ids, ref.genes, datasource) {
  ref_str <- paste("'", ref.genes,"'", sep = '', collapse = ',')
  sapply(gene.ids, function(gene.id) {
    query <- paste0("MATCH (g1:Gene {id: '", gene.id, "', organism: 511145})-[r:STRINGdb]-(g2:Gene {organism: 511145}) WHERE g2.id IN [",ref_str,"] RETURN r.",datasource)
    res <- query %>% cq %>% .$r
    ifelse(is.null(res), 0, mean(replace_na(res, 0)))
  })
}
# test
gene.ids <- 'b0001'
gene.ids <- 'b0116'
ref.genes <- training.genes
datasource <- 'coexpression'
score(gene.ids, training.genes, datasource)
```

Application à l’ensemble des gènes
```R
scores <- tibble(candidate=candidates) %>%
  mutate(
    coexpression = score(gene.id=candidate, ref.genes=training.genes, datasource='coexpression'),
    experimental  = score(gene.id=candidate, ref.genes=training.genes, datasource='experimental'),
    neighborhood  = score(gene.id=candidate, ref.genes=training.genes, datasource='neighborhood'),
    textmining  = score(gene.id=candidate, ref.genes=training.genes, datasource='textmining'),
    database  = score(gene.id=candidate, ref.genes=training.genes, datasource='database')
         )
scores <- scores  %>%  replace(is.na(.), 0)
scores
```

Plots
```R
scores %>%
  ggplot(aes(x=coexpression, y=neighborhood)) +
  geom_point(alpha=0.2, aes(color=candidate %in% training.genes, shape=candidate %in% training.genes))
```

ACP pour l’affichage
```R
pca = scores %>%
  dplyr::select(-candidate) %>%
  prcomp
pca %>% summary
```

Visualisation
```R
library(factoextra)
pca %>% fviz_pca_ind(col.ind = scores$candidate %in% training.genes, addEllipses = T, alpha.ind = .2, geom='point')
```

On peut déjà observer une petite distinction entre les gènes de référence et les autres.

Contribution des sources de données
```R
pca %>% fviz_pca_biplot(col.ind = scores$candidate %in% training.genes)
```

Analyse discriminante linéaire
```R
library(MASS)
mat <- scores %>% dplyr::select(-candidate) %>% as.matrix
model <- lda(x=mat, grouping = scores$candidate %in% training.genes)
model
```

Remarque : attention lors de l’utilisation des librairies MASS et tidyverse, notamment à l’ordre de leur chargement car les 2 fournissent la fonction select. Par exemple ici, nous avons d’abord chargé tidyverse puis MASS. Les bouts de code précédents qui utilise select ne vont plus fonctionner sans préciser qu’il s’agit de dplyr::select que l’on souhaite appeler et non celle de MASS.

Distribution des scores sous forme de boîtes à moustache pour les 2 classes :
```R
mat %>% 
  scale( center=T, scale=F ) %*% model$scaling %>% 
  as_tibble %>%
  mutate(gene.id=scores$candidate) %>%
  ggplot(aes(x=LD1, y=gene.id %in% training.genes, color=gene.id %in% training.genes, shape=gene.id %in% training.genes)) +
  geom_violin() +
  geom_boxplot(varwidth = T) +
  geom_jitter(height = 0.2, alpha=0.1, color='grey') +
  theme_light()
```

Density
```R
mat %>% 
  scale( center=T, scale=F ) %*% model$scaling %>% 
  as_tibble %>%
  mutate(gene.id=scores$candidate) %>%
  ggplot(aes(LD1)) + #, color=gene.id %in% training.genes, shape=gene.id %in% training.genes)) +
  geom_density()
```

On remarquera que faire >4k requête (1 par gène candidat) prend du temps afin d’obtenir la moyenne (ou autre) entre le candidat et les gènes de référence.

Tentative d’optimisation.

Récupération des “scores” entre les candidats et les gènes de référence en une requête
```R
train.str <- paste("'", training.genes,"'", sep = '', collapse = ',')
datasources <- c('coexpression', 'experimental', 'neighborhood', 'textmining', 'database' )
ds.str <- paste(" r.", datasources, sep='', collapse=',')
tb <- paste0("MATCH (candidate:Gene {organism: 511145})-[r:STRINGdb]-(training:Gene {organism: 511145}) WHERE training.id IN [", train.str,"]  RETURN DISTINCT candidate.id, training.id, ", ds.str) %>% cq
tb
```

Replace NAs with 0s and compute mean (or other) by candidate
```R
scores <- tb %>% 
  dplyr::select(-training.id) %>%
  replace(is.na(.), 0) %>%
  pivot_longer(cols = paste0('r.',datasources), names_to='datasource') %>%
  group_by(candidate.id, datasource) %>%
  summarise(score=mean(value), .groups="drop") %>%
  pivot_wider(names_from = datasource, values_from = score) %>% 
  ungroup
scores
```

Analyse discriminante linéaire
```R
mat <- scores %>% dplyr::select(-candidate.id) %>% as.matrix
model <- lda(x=mat, grouping = scores$candidate.id %in% training.genes)
model
```

Distribution des scores sous forme de boîtes à moustache pour les 2 classes :
```R
mat %>% 
  scale( center=T, scale=F ) %*% model$scaling %>% 
  as_tibble %>%
  mutate(gene.id=scores$candidate.id, training=gene.id %in% training.genes) %>%
  ggplot(aes(x=LD1, y=training,  color=training, shape=training)) +
  geom_violin() +
  geom_jitter(alpha=.5, height = .2)
```
