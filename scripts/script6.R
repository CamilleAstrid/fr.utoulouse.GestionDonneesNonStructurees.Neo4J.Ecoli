#==========================================#
# SCRIPT6
#==========================================#
# Script pour interagir avec la base de données Neo4j


setwd("~/GitHub/fr.utoulouse.GestionDonneesNonStructurees.Neo4J.Ecoli")

#------------------------------------------#
# Chargement des librairies
#------------------------------------------#
library(neo2R)
library(tidyverse)
library(reticulate)
library(jsonlite)

# On nettoie le cache
rm(list=ls())
gc()

#------------------------------------------#
# Connexion à la base de données Neo4j
#------------------------------------------#
neodb <- startGraph(
  "http://localhost:7474",
  check = FALSE,
  username = "neo4j", 
  password = "omybioinfo",
  importPath = paste0(getwd(), "../neo4j.import"),
  .opts = list(ssl_verifypeer=0)
)

#------------------------------------------#
# Création des fonctions utilitaires
#------------------------------------------#
cq = function(query, neo4j=neodb, ...) cypher(neo4j, query, arraysAsStrings = F, ...)
cq2tb = function(res) res$n %>% simplify2array %>% t
kill_all_transactions <- function() {
  txs <- "
  SHOW TRANSACTIONS
  YIELD transactionId
  RETURN transactionId
  " %>% cq
  
  tx_ids <- txs$transactionId
  
  for (id in tx_ids) {
    cmd <- sprintf("TERMINATE TRANSACTION '%s'", id)
    cat(">>", cmd, "\n")
    try(cmd %>% cq, silent = TRUE)
  }
}


#------------------------------------------#
# NETTOYAGE
#------------------------------------------#
# Suppression des transactions actives
kill_all_transactions()
'SHOW TRANSACTIONS' %>% cq
# Suppression des nodes et des relations existants
'MATCH (n) DETACH DELETE n' %>% cq
'MATCH ()-[r]-() RETURN count(r)' %>% cq %>% unlist
'MATCH (n) RETURN count(n)' %>% cq %>% unlist
# Suppression des index existants
'DROP INDEX node_gene_id IF EXISTS' %>% cq
'DROP INDEX node_keyword_id IF EXISTS' %>% cq
'DROP INDEX node_pubmed_id IF EXISTS' %>% cq
'DROP INDEX node_interpro_id IF EXISTS' %>% cq
'DROP INDEX node_pathway_id IF EXISTS' %>% cq
'DROP INDEX node_tu_id IF EXISTS' %>% cq
'DROP INDEX node_go_id IF EXISTS' %>% cq

#==========================================#
# Alimentation de la base de données
#==========================================#

#------------------------------------------#
# Genes
#------------------------------------------#
# Import
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

# Vérification
'MATCH (n:Gene) RETURN count(n)' %>% cq %>% unlist
# sous forme de tibble
"MATCH (n:Gene) RETURN n LIMIT 2" %>% cq %>% cq2tb

# Création d'un index
'CREATE INDEX node_gene_id FOR (n:Gene) ON (n.id)' %>% cq

#------------------------------------------#
# Keywords et liens Keyword → Gene
#------------------------------------------#
# Import
'LOAD CSV WITH HEADERS FROM "file:///uniprot.keywords.csv" AS row 
CREATE (n:Keyword)
SET n = row,
 n.id = row.keyword
'  %>% cq
'CREATE INDEX node_keyword_id FOR (n:Keyword) ON (n.id)' %>% cq

# Vérification
'MATCH (n:Keyword) RETURN count(n)' %>% cq %>% unlist

# Import du CSV
"
LOAD CSV WITH HEADERS FROM 'file:///uniprot.keywords.genes.csv' AS line 
MATCH (k:Keyword),(g:Gene) 
WHERE k.id=line.keyword AND g.id=line.bnumber
WITH k,g 
MERGE (k)-[:describes]->(g)
" %>% cq

# Vérification
"MATCH (:Keyword)-[r:describes]->(:Gene) RETURN count(r)" %>% cq %>% unlist

#------------------------------------------#
# PubMed
#------------------------------------------#
# Import
pmids <- read_tsv("generated.data/bnumber.PMID.tsv", col_types = "cc")
pmids %>% select(PMID) %>% unique %>% write_csv("neo4j.import/ncbi.pmid.csv")
pmids %>% rename(gene_id = bnumber) %>% write_csv("neo4j.import/ncbi.pmid.genes.csv")

'LOAD CSV WITH HEADERS FROM "file:///ncbi.pmid.csv" AS row 
CREATE (n:PubMed)
SET n = row,
 n.id = row.PMID
'  %>% cq
'CREATE INDEX node_pubmed_id FOR (n:PubMed) ON (n.id)' %>% cq

# Vérification
'MATCH (n:PubMed) RETURN count(n)' %>% cq %>% unlist

# Import du fichier CSV
# modification pour gérer les grosses quantités de données
"
CALL {
  LOAD CSV WITH HEADERS FROM 'file:///ncbi.pmid.genes.csv' AS line
  MATCH (p:PubMed {id: line.PMID}),
        (g:Gene {id: line.gene_id})
  MERGE (p)-[:cites]->(g)
} IN TRANSACTIONS OF 10000 ROWS
" %>% cq

# Vérification
"MATCH (:PubMed)-[r:cites]->(:Gene) RETURN count(r)" %>% cq %>% unlist

# Modification des sommets dans neo4j
'
LOAD CSV WITH HEADERS FROM "file:///pmid.selected.tsv" AS line FIELDTERMINATOR "\t" 
MATCH (p:PubMed) 
WHERE p.id = line.PMID
WITH p, line
SET p.name = line.title,
  p.title = line.title,
  p.year = toInteger(line.year),
  p.journal = line.journal,
  p.authors = line.authors,
  p.url = "https://pubmed.ncbi.nlm.nih.gov/"+line.PMID
' %>% cq

#------------------------------------------#
# InterPro : nodes+links InterPro → Gene
#------------------------------------------#

# Préparation des CSV à partir du JSON de référence
interpro_sets <- fromJSON("reference.sets/uniprot.InterPro.sets.json")
interpro_sets <- as_tibble(interpro_sets)

interpro_nodes <- interpro_sets %>%
  select(id, desc) %>%
  distinct()

interpro_edges <- interpro_sets %>%
  select(id, elements) %>%
  unnest_longer(elements) %>%
  rename(bnumber = elements) %>%
  distinct()

write_csv(interpro_nodes, "neo4j.import/uniprot.interpro.csv")
write_csv(interpro_edges, "neo4j.import/uniprot.interpro.genes.csv")

# Création des nodes InterPro
'LOAD CSV WITH HEADERS FROM "file:///uniprot.interpro.csv" AS row 
CREATE (n:InterPro)
SET n = row,
 n.id = row.id,
 n.desc = row.desc
' %>% cq
'CREATE INDEX node_interpro_id FOR (n:InterPro) ON (n.id)' %>% cq

# Vérification
"MATCH (n:InterPro) RETURN count(n)" %>% cq %>% unlist

# Création des liens InterPro → Gene
"
LOAD CSV WITH HEADERS FROM 'file:///uniprot.interpro.genes.csv' AS line 
MATCH (i:InterPro),(g:Gene) 
WHERE i.id = line.id AND g.id = line.bnumber
WITH i,g 
MERGE (i)-[:harbored_by]->(g)
" %>% cq

# Vérification
"MATCH (n:InterPro)-[r:harbored_by]->(:Gene) RETURN count(r)" %>% cq %>% unlist

#------------------------------------------#
# Pathways : nodes+links Pathway → Gene
#------------------------------------------#

# Préparation des CSV à partir du JSON de référence
pathway_sets <- fromJSON("reference.sets/ecocyc.Pathways.sets.json")
pathway_sets <- as_tibble(pathway_sets)

pathway_nodes <- pathway_sets %>%
  select(id, desc) %>%
  distinct()

pathway_edges <- pathway_sets %>%
  select(id, elements) %>%
  unnest_longer(elements) %>%
  rename(bnumber = elements) %>%
  distinct()

write_csv(pathway_nodes, "neo4j.import/ecocyc.pathways.csv")
write_csv(pathway_edges, "neo4j.import/ecocyc.pathways.genes.csv")

# Création des nodes Pathway
'LOAD CSV WITH HEADERS FROM "file:///ecocyc.pathways.csv" AS row 
CREATE (n:Pathway)
SET n = row,
 n.id = row.id,
 n.desc = row.desc
' %>% cq
'CREATE INDEX node_pathway_id FOR (n:Pathway) ON (n.id)' %>% cq

# Vérification
"MATCH (n:Pathway) RETURN count(n)" %>% cq %>% unlist

# Création des liens Pathway → Gene
"
LOAD CSV WITH HEADERS FROM 'file:///ecocyc.pathways.genes.csv' AS line 
MATCH (p:Pathway),(g:Gene) 
WHERE p.id = line.id AND g.id = line.bnumber
WITH p,g 
MERGE (p)-[:requires]->(g)
" %>% cq

# Vérification
"MATCH (n:Pathway)-[r:requires]->(:Gene) RETURN count(r)" %>% cq %>% unlist

#------------------------------------------#
# TUs : nodes+links TU → Gene
#------------------------------------------#

# Préparation des CSV à partir du JSON de référence
tu_sets <- fromJSON("reference.sets/ecocyc.TUs.sets.json")
tu_sets <- as_tibble(tu_sets)

tu_nodes <- tu_sets %>%
  select(id, desc) %>%
  distinct()

tu_edges <- tu_sets %>%
  select(id, elements) %>%
  unnest_longer(elements) %>%
  rename(bnumber = elements) %>%
  distinct()

write_csv(tu_nodes, "neo4j.import/ecocyc.tus.csv")
write_csv(tu_edges, "neo4j.import/ecocyc.tus.genes.csv")

# Création des nodes TU
'LOAD CSV WITH HEADERS FROM "file:///ecocyc.tus.csv" AS row 
CREATE (n:TU)
SET n = row,
 n.id = row.id,
 n.desc = row.desc
' %>% cq
'CREATE INDEX node_tu_id FOR (n:TU) ON (n.id)' %>% cq

# Vérification
"MATCH (n:TU) RETURN count(n)" %>% cq %>% unlist

# Création des liens TU → Gene
"
LOAD CSV WITH HEADERS FROM 'file:///ecocyc.tus.genes.csv' AS line 
MATCH (t:TU),(g:Gene) 
WHERE t.id = line.id AND g.id = line.bnumber
WITH t,g 
MERGE (t)-[:harbors]->(g)
" %>% cq

# Vérification
"MATCH (n:TU)-[r:harbors]->(:Gene) RETURN count(r)" %>% cq %>% unlist


#==========================================#
# Ajout de la Gene Ontology
#==========================================#

#------------------------------------------#
# Sommets
#------------------------------------------#
# Import
"
LOAD CSV WITH HEADERS FROM 'file:///go.nodes.tsv' AS row FIELDTERMINATOR '\t' 
CREATE (n:GOTerm)
SET n.id = row.id,
n.name = row.desc,
n.desc  = row.def,
n.namespace = row.namespace
"  %>% cq

# Vérification
"MATCH (n:GOTerm) RETURN count(n)" %>% cq %>% unlist

# Index
"CREATE INDEX node_go_id FOR (n:GOTerm) ON (n.id)" %>% cq

# Extait du contenu :
"MATCH (n:GOTerm) RETURN n LIMIT 6" %>% cq %>% cq2tb

#------------------------------------------#
# Relations is a
#------------------------------------------#"
"
CALL {
  LOAD CSV WITH HEADERS FROM 'file:///go.is_a.edges.tsv' AS line FIELDTERMINATOR '\t' 
  MATCH (t1:GOTerm),(t2:GOTerm) 
  WHERE t1.id=line.term1 AND t2.id=line.term2
  WITH t1,t2 
  MERGE (t2)-[r:is_a]->(t1)
} IN TRANSACTIONS OF 10000 ROWS
" %>% cq

# Vérification
"MATCH ()-[r:is_a]->() RETURN count(r)" %>% cq %>% unlist

#------------------------------------------#
# Relation part of
#------------------------------------------#
"
LOAD CSV WITH HEADERS FROM 'file:///go.part_of.edges.tsv' AS line FIELDTERMINATOR '\t' 
MATCH (t1:GOTerm),(t2:GOTerm) 
WHERE t1.id=line.term1 AND t2.id=line.term2
WITH t1,t2 
MERGE (t2)-[r:part_of]->(t1)
" %>% cq

# Vérification
"MATCH ()-[r:part_of]->() RETURN count(r)" %>% cq %>% unlist

#------------------------------------------#
# Merge into graph
#------------------------------------------#
"
LOAD CSV WITH HEADERS FROM 'file:///uniprot.GOTerm.bnumber.csv' AS line
MATCH (t:GOTerm),(g:Gene) 
WHERE t.id=line.GOTerm AND g.id=line.bnumber
WITH t,g
MERGE (t)-[:annotates]->(g)
" %>% cq

# Vérification
"MATCH ()-[r:annotates]->() RETURN count(r)" %>% cq %>% unlist

#==========================================#
# GeneProduct : nodes+links GeneProduct → Gene
#==========================================#

# Chargement du protéome UniProt (même fichier que dans les scripts précédents)
uniprot <- read_tsv("downloads/uniprotkb_proteome_UP000000625.tsv.gz", show_col_types = FALSE)

# Mapping UniProt ↔ bnumber (Gene)
gp_mapping <- uniprot %>%
  select(uniprotID = Entry, names = `Gene Names (ordered locus)`) %>%
  mutate(bnumber = stringr::str_extract(names, "b\\d+")) %>%
  filter(!is.na(bnumber)) %>%
  distinct(uniprotID, bnumber)

# Nodes GeneProduct
gp_nodes <- gp_mapping %>%
  distinct(uniprotID) %>%
  transmute(gene_product_id = uniprotID)

# Relations GeneProduct → Gene (encoded_by)
gp_edges <- gp_mapping %>%
  transmute(gene_product_id = uniprotID,
                   bnumber = bnumber)

# Écriture des CSV pour Neo4j
write_csv(gp_nodes, "neo4j.import/uniprot.geneproduct.nodes.csv")
write_csv(gp_edges, "neo4j.import/uniprot.geneproduct.encoded_by.csv")

# Création des nodes GeneProduct
"
LOAD CSV WITH HEADERS FROM 'file:///uniprot.geneproduct.nodes.csv' AS row
CREATE (n:GeneProduct)
SET n = row,
    n.id = row.gene_product_id
" %>% cq

# Index sur GeneProduct.id (avec IF NOT EXISTS pour éviter les erreurs aux relances)
"CREATE INDEX node_geneproduct_id IF NOT EXISTS FOR (n:GeneProduct) ON (n.id)" %>% cq

# Vérification
"MATCH (n:GeneProduct) RETURN count(n)" %>% cq %>% unlist

# Création des liens GeneProduct → Gene via :encoded_by
"
LOAD CSV WITH HEADERS FROM 'file:///uniprot.geneproduct.encoded_by.csv' AS line
MATCH (gp:GeneProduct {id: line.gene_product_id}),
      (g:Gene        {id: line.bnumber})
MERGE (gp)-[:encoded_by]->(g)
" %>% cq

# Vérification
"MATCH (:GeneProduct)-[r:encoded_by]->(:Gene) RETURN count(r)" %>% cq %>% unlist


#==========================================#
# END SCRIPT6
#==========================================#
