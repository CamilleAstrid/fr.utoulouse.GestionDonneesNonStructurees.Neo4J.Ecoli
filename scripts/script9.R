#!/usr/bin/env Rscript

#==========================================#
# SCRIPT9
#==========================================#
# Script pour calculer la distance semantique entre GOTerms sur-representes
# et pour générer des graphes de coexpression à partir des données de STRINGdb


setwd("~/GitHub/fr.utoulouse.GestionDonneesNonStructurees.Neo4J.Ecoli")

#------------------------------------------#
# Chargement des librairies
#------------------------------------------#
library(neo2R)
library(tidyverse)
library(reticulate)
library(jsonlite)
if (!require("BiocManager", quietly = TRUE)){install.packages("BiocManager", repos='https://mirror.ibcp.fr/pub/CRAN/')}
if (!require("rrvgo", quietly = TRUE)){BiocManager::install("rrvgo", update = FALSE)}
library(rrvgo)
if (!require("ontologyIndex")) install.packages("ontologyIndex", repos='https://mirror.ibcp.fr/pub/CRAN/')
library(ontologyIndex)
if (!require("GOSemSim")) BiocManager::install("GOSemSim", update = FALSE)
library(GOSemSim)
if (!require("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db", update = FALSE)
library(org.Hs.eg.db)
if (!require("STRINGdb")) BiocManager::install("STRINGdb", update = FALSE)
library(STRINGdb)
if (!require("kableExtra")) install.packages("kableExtra", repos='https://mirror.ibcp.fr/pub/CRAN/')
library(kableExtra)
if (!require("htmlwidgets")) install.packages("htmlwidgets", repos='https://mirror.ibcp.fr/pub/CRAN/')
library(htmlwidgets)

#------------------------------------------#
# Connexion à la base de données Neo4j
#------------------------------------------#
neodb <- startGraph(
    "http://localhost:7474",
    check = FALSE,
    username = "neo4j", 
    password = "omybioinfo",
    importPath = paste0(getwd(), "/neo4j.import"),
    .opts = list(ssl_verifypeer=0)
)

#------------------------------------------#
# Création des fonctions utilitaires
#------------------------------------------#
cq = function(query, neo4j=neodb, ...) cypher(neo4j, query, arraysAsStrings = F, ...)

#------------------------------------------#
# Chargement des données
#------------------------------------------#
# args <- commandArgs(trailingOnly = TRUE)
# if (length(args) < 1) {
#     args[1]="set.01"
# }
# input_file <- paste0("tmp/",args[1],".out")
# cat("Fichier fourni :", input_file, "\n")
# # Lire le fichier
# set1 <- read_tsv(input_file, col_names=F)

go <- get_ontology("downloads/go-basic.obo", extract_tags="everything")
anno <- read_csv("neo4j.import/uniprot.GOTerm.bnumber.csv")
keywords <- read_csv("neo4j.import/uniprot.keywords.genes.csv")

############################################
# set.01
############################################

set1 <- read_tsv("tmp/set.01.out", col_names=F)
s1 <- scan('query.sets/set.01.txt', character()) 
pwd <- "set.01"

#------------------------------------------#
# Calcul des similarités entre les GOTerms
#------------------------------------------#
simBP <- calculateSimMatrix(
    set1$X1,
    orgdb  = "org.Hs.eg.db", # proxy humain
    ont    = "BP",
    method = "Resnik"
)

# simBP
str(simBP)

scores <- setNames(-log10(set1$X2), set1$X1)
# scores <- setNames(-log10(set1$pvalue), set1$id)
scores

# Clustering des GOTerms
# reducedTerms <- reduceSimMatrix(simBP, scores, threshold=.7, orgdb="org.EcK12.eg.db")
reducedTerms <- reduceSimMatrix(
    simMatrix = simBP,
    scores    = scores,
    threshold = 0.7,
    orgdb     = "org.Hs.eg.db"
)

# heatmapPlot(simBP, reducedTerms, annotateParent=T, annotationLabel="parentTerm", fontsize=6)
# heatmapPlot(simBP, reducedTerms, annotateParent=T, annotationLabel="parentTerm", fontsize=12)
savePath <- paste0("generated.data/",pwd,"/heatmap_GO.png")
png(savePath, width = 1800, height = 1600, res = 200)
heatmapPlot(simBP, reducedTerms, annotateParent=TRUE, annotationLabel="parentTerm", fontsize=9)
dev.off()

# Visualisation avec un nuage de points dans l’espace sémantique
savePath <- paste0("generated.data/",pwd,"/scatterPlot_GO.png")
png(savePath, width = 1800, height = 1600, res = 200)
scatterPlot(simBP, reducedTerms)
dev.off()

# Visualisation avec une treemap
savePath <- paste0("generated.data/",pwd,"/treemapPlot_GO.png")
png(savePath, width = 1800, height = 1600, res = 200)
treemapPlot(reducedTerms)
dev.off()

# Visualisation avec un wordcloud
savePath <- paste0("generated.data/",pwd,"/wordcloudPlot_GO.png")
png(savePath, width = 1800, height = 1600, res = 200)
wordcloudPlot(reducedTerms, min.freq=1)
dev.off()


#==========================================#
# Approches par graphe ou matrice de distances/dissemblances/similarités/probabilités
#==========================================#

#------------------------------------------#
# Génération de matrices de similarité
#------------------------------------------#
# Ajout d’une colonne avec des 1 pour les association présentes dans les données de départ
keyword.mat <- keywords %>% 
    unique %>% 
    mutate(asso=1)
keyword.mat

# Génération d’une table pivot (les valeurs de la colonne keyword deviennent des colonnes)
keyword.tab <- keyword.mat %>% 
    pivot_wider(names_from = keyword, values_from = asso, values_fill = 0) %>% 
    as.data.frame
keyword.tab %>% head

# Transfert des identifiants de gènes de la colonne bnumber aux noms de lignes de la matrice de distance
rownames(keyword.tab) <- keyword.tab$bnumber
keyword.tab %>% head

# Calcul de la matrice de distance (entre paires de lignes/gènes) avec la méthode binary
keyword.tab <- keyword.tab %>% dplyr::select(-bnumber)
keyword.dist <- keyword.tab %>% dist(method='binary')
keyword.dist[1:10]

# Sous forme de scores à la StringDB
keyword.scores <- round( (1-keyword.dist) * 1000) %>% as.matrix
keyword.scores[1:10,1:10]

#------------------------------------------#
# Client R à STRINGdb
#------------------------------------------#
kabex <- function(o) o %>% kable(format="html", escape=F) %>% kable_styling(bootstrap_options = c("striped"))
library(DT) 
dted <- function(o) o %>% datatable(rownames=F, filter="top", options=list(pageLength = 10), escape=F)

# Connexion au serveur de STRING (NCBI taxon id de Ecoli est 511145)
# https://string-db.org/cgi/input.pl?input_page_active_form=organism
confidence_threshold <- 333
stringdb <- STRINGdb$new(version='12', species=511145, score_threshold=confidence_threshold, input_directory='downloads')

# Téléchargement du génome/protéome
proteins <- stringdb$get_proteins() %>% as_tibble
dt <- proteins %>% dted()
savePath <- paste0("generated.data/",pwd,"/proteins_table.html")
saveWidget(dt, savePath, selfcontained = TRUE)

# Mapping des identifiants vers ceux de StringDB
s1.mapped <- stringdb$mp(s1)

# Plot du graphe correspondant
stringdb$plot_network(s1.mapped)

# Enrichment
enrichment <- stringdb$get_enrichment(s1.mapped)
dtenrichment <- enrichment %>% dted()
savePath <- paste0("generated.data/",pwd,"/proteins_enrichment_table.html")
saveWidget(dtenrichment, savePath, selfcontained = TRUE)

# Sommets voisins d’un ou plusieurs sommets donnés
stringdb$get_neighbors(s1.mapped[1:3])

# Interactions entre un ensemble de sommets donnés
stringdb$get_interactions(s1.mapped)


#==========================================#
# Données de coexpression à partir des données complètes
#==========================================#

#------------------------------------------#
# STRINGdb detailed links → neo4j
#------------------------------------------#
# Chargement sous forme de tibble
links.detailed <- read_delim("downloads/511145.protein.links.detailed.v12.0.txt.gz", delim=" ", col_types = "ccnnnnnnnn")

# Génération des fichiers CSV pour l’import
links.detailed %>%
    filter(coexpression>0 & protein1 < protein2) %>%
    dplyr::select(protein1, protein2, coexpression) %>%
    mutate(organism=str_extract(protein1, '^\\d+'), id1=str_extract(protein1, 'b\\d+'), id2=str_extract(protein2, 'b\\d+')) %>%
    dplyr::select(organism:id2,coexpression) %>%
    write_csv("neo4j.import/string.coexpression.csv")

# Import dans neo4j
'
LOAD CSV WITH HEADERS FROM "file:///string.coexpression.csv" AS line
MATCH (g1:Gene),(g2:Gene)
 WHERE g1.id=line.id1 AND g2.id=line.id2
WITH g1,g2, toInteger(line.coexpression) AS value 
MERGE (g1)-[r:STRINGdb]-(g2)
 SET r.coexpression=value
'  %>% cq

# Vérification
'MATCH ()-[r:STRINGdb]-() RETURN count(r)' %>% cq %>% unlist

# Extraction des liens de coexpression d’au moins 0.4
coex <- "MATCH (g1:Gene)-[r:STRINGdb]-(g2:Gene) WHERE r.coexpression>=995 RETURN g1.id, g2.id, r.coexpression" %>% cq 
tibble(id1=coex$g1.id, id2=coex$g2.id, coexpression=coex$r.coexpression)

# Extraction d’un sous-graphe : on utilise ici neo4r et la fonction fournie (call_neo4j) avec le paramètre type="graph" afin de récupérer un sous-graphe.
g.coexpr <- 'MATCH p = ()-[r:STRINGdb]-() WHERE r.coexpression>=995 RETURN p'  %>% cq(result="graph")

# Manipulation pour igraph pour récupérer les propriétés des sommets sous forme de tibble
g.coexpr$nodes <-  g.coexpr$nodes |> 
    map(\(node) c(neo_id=node$id, node$properties)) %>% 
    bind_rows()
dtcoexprnodes <- g.coexpr$nodes %>% dted()
savePath <- paste0("generated.data/",pwd,"/proteins_coexpression_nodes_table.html")
saveWidget(dtcoexprnodes, savePath, selfcontained = TRUE)

# De même pour les arcs/arêtes du graphe à l’aide de la fonction unnest_relationships pour le passage dans igraph : besoin de réordonner les liens
g.coexpr$relationships <- g.coexpr$relationships |> 
    map(\(edge) c(startNode=edge$startNode, endNode=edge$endNode, type=edge$type, edge$properties)) %>% 
    bind_rows
dtcoexprrelation <- g.coexpr$relationships %>% dted()
savePath <- paste0("generated.data/",pwd,"/proteins_coexpression_relations_table.html")
saveWidget(dtcoexprrelation, savePath, selfcontained = TRUE)


#==========================================#
# Graphe
#==========================================#
# Chargement de la librairies et positionnement des valeurs par défaut pour certains paramètres
if (!require("igraph", quietly = TRUE)){install.packages("igraph", repos='https://mirror.ibcp.fr/pub/CRAN/')}
library(igraph)

igraph.options(vertex.color=NA)
igraph.options(vertex.label.cex=.6) # font size
igraph.options(vertex.label.family='sans')
igraph.options(vertex.size=2)
igraph.options(edge.label.cex=.6)
igraph.options(edge.label.family='sans')

# Création du graphe dans igraph à partir des sommets et des relations
# remarque : la première colonne du df passé comme vertices est considérée comme l’identifiant des sommes, donc neo_id ici
g.coexpr = graph_from_data_frame(d=g.coexpr$relationships, directed=FALSE, vertices = g.coexpr$nodes)

# Plot
savePath <- paste0("generated.data/",pwd,"/coexpression.png")
png(savePath, width = 1800, height = 1600, res = 200)
plot(g.coexpr, vertex.label=NA, main='coexpression')
dev.off()


#==========================================#
# Autres scores StringDB à intégrer
#==========================================#

#------------------------------------------#
# experimental
#------------------------------------------#
# Génération des fichiers CSV pour l’import
links.detailed %>%
    filter(experimental > 0 & protein1 < protein2) %>%
    dplyr::select(protein1, protein2, experimental) %>%
    mutate(
        organism = stringr::str_extract(protein1, '^\\d+'),
        id1 = stringr::str_extract(protein1, 'b\\d+'),
        id2 = stringr::str_extract(protein2, 'b\\d+')
    ) %>%
    dplyr::select(organism:id2, experimental) %>%
    write_csv("neo4j.import/string.experimental.csv")

# Import dans neo4j
'
LOAD CSV WITH HEADERS FROM "file:///string.experimental.csv" AS line
MATCH (g1:Gene),(g2:Gene)
    WHERE g1.id = line.id1 AND g2.id = line.id2
WITH g1,g2, toInteger(line.experimental) AS value 
MERGE (g1)-[r:STRINGdb_experimental]-(g2)
    SET r.experimental = value
' %>% cq

# Vérification
'MATCH ()-[r:STRINGdb_experimental]-() RETURN count(r)' %>% cq %>% unlist

#------------------------------------------#
# neighborhood
#------------------------------------------#
# Génération des fichiers CSV pour l’import
links.detailed %>%
    filter(neighborhood > 0 & protein1 < protein2) %>%
    dplyr::select(protein1, protein2, neighborhood) %>%
    mutate(
        organism = stringr::str_extract(protein1, '^\\d+'),
        id1 = stringr::str_extract(protein1, 'b\\d+'),
        id2 = stringr::str_extract(protein2, 'b\\d+')
    ) %>%
    dplyr::select(organism:id2, neighborhood) %>%
    write_csv("neo4j.import/string.neighborhood.csv")

# Import dans neo4j
'
LOAD CSV WITH HEADERS FROM "file:///string.neighborhood.csv" AS line
MATCH (g1:Gene),(g2:Gene)
    WHERE g1.id = line.id1 AND g2.id = line.id2
WITH g1,g2, toInteger(line.neighborhood) AS value 
MERGE (g1)-[r:STRINGdb_neighborhood]-(g2)
    SET r.neighborhood = value
' %>% cq

# Vérification
'MATCH ()-[r:STRINGdb_neighborhood]-() RETURN count(r)' %>% cq %>% unlist

#------------------------------------------#
# textmining
#------------------------------------------#
# Génération des fichiers CSV pour l’import
links.detailed %>%
    filter(textmining > 0 & protein1 < protein2) %>%
    dplyr::select(protein1, protein2, textmining) %>%
    mutate(
        organism = stringr::str_extract(protein1, '^\\d+'),
        id1 = stringr::str_extract(protein1, 'b\\d+'),
        id2 = stringr::str_extract(protein2, 'b\\d+')
    ) %>%
    dplyr::select(organism:id2, textmining) %>%
    write_csv("neo4j.import/string.textmining.csv")

# Import dans neo4j
'
LOAD CSV WITH HEADERS FROM "file:///string.textmining.csv" AS line
MATCH (g1:Gene),(g2:Gene)
    WHERE g1.id = line.id1 AND g2.id = line.id2
WITH g1,g2, toInteger(line.textmining) AS value 
MERGE (g1)-[r:STRINGdb_textmining]-(g2)
    SET r.textmining = value
' %>% cq

# Vérification
'MATCH ()-[r:STRINGdb_textmining]-() RETURN count(r)' %>% cq %>% unlist

#------------------------------------------#
# database
#------------------------------------------#
# Génération des fichiers CSV pour l’import
links.detailed %>%
    filter(database > 0 & protein1 < protein2) %>%
    dplyr::select(protein1, protein2, database) %>%
    mutate(
        organism = stringr::str_extract(protein1, '^\\d+'),
        id1 = stringr::str_extract(protein1, 'b\\d+'),
        id2 = stringr::str_extract(protein2, 'b\\d+')
    ) %>%
    dplyr::select(organism:id2, database) %>%
    write_csv("neo4j.import/string.database.csv")

# Import dans neo4j
'
LOAD CSV WITH HEADERS FROM "file:///string.database.csv" AS line
MATCH (g1:Gene),(g2:Gene)
    WHERE g1.id = line.id1 AND g2.id = line.id2
WITH g1,g2, toInteger(line.database) AS value 
MERGE (g1)-[r:STRINGdb_database]-(g2)
    SET r.database = value
' %>% cq

# Vérification
'MATCH ()-[r:STRINGdb_database]-() RETURN count(r)' %>% cq %>% unlist

#------------------------------------------#
# combined_score
#------------------------------------------#
# Génération des fichiers CSV pour l’import
links.detailed %>%
    filter(combined_score > 0 & protein1 < protein2) %>%
    dplyr::select(protein1, protein2, combined_score) %>%
    mutate(
        organism = stringr::str_extract(protein1, '^\\d+'),
        id1 = stringr::str_extract(protein1, 'b\\d+'),
        id2 = stringr::str_extract(protein2, 'b\\d+')
    ) %>%
    dplyr::select(organism:id2, combined_score) %>%
    write_csv("neo4j.import/string.combined_score.csv")

# Import dans neo4j
'
LOAD CSV WITH HEADERS FROM "file:///string.combined_score.csv" AS line
MATCH (g1:Gene),(g2:Gene)
    WHERE g1.id = line.id1 AND g2.id = line.id2
WITH g1,g2, toInteger(line.combined_score) AS value 
MERGE (g1)-[r:STRINGdb_combined]-(g2)
    SET r.combined_score = value
' %>% cq

# Vérification
'MATCH ()-[r:STRINGdb_combined]-() RETURN count(r)' %>% cq %>% unlist


#==========================================#
# Combinaison de matrices de similarité
#==========================================#
# Récupération du graphe g avec les liens de coexpression, neighborhood et experiment à partir de neo4j
# STRINGdb : r.coexpression
# STRINGdb_experimental : r.experimental (noté ppi)
# STRINGdb_neighborhood : r.neighborhood
# STRINGdb_textmining : r.textmining
# STRINGdb_database : r.database
# STRINGdb_combined : r.combined_score

# On fusionne tout ça par paire de gènes (id1,id2)
edges_raw <- "
MATCH (g1:Gene)-[r:STRINGdb]-(g2:Gene)
OPTIONAL MATCH (g1)-[re:STRINGdb_experimental]-(g2)
OPTIONAL MATCH (g1)-[rn:STRINGdb_neighborhood]-(g2)
OPTIONAL MATCH (g1)-[rt:STRINGdb_textmining]-(g2)
OPTIONAL MATCH (g1)-[rd:STRINGdb_database]-(g2)
OPTIONAL MATCH (g1)-[rc:STRINGdb_combined_score]-(g2)
WHERE g1.id < g2.id   // pour ne garder chaque paire qu'une seule fois
RETURN g1.id AS id1,
       g2.id AS id2,
       r.coexpression AS coexpression,
       re.experimental AS ppi,
       rn.neighborhood AS neighborhood,
       rt.textmining AS textmining,
       rd.database AS database,
       rc.combined_score AS combined_score
" %>% cq
# remplissage des NA par 0
edges_raw[is.na(edges_raw)] <- 0

# Mise en forme en tibble
edges_df <- tibble(
    id1 = edges_raw$id1,
    id2 = edges_raw$id2,
    coexpression = as.numeric(edges_raw$coexpression),
    ppi = as.numeric(edges_raw$ppi),
    neighborhood = as.numeric(edges_raw$neighborhood),
    textmining = as.numeric(edges_raw$textmining),
    database = as.numeric(edges_raw$database),
    combined_score = as.numeric(edges_raw$combined_score)
)

# Création du graphe non orienté dans igraph
# Les colonnes id1 / id2 servent d’extrémités des arêtes,
# et les autres colonnes deviennent des attributs d’arêtes.
g <- graph_from_data_frame(
    d = edges_df,
    directed = FALSE,
    vertices = NULL # les sommets sont déduits de id1 / id2
)

# Calcul du score combiné par STRINGdb (source: https://string-db.org/help/faq/#how-are-the-scores-computed)
prior <- 0.041
no_prior <- function(x, prior = 0.041) (ifelse(is.na(x), 0, x) / 1000 - prior) / (1-prior)
s_coexp_nop <- no_prior(E(g)$coexpression)
s_ppi_nop <- no_prior(E(g)$ppi)
s_neighborhood_nop <- no_prior(E(g)$neighborhood)
s_tot_nop <- 1 - (1 - s_coexp_nop) * (1 - s_ppi_nop) * (1 - s_neighborhood_nop)
E(g)$combined_score <- round(1000 * (s_tot_nop + prior *(1 - s_tot_nop)))


#==========================================#
# Prioritisation de gènes
#==========================================#
# Gènes “training”.
# Choix du cycle des citrates (pathway TCA).
# Récupération des gènes impliqués dans le pathway
training.genes <- "MATCH (Pathway {id: 'TCA'})-[:requires]->(g: Gene {organism: 511145}) return g.bnumber " %>% cq %>% .$g %>% unlist

# Gènes “candidats” (ensemble du génome)
candidates <- "MATCH (g:Gene {organism: 511145}) RETURN g.bnumber" %>% cq %>% .$g 

# Distance d’un gène à l’ensemble de référence
score <- function(gene.ids, ref.genes, datasource) {
    # mapping datasource -> type de relation dans Neo4j
    rel_type <- switch(
        datasource,
        coexpression = "STRINGdb",
        experimental = "STRINGdb_experimental",
        neighborhood = "STRINGdb_neighborhood",
        textmining   = "STRINGdb_textmining",
        database     = "STRINGdb_database",
        stop(paste("Datasource inconnue :", datasource))
    )
    ref_str <- paste("'", ref.genes, "'", sep = "", collapse = ",")
    sapply(gene.ids, function(gene.id) {
        query <- paste0(
        "MATCH (g1:Gene {id: '", gene.id, "', organism: 511145})",
        "-[r:", rel_type, "]-",
        "(g2:Gene {organism: 511145}) ",
        "WHERE g2.id IN [", ref_str, "] ",
        "RETURN r.", datasource, " AS score"
        )
        res <- query %>% cq %>% .$score
        ifelse(is.null(res), 0, mean(replace_na(res, 0)))
    })
}

# test
gene.ids <- 'b0001'
gene.ids <- 'b0116'
ref.genes <- training.genes
datasource <- 'coexpression'
score(gene.ids, training.genes, datasource)

# Application à l’ensemble des gènes
scores <- tibble(candidate=candidates) %>%
  mutate(
    coexpression = score(gene.id=candidate, ref.genes=training.genes, datasource='coexpression'),
    experimental = score(gene.id=candidate, ref.genes=training.genes, datasource='experimental'),
    neighborhood = score(gene.id=candidate, ref.genes=training.genes, datasource='neighborhood'),
    textmining = score(gene.id=candidate, ref.genes=training.genes, datasource='textmining'),
    database = score(gene.id=candidate, ref.genes=training.genes, datasource='database')
    )
scores <- scores %>% replace(is.na(.), 0)

# Plots
savePath <- paste0("generated.data/",pwd,"/priorisation_plot.png")
png(savePath, width = 1800, height = 1600, res = 200)
scores %>%
    ggplot(aes(x=coexpression, y=neighborhood)) +
    geom_point(alpha=0.2, aes(color=candidate %in% training.genes, shape=candidate %in% training.genes))
dev.off()

# ACP pour l’affichage
pca = scores %>%
    dplyr::select(-candidate) %>%
    prcomp
pca %>% summary

# Visualisation
library(factoextra)
savePath <- paste0("generated.data/",pwd,"/priorisation_plot-pca.png")
png(savePath, width = 1800, height = 1600, res = 200)
pca %>% fviz_pca_ind(col.ind = scores$candidate %in% training.genes, addEllipses = T, alpha.ind = .2, geom='point')
dev.off()

# Contribution des sources de données
savePath <- paste0("generated.data/",pwd,"/priorisation_plot-pca-contribution.png")
png(savePath, width = 1800, height = 1600, res = 200)
pca %>% fviz_pca_biplot(col.ind = scores$candidate %in% training.genes)
dev.off()

# Analyse discriminante linéaire
if (!require("MASS", quietly = TRUE)){install.packages("MASS", repos='https://mirror.ibcp.fr/pub/CRAN/')}
library(MASS)
mat <- scores %>% dplyr::select(-candidate) %>% as.matrix
model <- lda(x=mat, grouping = scores$candidate %in% training.genes)

# Distribution des scores sous forme de boîtes à moustache pour les 2 classes
savePath <- paste0("generated.data/",pwd,"/priorisation_plot-LDA-boxplot.png")
png(savePath, width = 1800, height = 1600, res = 200)
mat %>% 
    scale( center=T, scale=F ) %*% model$scaling %>% 
    as_tibble %>%
    mutate(gene.id=scores$candidate) %>%
    ggplot(aes(x=LD1, y=gene.id %in% training.genes, color=gene.id %in% training.genes, shape=gene.id %in% training.genes)) +
    geom_violin() +
    geom_boxplot(varwidth = T) +
    geom_jitter(height = 0.2, alpha=0.1, color='grey') +
    theme_light()
dev.off()

# Density
savePath <- paste0("generated.data/",pwd,"/priorisation_plot-density.png")
png(savePath, width = 1800, height = 1600, res = 200)
mat %>% 
    scale( center=T, scale=F ) %*% model$scaling %>% 
    as_tibble %>%
    mutate(gene.id=scores$candidate) %>%
    ggplot(aes(LD1)) + #, color=gene.id %in% training.genes, shape=gene.id %in% training.genes)) +
    geom_density()
dev.off()

#------------------------------------------#
# Tentative d’optimisation
#------------------------------------------#

# Récupération des “scores” entre les candidats et les gènes de référence en une requête
train.str <- paste("'", training.genes,"'", sep = '', collapse = ',')
datasources <- c('coexpression', 'experimental', 'neighborhood', 'textmining', 'database' )
ds.str <- paste(" r.", datasources, sep='', collapse=',')

tb <- paste0("
MATCH (candidate:Gene {organism: 511145}), (training:Gene {organism: 511145})
WHERE training.id IN [", train.str, "] AND candidate.id <> training.id

OPTIONAL MATCH (candidate)-[rc:STRINGdb]-(training)
OPTIONAL MATCH (candidate)-[re:STRINGdb_experimental]-(training)
OPTIONAL MATCH (candidate)-[rn:STRINGdb_neighborhood]-(training)
OPTIONAL MATCH (candidate)-[rt:STRINGdb_textmining]-(training)
OPTIONAL MATCH (candidate)-[rd:STRINGdb_database]-(training)

RETURN DISTINCT
    candidate.id AS candidate_id,
    training.id  AS training_id,
    rc.coexpression  AS coexpression,
    re.experimental  AS experimental,
    rn.neighborhood  AS neighborhood,
    rt.textmining    AS textmining,
    rd.database      AS database
") %>% cq

# Replace NAs with 0s and compute mean (or other) by candidate
scores <- tb %>%
    dplyr::select(-training_id) %>%
    replace(is.na(.), 0) %>%
    group_by(candidate_id) %>%
    summarise(across(all_of(datasources), mean), .groups = "drop")

# Suppression des colonnes entièrement nulles (aucune info)
scores <- scores %>%
    dplyr::select(candidate_id, where(~ any(. != 0)))

# Analyse discriminante linéaire
mat <- scores %>% dplyr::select(-candidate_id) %>% as.matrix
model <- lda(x=mat, grouping = scores$candidate_id %in% training.genes)

# Distribution des scores sous forme de boîtes à moustache pour les 2 classes
savePath <- paste0("generated.data/",pwd,"/priorisation_plot-LDA-boxplot-optimisation.png")
png(savePath, width = 1800, height = 1600, res = 200)
mat %>% 
    scale( center=T, scale=F ) %*% model$scaling %>% 
    as_tibble %>%
    mutate(gene.id=scores$candidate_id, training=gene.id %in% training.genes) %>%
    ggplot(aes(x=LD1, y=training, color=training, shape=training)) +
    geom_violin() +
    geom_jitter(alpha=.5, height = .2)
dev.off()


############################################
# set.02
############################################

set1 <- read_tsv("tmp/set.02.out", col_names=F)
s1 <- scan('query.sets/set.02.txt', character()) 
pwd <- "set.02"

#------------------------------------------#
# Calcul des similarités entre les GOTerms
#------------------------------------------#
simBP <- calculateSimMatrix(
    set1$X1,
    orgdb  = "org.Hs.eg.db", # proxy humain
    ont    = "BP",
    method = "Resnik"
)

# simBP
str(simBP)

scores <- setNames(-log10(set1$X2), set1$X1)
# scores <- setNames(-log10(set1$pvalue), set1$id)
scores

# Clustering des GOTerms
# reducedTerms <- reduceSimMatrix(simBP, scores, threshold=.7, orgdb="org.EcK12.eg.db")
reducedTerms <- reduceSimMatrix(
    simMatrix = simBP,
    scores    = scores,
    threshold = 0.7,
    orgdb     = "org.Hs.eg.db"
)

# heatmapPlot(simBP, reducedTerms, annotateParent=T, annotationLabel="parentTerm", fontsize=6)
# heatmapPlot(simBP, reducedTerms, annotateParent=T, annotationLabel="parentTerm", fontsize=12)
savePath <- paste0("generated.data/",pwd,"/heatmap_GO.png")
png(savePath, width = 1800, height = 1600, res = 200)
heatmapPlot(simBP, reducedTerms, annotateParent=TRUE, annotationLabel="parentTerm", fontsize=9)
dev.off()

# Visualisation avec un nuage de points dans l’espace sémantique
savePath <- paste0("generated.data/",pwd,"/scatterPlot_GO.png")
png(savePath, width = 1800, height = 1600, res = 200)
scatterPlot(simBP, reducedTerms)
dev.off()

# Visualisation avec une treemap
savePath <- paste0("generated.data/",pwd,"/treemapPlot_GO.png")
png(savePath, width = 1800, height = 1600, res = 200)
treemapPlot(reducedTerms)
dev.off()

# Visualisation avec un wordcloud
savePath <- paste0("generated.data/",pwd,"/wordcloudPlot_GO.png")
png(savePath, width = 1800, height = 1600, res = 200)
wordcloudPlot(reducedTerms, min.freq=1)
dev.off()


#==========================================#
# Approches par graphe ou matrice de distances/dissemblances/similarités/probabilités
#==========================================#

#------------------------------------------#
# Génération de matrices de similarité
#------------------------------------------#
# Ajout d’une colonne avec des 1 pour les association présentes dans les données de départ
keyword.mat <- keywords %>% 
    unique %>% 
    mutate(asso=1)
keyword.mat

# Génération d’une table pivot (les valeurs de la colonne keyword deviennent des colonnes)
keyword.tab <- keyword.mat %>% 
    pivot_wider(names_from = keyword, values_from = asso, values_fill = 0) %>% 
    as.data.frame
keyword.tab %>% head

# Transfert des identifiants de gènes de la colonne bnumber aux noms de lignes de la matrice de distance
rownames(keyword.tab) <- keyword.tab$bnumber
keyword.tab %>% head

# Calcul de la matrice de distance (entre paires de lignes/gènes) avec la méthode binary
keyword.tab <- keyword.tab %>% dplyr::select(-bnumber)
keyword.dist <- keyword.tab %>% dist(method='binary')
keyword.dist[1:10]

# Sous forme de scores à la StringDB
keyword.scores <- round( (1-keyword.dist) * 1000) %>% as.matrix
keyword.scores[1:10,1:10]

#------------------------------------------#
# Client R à STRINGdb
#------------------------------------------#
kabex <- function(o) o %>% kable(format="html", escape=F) %>% kable_styling(bootstrap_options = c("striped"))
library(DT) 
dted <- function(o) o %>% datatable(rownames=F, filter="top", options=list(pageLength = 10), escape=F)

# Connexion au serveur de STRING (NCBI taxon id de Ecoli est 511145)
# https://string-db.org/cgi/input.pl?input_page_active_form=organism
confidence_threshold <- 333
stringdb <- STRINGdb$new(version='12', species=511145, score_threshold=confidence_threshold, input_directory='downloads')

# Téléchargement du génome/protéome
proteins <- stringdb$get_proteins() %>% as_tibble
dt <- proteins %>% dted()
savePath <- paste0("generated.data/",pwd,"/proteins_table.html")
saveWidget(dt, savePath, selfcontained = TRUE)

# Mapping des identifiants vers ceux de StringDB
s1.mapped <- stringdb$mp(s1)

# Plot du graphe correspondant
stringdb$plot_network(s1.mapped)

# Enrichment
enrichment <- stringdb$get_enrichment(s1.mapped)
dtenrichment <- enrichment %>% dted()
savePath <- paste0("generated.data/",pwd,"/proteins_enrichment_table.html")
saveWidget(dtenrichment, savePath, selfcontained = TRUE)

# Sommets voisins d’un ou plusieurs sommets donnés
stringdb$get_neighbors(s1.mapped[1:3])

# Interactions entre un ensemble de sommets donnés
stringdb$get_interactions(s1.mapped)


#==========================================#
# Données de coexpression à partir des données complètes
#==========================================#

#------------------------------------------#
# STRINGdb detailed links → neo4j
#------------------------------------------#
# Chargement sous forme de tibble
links.detailed <- read_delim("downloads/511145.protein.links.detailed.v12.0.txt.gz", delim=" ", col_types = "ccnnnnnnnn")

# Génération des fichiers CSV pour l’import
links.detailed %>%
    filter(coexpression>0 & protein1 < protein2) %>%
    dplyr::select(protein1, protein2, coexpression) %>%
    mutate(organism=str_extract(protein1, '^\\d+'), id1=str_extract(protein1, 'b\\d+'), id2=str_extract(protein2, 'b\\d+')) %>%
    dplyr::select(organism:id2,coexpression) %>%
    write_csv("neo4j.import/string.coexpression.csv")

# Import dans neo4j
'
LOAD CSV WITH HEADERS FROM "file:///string.coexpression.csv" AS line
MATCH (g1:Gene),(g2:Gene)
 WHERE g1.id=line.id1 AND g2.id=line.id2
WITH g1,g2, toInteger(line.coexpression) AS value 
MERGE (g1)-[r:STRINGdb]-(g2)
 SET r.coexpression=value
'  %>% cq

# Vérification
'MATCH ()-[r:STRINGdb]-() RETURN count(r)' %>% cq %>% unlist

# Extraction des liens de coexpression d’au moins 0.4
coex <- "MATCH (g1:Gene)-[r:STRINGdb]-(g2:Gene) WHERE r.coexpression>=995 RETURN g1.id, g2.id, r.coexpression" %>% cq 
tibble(id1=coex$g1.id, id2=coex$g2.id, coexpression=coex$r.coexpression)

# Extraction d’un sous-graphe : on utilise ici neo4r et la fonction fournie (call_neo4j) avec le paramètre type="graph" afin de récupérer un sous-graphe.
g.coexpr <- 'MATCH p = ()-[r:STRINGdb]-() WHERE r.coexpression>=995 RETURN p'  %>% cq(result="graph")

# Manipulation pour igraph pour récupérer les propriétés des sommets sous forme de tibble
g.coexpr$nodes <-  g.coexpr$nodes |> 
    map(\(node) c(neo_id=node$id, node$properties)) %>% 
    bind_rows()
dtcoexprnodes <- g.coexpr$nodes %>% dted()
savePath <- paste0("generated.data/",pwd,"/proteins_coexpression_nodes_table.html")
saveWidget(dtcoexprnodes, savePath, selfcontained = TRUE)

# De même pour les arcs/arêtes du graphe à l’aide de la fonction unnest_relationships pour le passage dans igraph : besoin de réordonner les liens
g.coexpr$relationships <- g.coexpr$relationships |> 
    map(\(edge) c(startNode=edge$startNode, endNode=edge$endNode, type=edge$type, edge$properties)) %>% 
    bind_rows
dtcoexprrelation <- g.coexpr$relationships %>% dted()
savePath <- paste0("generated.data/",pwd,"/proteins_coexpression_relations_table.html")
saveWidget(dtcoexprrelation, savePath, selfcontained = TRUE)


#==========================================#
# Graphe
#==========================================#
# Chargement de la librairies et positionnement des valeurs par défaut pour certains paramètres
if (!require("igraph", quietly = TRUE)){install.packages("igraph", repos='https://mirror.ibcp.fr/pub/CRAN/')}
library(igraph)

igraph.options(vertex.color=NA)
igraph.options(vertex.label.cex=.6) # font size
igraph.options(vertex.label.family='sans')
igraph.options(vertex.size=2)
igraph.options(edge.label.cex=.6)
igraph.options(edge.label.family='sans')

# Création du graphe dans igraph à partir des sommets et des relations
# remarque : la première colonne du df passé comme vertices est considérée comme l’identifiant des sommes, donc neo_id ici
g.coexpr = graph_from_data_frame(d=g.coexpr$relationships, directed=FALSE, vertices = g.coexpr$nodes)

# Plot
savePath <- paste0("generated.data/",pwd,"/coexpression.png")
png(savePath, width = 1800, height = 1600, res = 200)
plot(g.coexpr, vertex.label=NA, main='coexpression')
dev.off()


#==========================================#
# Autres scores StringDB à intégrer
#==========================================#

#------------------------------------------#
# experimental
#------------------------------------------#
# Génération des fichiers CSV pour l’import
links.detailed %>%
    filter(experimental > 0 & protein1 < protein2) %>%
    dplyr::select(protein1, protein2, experimental) %>%
    mutate(
        organism = stringr::str_extract(protein1, '^\\d+'),
        id1 = stringr::str_extract(protein1, 'b\\d+'),
        id2 = stringr::str_extract(protein2, 'b\\d+')
    ) %>%
    dplyr::select(organism:id2, experimental) %>%
    write_csv("neo4j.import/string.experimental.csv")

# Import dans neo4j
'
LOAD CSV WITH HEADERS FROM "file:///string.experimental.csv" AS line
MATCH (g1:Gene),(g2:Gene)
    WHERE g1.id = line.id1 AND g2.id = line.id2
WITH g1,g2, toInteger(line.experimental) AS value 
MERGE (g1)-[r:STRINGdb_experimental]-(g2)
    SET r.experimental = value
' %>% cq

# Vérification
'MATCH ()-[r:STRINGdb_experimental]-() RETURN count(r)' %>% cq %>% unlist

#------------------------------------------#
# neighborhood
#------------------------------------------#
# Génération des fichiers CSV pour l’import
links.detailed %>%
    filter(neighborhood > 0 & protein1 < protein2) %>%
    dplyr::select(protein1, protein2, neighborhood) %>%
    mutate(
        organism = stringr::str_extract(protein1, '^\\d+'),
        id1 = stringr::str_extract(protein1, 'b\\d+'),
        id2 = stringr::str_extract(protein2, 'b\\d+')
    ) %>%
    dplyr::select(organism:id2, neighborhood) %>%
    write_csv("neo4j.import/string.neighborhood.csv")

# Import dans neo4j
'
LOAD CSV WITH HEADERS FROM "file:///string.neighborhood.csv" AS line
MATCH (g1:Gene),(g2:Gene)
    WHERE g1.id = line.id1 AND g2.id = line.id2
WITH g1,g2, toInteger(line.neighborhood) AS value 
MERGE (g1)-[r:STRINGdb_neighborhood]-(g2)
    SET r.neighborhood = value
' %>% cq

# Vérification
'MATCH ()-[r:STRINGdb_neighborhood]-() RETURN count(r)' %>% cq %>% unlist

#------------------------------------------#
# textmining
#------------------------------------------#
# Génération des fichiers CSV pour l’import
links.detailed %>%
    filter(textmining > 0 & protein1 < protein2) %>%
    dplyr::select(protein1, protein2, textmining) %>%
    mutate(
        organism = stringr::str_extract(protein1, '^\\d+'),
        id1 = stringr::str_extract(protein1, 'b\\d+'),
        id2 = stringr::str_extract(protein2, 'b\\d+')
    ) %>%
    dplyr::select(organism:id2, textmining) %>%
    write_csv("neo4j.import/string.textmining.csv")

# Import dans neo4j
'
LOAD CSV WITH HEADERS FROM "file:///string.textmining.csv" AS line
MATCH (g1:Gene),(g2:Gene)
    WHERE g1.id = line.id1 AND g2.id = line.id2
WITH g1,g2, toInteger(line.textmining) AS value 
MERGE (g1)-[r:STRINGdb_textmining]-(g2)
    SET r.textmining = value
' %>% cq

# Vérification
'MATCH ()-[r:STRINGdb_textmining]-() RETURN count(r)' %>% cq %>% unlist

#------------------------------------------#
# database
#------------------------------------------#
# Génération des fichiers CSV pour l’import
links.detailed %>%
    filter(database > 0 & protein1 < protein2) %>%
    dplyr::select(protein1, protein2, database) %>%
    mutate(
        organism = stringr::str_extract(protein1, '^\\d+'),
        id1 = stringr::str_extract(protein1, 'b\\d+'),
        id2 = stringr::str_extract(protein2, 'b\\d+')
    ) %>%
    dplyr::select(organism:id2, database) %>%
    write_csv("neo4j.import/string.database.csv")

# Import dans neo4j
'
LOAD CSV WITH HEADERS FROM "file:///string.database.csv" AS line
MATCH (g1:Gene),(g2:Gene)
    WHERE g1.id = line.id1 AND g2.id = line.id2
WITH g1,g2, toInteger(line.database) AS value 
MERGE (g1)-[r:STRINGdb_database]-(g2)
    SET r.database = value
' %>% cq

# Vérification
'MATCH ()-[r:STRINGdb_database]-() RETURN count(r)' %>% cq %>% unlist

#------------------------------------------#
# combined_score
#------------------------------------------#
# Génération des fichiers CSV pour l’import
links.detailed %>%
    filter(combined_score > 0 & protein1 < protein2) %>%
    dplyr::select(protein1, protein2, combined_score) %>%
    mutate(
        organism = stringr::str_extract(protein1, '^\\d+'),
        id1 = stringr::str_extract(protein1, 'b\\d+'),
        id2 = stringr::str_extract(protein2, 'b\\d+')
    ) %>%
    dplyr::select(organism:id2, combined_score) %>%
    write_csv("neo4j.import/string.combined_score.csv")

# Import dans neo4j
'
LOAD CSV WITH HEADERS FROM "file:///string.combined_score.csv" AS line
MATCH (g1:Gene),(g2:Gene)
    WHERE g1.id = line.id1 AND g2.id = line.id2
WITH g1,g2, toInteger(line.combined_score) AS value 
MERGE (g1)-[r:STRINGdb_combined]-(g2)
    SET r.combined_score = value
' %>% cq

# Vérification
'MATCH ()-[r:STRINGdb_combined]-() RETURN count(r)' %>% cq %>% unlist


#==========================================#
# Combinaison de matrices de similarité
#==========================================#
# Récupération du graphe g avec les liens de coexpression, neighborhood et experiment à partir de neo4j
# STRINGdb : r.coexpression
# STRINGdb_experimental : r.experimental (noté ppi)
# STRINGdb_neighborhood : r.neighborhood
# STRINGdb_textmining : r.textmining
# STRINGdb_database : r.database
# STRINGdb_combined : r.combined_score

# On fusionne tout ça par paire de gènes (id1,id2)
edges_raw <- "
MATCH (g1:Gene)-[r:STRINGdb]-(g2:Gene)
OPTIONAL MATCH (g1)-[re:STRINGdb_experimental]-(g2)
OPTIONAL MATCH (g1)-[rn:STRINGdb_neighborhood]-(g2)
OPTIONAL MATCH (g1)-[rt:STRINGdb_textmining]-(g2)
OPTIONAL MATCH (g1)-[rd:STRINGdb_database]-(g2)
OPTIONAL MATCH (g1)-[rc:STRINGdb_combined_score]-(g2)
WHERE g1.id < g2.id   // pour ne garder chaque paire qu'une seule fois
RETURN g1.id AS id1,
       g2.id AS id2,
       r.coexpression AS coexpression,
       re.experimental AS ppi,
       rn.neighborhood AS neighborhood,
       rt.textmining AS textmining,
       rd.database AS database,
       rc.combined_score AS combined_score
" %>% cq
# remplissage des NA par 0
edges_raw[is.na(edges_raw)] <- 0

# Mise en forme en tibble
edges_df <- tibble(
    id1 = edges_raw$id1,
    id2 = edges_raw$id2,
    coexpression = as.numeric(edges_raw$coexpression),
    ppi = as.numeric(edges_raw$ppi),
    neighborhood = as.numeric(edges_raw$neighborhood),
    textmining = as.numeric(edges_raw$textmining),
    database = as.numeric(edges_raw$database),
    combined_score = as.numeric(edges_raw$combined_score)
)

# Création du graphe non orienté dans igraph
# Les colonnes id1 / id2 servent d’extrémités des arêtes,
# et les autres colonnes deviennent des attributs d’arêtes.
g <- graph_from_data_frame(
    d = edges_df,
    directed = FALSE,
    vertices = NULL # les sommets sont déduits de id1 / id2
)

# Calcul du score combiné par STRINGdb (source: https://string-db.org/help/faq/#how-are-the-scores-computed)
prior <- 0.041
no_prior <- function(x, prior = 0.041) (ifelse(is.na(x), 0, x) / 1000 - prior) / (1-prior)
s_coexp_nop <- no_prior(E(g)$coexpression)
s_ppi_nop <- no_prior(E(g)$ppi)
s_neighborhood_nop <- no_prior(E(g)$neighborhood)
s_tot_nop <- 1 - (1 - s_coexp_nop) * (1 - s_ppi_nop) * (1 - s_neighborhood_nop)
E(g)$combined_score <- round(1000 * (s_tot_nop + prior *(1 - s_tot_nop)))


#==========================================#
# Prioritisation de gènes
#==========================================#
# Gènes “training”.
# Choix du cycle des citrates (pathway TCA).
# Récupération des gènes impliqués dans le pathway
training.genes <- "MATCH (Pathway {id: 'TCA'})-[:requires]->(g: Gene {organism: 511145}) return g.bnumber " %>% cq %>% .$g %>% unlist

# Gènes “candidats” (ensemble du génome)
candidates <- "MATCH (g:Gene {organism: 511145}) RETURN g.bnumber" %>% cq %>% .$g 

# Distance d’un gène à l’ensemble de référence
score <- function(gene.ids, ref.genes, datasource) {
    # mapping datasource -> type de relation dans Neo4j
    rel_type <- switch(
        datasource,
        coexpression = "STRINGdb",
        experimental = "STRINGdb_experimental",
        neighborhood = "STRINGdb_neighborhood",
        textmining   = "STRINGdb_textmining",
        database     = "STRINGdb_database",
        stop(paste("Datasource inconnue :", datasource))
    )
    ref_str <- paste("'", ref.genes, "'", sep = "", collapse = ",")
    sapply(gene.ids, function(gene.id) {
        query <- paste0(
        "MATCH (g1:Gene {id: '", gene.id, "', organism: 511145})",
        "-[r:", rel_type, "]-",
        "(g2:Gene {organism: 511145}) ",
        "WHERE g2.id IN [", ref_str, "] ",
        "RETURN r.", datasource, " AS score"
        )
        res <- query %>% cq %>% .$score
        ifelse(is.null(res), 0, mean(replace_na(res, 0)))
    })
}

# test
gene.ids <- 'b0001'
gene.ids <- 'b0116'
ref.genes <- training.genes
datasource <- 'coexpression'
score(gene.ids, training.genes, datasource)

# Application à l’ensemble des gènes
scores <- tibble(candidate=candidates) %>%
  mutate(
    coexpression = score(gene.id=candidate, ref.genes=training.genes, datasource='coexpression'),
    experimental = score(gene.id=candidate, ref.genes=training.genes, datasource='experimental'),
    neighborhood = score(gene.id=candidate, ref.genes=training.genes, datasource='neighborhood'),
    textmining = score(gene.id=candidate, ref.genes=training.genes, datasource='textmining'),
    database = score(gene.id=candidate, ref.genes=training.genes, datasource='database')
    )
scores <- scores %>% replace(is.na(.), 0)

# Plots
savePath <- paste0("generated.data/",pwd,"/priorisation_plot.png")
png(savePath, width = 1800, height = 1600, res = 200)
scores %>%
    ggplot(aes(x=coexpression, y=neighborhood)) +
    geom_point(alpha=0.2, aes(color=candidate %in% training.genes, shape=candidate %in% training.genes))
dev.off()

# ACP pour l’affichage
pca = scores %>%
    dplyr::select(-candidate) %>%
    prcomp
pca %>% summary

# Visualisation
library(factoextra)
savePath <- paste0("generated.data/",pwd,"/priorisation_plot-pca.png")
png(savePath, width = 1800, height = 1600, res = 200)
pca %>% fviz_pca_ind(col.ind = scores$candidate %in% training.genes, addEllipses = T, alpha.ind = .2, geom='point')
dev.off()

# Contribution des sources de données
savePath <- paste0("generated.data/",pwd,"/priorisation_plot-pca-contribution.png")
png(savePath, width = 1800, height = 1600, res = 200)
pca %>% fviz_pca_biplot(col.ind = scores$candidate %in% training.genes)
dev.off()

# Analyse discriminante linéaire
if (!require("MASS", quietly = TRUE)){install.packages("MASS", repos='https://mirror.ibcp.fr/pub/CRAN/')}
library(MASS)
mat <- scores %>% dplyr::select(-candidate) %>% as.matrix
model <- lda(x=mat, grouping = scores$candidate %in% training.genes)

# Distribution des scores sous forme de boîtes à moustache pour les 2 classes
savePath <- paste0("generated.data/",pwd,"/priorisation_plot-LDA-boxplot.png")
png(savePath, width = 1800, height = 1600, res = 200)
mat %>% 
    scale( center=T, scale=F ) %*% model$scaling %>% 
    as_tibble %>%
    mutate(gene.id=scores$candidate) %>%
    ggplot(aes(x=LD1, y=gene.id %in% training.genes, color=gene.id %in% training.genes, shape=gene.id %in% training.genes)) +
    geom_violin() +
    geom_boxplot(varwidth = T) +
    geom_jitter(height = 0.2, alpha=0.1, color='grey') +
    theme_light()
dev.off()

# Density
savePath <- paste0("generated.data/",pwd,"/priorisation_plot-density.png")
png(savePath, width = 1800, height = 1600, res = 200)
mat %>% 
    scale( center=T, scale=F ) %*% model$scaling %>% 
    as_tibble %>%
    mutate(gene.id=scores$candidate) %>%
    ggplot(aes(LD1)) + #, color=gene.id %in% training.genes, shape=gene.id %in% training.genes)) +
    geom_density()
dev.off()

#------------------------------------------#
# Tentative d’optimisation
#------------------------------------------#

# Récupération des “scores” entre les candidats et les gènes de référence en une requête
train.str <- paste("'", training.genes,"'", sep = '', collapse = ',')
datasources <- c('coexpression', 'experimental', 'neighborhood', 'textmining', 'database' )
ds.str <- paste(" r.", datasources, sep='', collapse=',')

tb <- paste0("
MATCH (candidate:Gene {organism: 511145}), (training:Gene {organism: 511145})
WHERE training.id IN [", train.str, "] AND candidate.id <> training.id

OPTIONAL MATCH (candidate)-[rc:STRINGdb]-(training)
OPTIONAL MATCH (candidate)-[re:STRINGdb_experimental]-(training)
OPTIONAL MATCH (candidate)-[rn:STRINGdb_neighborhood]-(training)
OPTIONAL MATCH (candidate)-[rt:STRINGdb_textmining]-(training)
OPTIONAL MATCH (candidate)-[rd:STRINGdb_database]-(training)

RETURN DISTINCT
    candidate.id AS candidate_id,
    training.id  AS training_id,
    rc.coexpression  AS coexpression,
    re.experimental  AS experimental,
    rn.neighborhood  AS neighborhood,
    rt.textmining    AS textmining,
    rd.database      AS database
") %>% cq

# Replace NAs with 0s and compute mean (or other) by candidate
scores <- tb %>%
    dplyr::select(-training_id) %>%
    replace(is.na(.), 0) %>%
    group_by(candidate_id) %>%
    summarise(across(all_of(datasources), mean), .groups = "drop")

# Suppression des colonnes entièrement nulles (aucune info)
scores <- scores %>%
    dplyr::select(candidate_id, where(~ any(. != 0)))

# Analyse discriminante linéaire
mat <- scores %>% dplyr::select(-candidate_id) %>% as.matrix
model <- lda(x=mat, grouping = scores$candidate_id %in% training.genes)

# Distribution des scores sous forme de boîtes à moustache pour les 2 classes
savePath <- paste0("generated.data/",pwd,"/priorisation_plot-LDA-boxplot-optimisation.png")
png(savePath, width = 1800, height = 1600, res = 200)
mat %>% 
    scale( center=T, scale=F ) %*% model$scaling %>% 
    as_tibble %>%
    mutate(gene.id=scores$candidate_id, training=gene.id %in% training.genes) %>%
    ggplot(aes(x=LD1, y=training, color=training, shape=training)) +
    geom_violin() +
    geom_jitter(alpha=.5, height = .2)
dev.off()



############################################
# set.03
############################################

set1 <- read_tsv("tmp/set.03.out", col_names=F)
s1 <- scan('query.sets/set.03.txt', character()) 
pwd <- "set.03"

#------------------------------------------#
# Calcul des similarités entre les GOTerms
#------------------------------------------#
simBP <- calculateSimMatrix(
    set1$X1,
    orgdb  = "org.Hs.eg.db", # proxy humain
    ont    = "BP",
    method = "Resnik"
)

# simBP
str(simBP)

scores <- setNames(-log10(set1$X2), set1$X1)
# scores <- setNames(-log10(set1$pvalue), set1$id)
scores

# Clustering des GOTerms
# reducedTerms <- reduceSimMatrix(simBP, scores, threshold=.7, orgdb="org.EcK12.eg.db")
reducedTerms <- reduceSimMatrix(
    simMatrix = simBP,
    scores    = scores,
    threshold = 0.7,
    orgdb     = "org.Hs.eg.db"
)

# heatmapPlot(simBP, reducedTerms, annotateParent=T, annotationLabel="parentTerm", fontsize=6)
# heatmapPlot(simBP, reducedTerms, annotateParent=T, annotationLabel="parentTerm", fontsize=12)
savePath <- paste0("generated.data/",pwd,"/heatmap_GO.png")
png(savePath, width = 1800, height = 1600, res = 200)
heatmapPlot(simBP, reducedTerms, annotateParent=TRUE, annotationLabel="parentTerm", fontsize=9)
dev.off()

# Visualisation avec un nuage de points dans l’espace sémantique
savePath <- paste0("generated.data/",pwd,"/scatterPlot_GO.png")
png(savePath, width = 1800, height = 1600, res = 200)
scatterPlot(simBP, reducedTerms)
dev.off()

# Visualisation avec une treemap
savePath <- paste0("generated.data/",pwd,"/treemapPlot_GO.png")
png(savePath, width = 1800, height = 1600, res = 200)
treemapPlot(reducedTerms)
dev.off()

# Visualisation avec un wordcloud
savePath <- paste0("generated.data/",pwd,"/wordcloudPlot_GO.png")
png(savePath, width = 1800, height = 1600, res = 200)
wordcloudPlot(reducedTerms, min.freq=1)
dev.off()


#==========================================#
# Approches par graphe ou matrice de distances/dissemblances/similarités/probabilités
#==========================================#

#------------------------------------------#
# Génération de matrices de similarité
#------------------------------------------#
# Ajout d’une colonne avec des 1 pour les association présentes dans les données de départ
keyword.mat <- keywords %>% 
    unique %>% 
    mutate(asso=1)
keyword.mat

# Génération d’une table pivot (les valeurs de la colonne keyword deviennent des colonnes)
keyword.tab <- keyword.mat %>% 
    pivot_wider(names_from = keyword, values_from = asso, values_fill = 0) %>% 
    as.data.frame
keyword.tab %>% head

# Transfert des identifiants de gènes de la colonne bnumber aux noms de lignes de la matrice de distance
rownames(keyword.tab) <- keyword.tab$bnumber
keyword.tab %>% head

# Calcul de la matrice de distance (entre paires de lignes/gènes) avec la méthode binary
keyword.tab <- keyword.tab %>% dplyr::select(-bnumber)
keyword.dist <- keyword.tab %>% dist(method='binary')
keyword.dist[1:10]

# Sous forme de scores à la StringDB
keyword.scores <- round( (1-keyword.dist) * 1000) %>% as.matrix
keyword.scores[1:10,1:10]

#------------------------------------------#
# Client R à STRINGdb
#------------------------------------------#
kabex <- function(o) o %>% kable(format="html", escape=F) %>% kable_styling(bootstrap_options = c("striped"))
library(DT) 
dted <- function(o) o %>% datatable(rownames=F, filter="top", options=list(pageLength = 10), escape=F)

# Connexion au serveur de STRING (NCBI taxon id de Ecoli est 511145)
# https://string-db.org/cgi/input.pl?input_page_active_form=organism
confidence_threshold <- 333
stringdb <- STRINGdb$new(version='12', species=511145, score_threshold=confidence_threshold, input_directory='downloads')

# Téléchargement du génome/protéome
proteins <- stringdb$get_proteins() %>% as_tibble
dt <- proteins %>% dted()
savePath <- paste0("generated.data/",pwd,"/proteins_table.html")
saveWidget(dt, savePath, selfcontained = TRUE)

# Mapping des identifiants vers ceux de StringDB
s1.mapped <- stringdb$mp(s1)

# Plot du graphe correspondant
stringdb$plot_network(s1.mapped)

# Enrichment
enrichment <- stringdb$get_enrichment(s1.mapped)
dtenrichment <- enrichment %>% dted()
savePath <- paste0("generated.data/",pwd,"/proteins_enrichment_table.html")
saveWidget(dtenrichment, savePath, selfcontained = TRUE)

# Sommets voisins d’un ou plusieurs sommets donnés
stringdb$get_neighbors(s1.mapped[1:3])

# Interactions entre un ensemble de sommets donnés
stringdb$get_interactions(s1.mapped)


#==========================================#
# Données de coexpression à partir des données complètes
#==========================================#

#------------------------------------------#
# STRINGdb detailed links → neo4j
#------------------------------------------#
# Chargement sous forme de tibble
links.detailed <- read_delim("downloads/511145.protein.links.detailed.v12.0.txt.gz", delim=" ", col_types = "ccnnnnnnnn")

# Génération des fichiers CSV pour l’import
links.detailed %>%
    filter(coexpression>0 & protein1 < protein2) %>%
    dplyr::select(protein1, protein2, coexpression) %>%
    mutate(organism=str_extract(protein1, '^\\d+'), id1=str_extract(protein1, 'b\\d+'), id2=str_extract(protein2, 'b\\d+')) %>%
    dplyr::select(organism:id2,coexpression) %>%
    write_csv("neo4j.import/string.coexpression.csv")

# Import dans neo4j
'
LOAD CSV WITH HEADERS FROM "file:///string.coexpression.csv" AS line
MATCH (g1:Gene),(g2:Gene)
 WHERE g1.id=line.id1 AND g2.id=line.id2
WITH g1,g2, toInteger(line.coexpression) AS value 
MERGE (g1)-[r:STRINGdb]-(g2)
 SET r.coexpression=value
'  %>% cq

# Vérification
'MATCH ()-[r:STRINGdb]-() RETURN count(r)' %>% cq %>% unlist

# Extraction des liens de coexpression d’au moins 0.4
coex <- "MATCH (g1:Gene)-[r:STRINGdb]-(g2:Gene) WHERE r.coexpression>=995 RETURN g1.id, g2.id, r.coexpression" %>% cq 
tibble(id1=coex$g1.id, id2=coex$g2.id, coexpression=coex$r.coexpression)

# Extraction d’un sous-graphe : on utilise ici neo4r et la fonction fournie (call_neo4j) avec le paramètre type="graph" afin de récupérer un sous-graphe.
g.coexpr <- 'MATCH p = ()-[r:STRINGdb]-() WHERE r.coexpression>=995 RETURN p'  %>% cq(result="graph")

# Manipulation pour igraph pour récupérer les propriétés des sommets sous forme de tibble
g.coexpr$nodes <-  g.coexpr$nodes |> 
    map(\(node) c(neo_id=node$id, node$properties)) %>% 
    bind_rows()
dtcoexprnodes <- g.coexpr$nodes %>% dted()
savePath <- paste0("generated.data/",pwd,"/proteins_coexpression_nodes_table.html")
saveWidget(dtcoexprnodes, savePath, selfcontained = TRUE)

# De même pour les arcs/arêtes du graphe à l’aide de la fonction unnest_relationships pour le passage dans igraph : besoin de réordonner les liens
g.coexpr$relationships <- g.coexpr$relationships |> 
    map(\(edge) c(startNode=edge$startNode, endNode=edge$endNode, type=edge$type, edge$properties)) %>% 
    bind_rows
dtcoexprrelation <- g.coexpr$relationships %>% dted()
savePath <- paste0("generated.data/",pwd,"/proteins_coexpression_relations_table.html")
saveWidget(dtcoexprrelation, savePath, selfcontained = TRUE)


#==========================================#
# Graphe
#==========================================#
# Chargement de la librairies et positionnement des valeurs par défaut pour certains paramètres
if (!require("igraph", quietly = TRUE)){install.packages("igraph", repos='https://mirror.ibcp.fr/pub/CRAN/')}
library(igraph)

igraph.options(vertex.color=NA)
igraph.options(vertex.label.cex=.6) # font size
igraph.options(vertex.label.family='sans')
igraph.options(vertex.size=2)
igraph.options(edge.label.cex=.6)
igraph.options(edge.label.family='sans')

# Création du graphe dans igraph à partir des sommets et des relations
# remarque : la première colonne du df passé comme vertices est considérée comme l’identifiant des sommes, donc neo_id ici
g.coexpr = graph_from_data_frame(d=g.coexpr$relationships, directed=FALSE, vertices = g.coexpr$nodes)

# Plot
savePath <- paste0("generated.data/",pwd,"/coexpression.png")
png(savePath, width = 1800, height = 1600, res = 200)
plot(g.coexpr, vertex.label=NA, main='coexpression')
dev.off()


#==========================================#
# Autres scores StringDB à intégrer
#==========================================#

#------------------------------------------#
# experimental
#------------------------------------------#
# Génération des fichiers CSV pour l’import
links.detailed %>%
    filter(experimental > 0 & protein1 < protein2) %>%
    dplyr::select(protein1, protein2, experimental) %>%
    mutate(
        organism = stringr::str_extract(protein1, '^\\d+'),
        id1 = stringr::str_extract(protein1, 'b\\d+'),
        id2 = stringr::str_extract(protein2, 'b\\d+')
    ) %>%
    dplyr::select(organism:id2, experimental) %>%
    write_csv("neo4j.import/string.experimental.csv")

# Import dans neo4j
'
LOAD CSV WITH HEADERS FROM "file:///string.experimental.csv" AS line
MATCH (g1:Gene),(g2:Gene)
    WHERE g1.id = line.id1 AND g2.id = line.id2
WITH g1,g2, toInteger(line.experimental) AS value 
MERGE (g1)-[r:STRINGdb_experimental]-(g2)
    SET r.experimental = value
' %>% cq

# Vérification
'MATCH ()-[r:STRINGdb_experimental]-() RETURN count(r)' %>% cq %>% unlist

#------------------------------------------#
# neighborhood
#------------------------------------------#
# Génération des fichiers CSV pour l’import
links.detailed %>%
    filter(neighborhood > 0 & protein1 < protein2) %>%
    dplyr::select(protein1, protein2, neighborhood) %>%
    mutate(
        organism = stringr::str_extract(protein1, '^\\d+'),
        id1 = stringr::str_extract(protein1, 'b\\d+'),
        id2 = stringr::str_extract(protein2, 'b\\d+')
    ) %>%
    dplyr::select(organism:id2, neighborhood) %>%
    write_csv("neo4j.import/string.neighborhood.csv")

# Import dans neo4j
'
LOAD CSV WITH HEADERS FROM "file:///string.neighborhood.csv" AS line
MATCH (g1:Gene),(g2:Gene)
    WHERE g1.id = line.id1 AND g2.id = line.id2
WITH g1,g2, toInteger(line.neighborhood) AS value 
MERGE (g1)-[r:STRINGdb_neighborhood]-(g2)
    SET r.neighborhood = value
' %>% cq

# Vérification
'MATCH ()-[r:STRINGdb_neighborhood]-() RETURN count(r)' %>% cq %>% unlist

#------------------------------------------#
# textmining
#------------------------------------------#
# Génération des fichiers CSV pour l’import
links.detailed %>%
    filter(textmining > 0 & protein1 < protein2) %>%
    dplyr::select(protein1, protein2, textmining) %>%
    mutate(
        organism = stringr::str_extract(protein1, '^\\d+'),
        id1 = stringr::str_extract(protein1, 'b\\d+'),
        id2 = stringr::str_extract(protein2, 'b\\d+')
    ) %>%
    dplyr::select(organism:id2, textmining) %>%
    write_csv("neo4j.import/string.textmining.csv")

# Import dans neo4j
'
LOAD CSV WITH HEADERS FROM "file:///string.textmining.csv" AS line
MATCH (g1:Gene),(g2:Gene)
    WHERE g1.id = line.id1 AND g2.id = line.id2
WITH g1,g2, toInteger(line.textmining) AS value 
MERGE (g1)-[r:STRINGdb_textmining]-(g2)
    SET r.textmining = value
' %>% cq

# Vérification
'MATCH ()-[r:STRINGdb_textmining]-() RETURN count(r)' %>% cq %>% unlist

#------------------------------------------#
# database
#------------------------------------------#
# Génération des fichiers CSV pour l’import
links.detailed %>%
    filter(database > 0 & protein1 < protein2) %>%
    dplyr::select(protein1, protein2, database) %>%
    mutate(
        organism = stringr::str_extract(protein1, '^\\d+'),
        id1 = stringr::str_extract(protein1, 'b\\d+'),
        id2 = stringr::str_extract(protein2, 'b\\d+')
    ) %>%
    dplyr::select(organism:id2, database) %>%
    write_csv("neo4j.import/string.database.csv")

# Import dans neo4j
'
LOAD CSV WITH HEADERS FROM "file:///string.database.csv" AS line
MATCH (g1:Gene),(g2:Gene)
    WHERE g1.id = line.id1 AND g2.id = line.id2
WITH g1,g2, toInteger(line.database) AS value 
MERGE (g1)-[r:STRINGdb_database]-(g2)
    SET r.database = value
' %>% cq

# Vérification
'MATCH ()-[r:STRINGdb_database]-() RETURN count(r)' %>% cq %>% unlist

#------------------------------------------#
# combined_score
#------------------------------------------#
# Génération des fichiers CSV pour l’import
links.detailed %>%
    filter(combined_score > 0 & protein1 < protein2) %>%
    dplyr::select(protein1, protein2, combined_score) %>%
    mutate(
        organism = stringr::str_extract(protein1, '^\\d+'),
        id1 = stringr::str_extract(protein1, 'b\\d+'),
        id2 = stringr::str_extract(protein2, 'b\\d+')
    ) %>%
    dplyr::select(organism:id2, combined_score) %>%
    write_csv("neo4j.import/string.combined_score.csv")

# Import dans neo4j
'
LOAD CSV WITH HEADERS FROM "file:///string.combined_score.csv" AS line
MATCH (g1:Gene),(g2:Gene)
    WHERE g1.id = line.id1 AND g2.id = line.id2
WITH g1,g2, toInteger(line.combined_score) AS value 
MERGE (g1)-[r:STRINGdb_combined]-(g2)
    SET r.combined_score = value
' %>% cq

# Vérification
'MATCH ()-[r:STRINGdb_combined]-() RETURN count(r)' %>% cq %>% unlist


#==========================================#
# Combinaison de matrices de similarité
#==========================================#
# Récupération du graphe g avec les liens de coexpression, neighborhood et experiment à partir de neo4j
# STRINGdb : r.coexpression
# STRINGdb_experimental : r.experimental (noté ppi)
# STRINGdb_neighborhood : r.neighborhood
# STRINGdb_textmining : r.textmining
# STRINGdb_database : r.database
# STRINGdb_combined : r.combined_score

# On fusionne tout ça par paire de gènes (id1,id2)
edges_raw <- "
MATCH (g1:Gene)-[r:STRINGdb]-(g2:Gene)
OPTIONAL MATCH (g1)-[re:STRINGdb_experimental]-(g2)
OPTIONAL MATCH (g1)-[rn:STRINGdb_neighborhood]-(g2)
OPTIONAL MATCH (g1)-[rt:STRINGdb_textmining]-(g2)
OPTIONAL MATCH (g1)-[rd:STRINGdb_database]-(g2)
OPTIONAL MATCH (g1)-[rc:STRINGdb_combined_score]-(g2)
WHERE g1.id < g2.id   // pour ne garder chaque paire qu'une seule fois
RETURN g1.id AS id1,
       g2.id AS id2,
       r.coexpression AS coexpression,
       re.experimental AS ppi,
       rn.neighborhood AS neighborhood,
       rt.textmining AS textmining,
       rd.database AS database,
       rc.combined_score AS combined_score
" %>% cq
# remplissage des NA par 0
edges_raw[is.na(edges_raw)] <- 0

# Mise en forme en tibble
edges_df <- tibble(
    id1 = edges_raw$id1,
    id2 = edges_raw$id2,
    coexpression = as.numeric(edges_raw$coexpression),
    ppi = as.numeric(edges_raw$ppi),
    neighborhood = as.numeric(edges_raw$neighborhood),
    textmining = as.numeric(edges_raw$textmining),
    database = as.numeric(edges_raw$database),
    combined_score = as.numeric(edges_raw$combined_score)
)

# Création du graphe non orienté dans igraph
# Les colonnes id1 / id2 servent d’extrémités des arêtes,
# et les autres colonnes deviennent des attributs d’arêtes.
g <- graph_from_data_frame(
    d = edges_df,
    directed = FALSE,
    vertices = NULL # les sommets sont déduits de id1 / id2
)

# Calcul du score combiné par STRINGdb (source: https://string-db.org/help/faq/#how-are-the-scores-computed)
prior <- 0.041
no_prior <- function(x, prior = 0.041) (ifelse(is.na(x), 0, x) / 1000 - prior) / (1-prior)
s_coexp_nop <- no_prior(E(g)$coexpression)
s_ppi_nop <- no_prior(E(g)$ppi)
s_neighborhood_nop <- no_prior(E(g)$neighborhood)
s_tot_nop <- 1 - (1 - s_coexp_nop) * (1 - s_ppi_nop) * (1 - s_neighborhood_nop)
E(g)$combined_score <- round(1000 * (s_tot_nop + prior *(1 - s_tot_nop)))


#==========================================#
# Prioritisation de gènes
#==========================================#
# Gènes “training”.
# Choix du cycle des citrates (pathway TCA).
# Récupération des gènes impliqués dans le pathway
training.genes <- "MATCH (Pathway {id: 'TCA'})-[:requires]->(g: Gene {organism: 511145}) return g.bnumber " %>% cq %>% .$g %>% unlist

# Gènes “candidats” (ensemble du génome)
candidates <- "MATCH (g:Gene {organism: 511145}) RETURN g.bnumber" %>% cq %>% .$g 

# Distance d’un gène à l’ensemble de référence
score <- function(gene.ids, ref.genes, datasource) {
    # mapping datasource -> type de relation dans Neo4j
    rel_type <- switch(
        datasource,
        coexpression = "STRINGdb",
        experimental = "STRINGdb_experimental",
        neighborhood = "STRINGdb_neighborhood",
        textmining   = "STRINGdb_textmining",
        database     = "STRINGdb_database",
        stop(paste("Datasource inconnue :", datasource))
    )
    ref_str <- paste("'", ref.genes, "'", sep = "", collapse = ",")
    sapply(gene.ids, function(gene.id) {
        query <- paste0(
        "MATCH (g1:Gene {id: '", gene.id, "', organism: 511145})",
        "-[r:", rel_type, "]-",
        "(g2:Gene {organism: 511145}) ",
        "WHERE g2.id IN [", ref_str, "] ",
        "RETURN r.", datasource, " AS score"
        )
        res <- query %>% cq %>% .$score
        ifelse(is.null(res), 0, mean(replace_na(res, 0)))
    })
}

# test
gene.ids <- 'b0001'
gene.ids <- 'b0116'
ref.genes <- training.genes
datasource <- 'coexpression'
score(gene.ids, training.genes, datasource)

# Application à l’ensemble des gènes
scores <- tibble(candidate=candidates) %>%
  mutate(
    coexpression = score(gene.id=candidate, ref.genes=training.genes, datasource='coexpression'),
    experimental = score(gene.id=candidate, ref.genes=training.genes, datasource='experimental'),
    neighborhood = score(gene.id=candidate, ref.genes=training.genes, datasource='neighborhood'),
    textmining = score(gene.id=candidate, ref.genes=training.genes, datasource='textmining'),
    database = score(gene.id=candidate, ref.genes=training.genes, datasource='database')
    )
scores <- scores %>% replace(is.na(.), 0)

# Plots
savePath <- paste0("generated.data/",pwd,"/priorisation_plot.png")
png(savePath, width = 1800, height = 1600, res = 200)
scores %>%
    ggplot(aes(x=coexpression, y=neighborhood)) +
    geom_point(alpha=0.2, aes(color=candidate %in% training.genes, shape=candidate %in% training.genes))
dev.off()

# ACP pour l’affichage
pca = scores %>%
    dplyr::select(-candidate) %>%
    prcomp
pca %>% summary

# Visualisation
library(factoextra)
savePath <- paste0("generated.data/",pwd,"/priorisation_plot-pca.png")
png(savePath, width = 1800, height = 1600, res = 200)
pca %>% fviz_pca_ind(col.ind = scores$candidate %in% training.genes, addEllipses = T, alpha.ind = .2, geom='point')
dev.off()

# Contribution des sources de données
savePath <- paste0("generated.data/",pwd,"/priorisation_plot-pca-contribution.png")
png(savePath, width = 1800, height = 1600, res = 200)
pca %>% fviz_pca_biplot(col.ind = scores$candidate %in% training.genes)
dev.off()

# Analyse discriminante linéaire
if (!require("MASS", quietly = TRUE)){install.packages("MASS", repos='https://mirror.ibcp.fr/pub/CRAN/')}
library(MASS)
mat <- scores %>% dplyr::select(-candidate) %>% as.matrix
model <- lda(x=mat, grouping = scores$candidate %in% training.genes)

# Distribution des scores sous forme de boîtes à moustache pour les 2 classes
savePath <- paste0("generated.data/",pwd,"/priorisation_plot-LDA-boxplot.png")
png(savePath, width = 1800, height = 1600, res = 200)
mat %>% 
    scale( center=T, scale=F ) %*% model$scaling %>% 
    as_tibble %>%
    mutate(gene.id=scores$candidate) %>%
    ggplot(aes(x=LD1, y=gene.id %in% training.genes, color=gene.id %in% training.genes, shape=gene.id %in% training.genes)) +
    geom_violin() +
    geom_boxplot(varwidth = T) +
    geom_jitter(height = 0.2, alpha=0.1, color='grey') +
    theme_light()
dev.off()

# Density
savePath <- paste0("generated.data/",pwd,"/priorisation_plot-density.png")
png(savePath, width = 1800, height = 1600, res = 200)
mat %>% 
    scale( center=T, scale=F ) %*% model$scaling %>% 
    as_tibble %>%
    mutate(gene.id=scores$candidate) %>%
    ggplot(aes(LD1)) + #, color=gene.id %in% training.genes, shape=gene.id %in% training.genes)) +
    geom_density()
dev.off()

#------------------------------------------#
# Tentative d’optimisation
#------------------------------------------#

# Récupération des “scores” entre les candidats et les gènes de référence en une requête
train.str <- paste("'", training.genes,"'", sep = '', collapse = ',')
datasources <- c('coexpression', 'experimental', 'neighborhood', 'textmining', 'database' )
ds.str <- paste(" r.", datasources, sep='', collapse=',')

tb <- paste0("
MATCH (candidate:Gene {organism: 511145}), (training:Gene {organism: 511145})
WHERE training.id IN [", train.str, "] AND candidate.id <> training.id

OPTIONAL MATCH (candidate)-[rc:STRINGdb]-(training)
OPTIONAL MATCH (candidate)-[re:STRINGdb_experimental]-(training)
OPTIONAL MATCH (candidate)-[rn:STRINGdb_neighborhood]-(training)
OPTIONAL MATCH (candidate)-[rt:STRINGdb_textmining]-(training)
OPTIONAL MATCH (candidate)-[rd:STRINGdb_database]-(training)

RETURN DISTINCT
    candidate.id AS candidate_id,
    training.id  AS training_id,
    rc.coexpression  AS coexpression,
    re.experimental  AS experimental,
    rn.neighborhood  AS neighborhood,
    rt.textmining    AS textmining,
    rd.database      AS database
") %>% cq

# Replace NAs with 0s and compute mean (or other) by candidate
scores <- tb %>%
    dplyr::select(-training_id) %>%
    replace(is.na(.), 0) %>%
    group_by(candidate_id) %>%
    summarise(across(all_of(datasources), mean), .groups = "drop")

# Suppression des colonnes entièrement nulles (aucune info)
scores <- scores %>%
    dplyr::select(candidate_id, where(~ any(. != 0)))

# Analyse discriminante linéaire
mat <- scores %>% dplyr::select(-candidate_id) %>% as.matrix
model <- lda(x=mat, grouping = scores$candidate_id %in% training.genes)

# Distribution des scores sous forme de boîtes à moustache pour les 2 classes
savePath <- paste0("generated.data/",pwd,"/priorisation_plot-LDA-boxplot-optimisation.png")
png(savePath, width = 1800, height = 1600, res = 200)
mat %>% 
    scale( center=T, scale=F ) %*% model$scaling %>% 
    as_tibble %>%
    mutate(gene.id=scores$candidate_id, training=gene.id %in% training.genes) %>%
    ggplot(aes(x=LD1, y=training, color=training, shape=training)) +
    geom_violin() +
    geom_jitter(alpha=.5, height = .2)
dev.off()



############################################
# set.M2.8
############################################

set1 <- read_tsv("tmp/set.M2.8.out", col_names=F)
s1 <- scan('query.sets/set.M2.8.txt', character()) 
pwd <- "set.M2.8"

#------------------------------------------#
# Calcul des similarités entre les GOTerms
#------------------------------------------#
simBP <- calculateSimMatrix(
    set1$X1,
    orgdb  = "org.Hs.eg.db", # proxy humain
    ont    = "BP",
    method = "Resnik"
)

# simBP
str(simBP)

scores <- setNames(-log10(set1$X2), set1$X1)
# scores <- setNames(-log10(set1$pvalue), set1$id)
scores

# Clustering des GOTerms
# reducedTerms <- reduceSimMatrix(simBP, scores, threshold=.7, orgdb="org.EcK12.eg.db")
reducedTerms <- reduceSimMatrix(
    simMatrix = simBP,
    scores    = scores,
    threshold = 0.7,
    orgdb     = "org.Hs.eg.db"
)

# heatmapPlot(simBP, reducedTerms, annotateParent=T, annotationLabel="parentTerm", fontsize=6)
# heatmapPlot(simBP, reducedTerms, annotateParent=T, annotationLabel="parentTerm", fontsize=12)
savePath <- paste0("generated.data/",pwd,"/heatmap_GO.png")
png(savePath, width = 1800, height = 1600, res = 200)
heatmapPlot(simBP, reducedTerms, annotateParent=TRUE, annotationLabel="parentTerm", fontsize=9)
dev.off()

# Visualisation avec un nuage de points dans l’espace sémantique
#savePath <- paste0("generated.data/",pwd,"/scatterPlot_GO.png")
#png(savePath, width = 1800, height = 1600, res = 200)
#scatterPlot(simBP, reducedTerms)
#dev.off()
#Erreur dans cmdscale(as.matrix(as.dist(1 - simMatrix)), eig = TRUE, k = 2) : 
#'k' doit être dans {1, 2, ..  n - 1}

# Visualisation avec une treemap
savePath <- paste0("generated.data/",pwd,"/treemapPlot_GO.png")
png(savePath, width = 1800, height = 1600, res = 200)
treemapPlot(reducedTerms)
dev.off()

# Visualisation avec un wordcloud
savePath <- paste0("generated.data/",pwd,"/wordcloudPlot_GO.png")
png(savePath, width = 1800, height = 1600, res = 200)
wordcloudPlot(reducedTerms, min.freq=1)
dev.off()


#==========================================#
# Approches par graphe ou matrice de distances/dissemblances/similarités/probabilités
#==========================================#

#------------------------------------------#
# Génération de matrices de similarité
#------------------------------------------#
# Ajout d’une colonne avec des 1 pour les association présentes dans les données de départ
keyword.mat <- keywords %>% 
    unique %>% 
    mutate(asso=1)
keyword.mat

# Génération d’une table pivot (les valeurs de la colonne keyword deviennent des colonnes)
keyword.tab <- keyword.mat %>% 
    pivot_wider(names_from = keyword, values_from = asso, values_fill = 0) %>% 
    as.data.frame
keyword.tab %>% head

# Transfert des identifiants de gènes de la colonne bnumber aux noms de lignes de la matrice de distance
rownames(keyword.tab) <- keyword.tab$bnumber
keyword.tab %>% head

# Calcul de la matrice de distance (entre paires de lignes/gènes) avec la méthode binary
keyword.tab <- keyword.tab %>% dplyr::select(-bnumber)
keyword.dist <- keyword.tab %>% dist(method='binary')
keyword.dist[1:10]

# Sous forme de scores à la StringDB
keyword.scores <- round( (1-keyword.dist) * 1000) %>% as.matrix
keyword.scores[1:10,1:10]

#------------------------------------------#
# Client R à STRINGdb
#------------------------------------------#
kabex <- function(o) o %>% kable(format="html", escape=F) %>% kable_styling(bootstrap_options = c("striped"))
library(DT) 
dted <- function(o) o %>% datatable(rownames=F, filter="top", options=list(pageLength = 10), escape=F)

# Connexion au serveur de STRING (NCBI taxon id de Ecoli est 511145)
# https://string-db.org/cgi/input.pl?input_page_active_form=organism
confidence_threshold <- 333
stringdb <- STRINGdb$new(version='12', species=511145, score_threshold=confidence_threshold, input_directory='downloads')

# Téléchargement du génome/protéome
proteins <- stringdb$get_proteins() %>% as_tibble
dt <- proteins %>% dted()
savePath <- paste0("generated.data/",pwd,"/proteins_table.html")
saveWidget(dt, savePath, selfcontained = TRUE)

# Mapping des identifiants vers ceux de StringDB
s1.mapped <- stringdb$mp(s1)

# Plot du graphe correspondant
stringdb$plot_network(s1.mapped)

# Enrichment
enrichment <- stringdb$get_enrichment(s1.mapped)
dtenrichment <- enrichment %>% dted()
savePath <- paste0("generated.data/",pwd,"/proteins_enrichment_table.html")
saveWidget(dtenrichment, savePath, selfcontained = TRUE)

# Sommets voisins d’un ou plusieurs sommets donnés
stringdb$get_neighbors(s1.mapped[1:3])

# Interactions entre un ensemble de sommets donnés
stringdb$get_interactions(s1.mapped)


#==========================================#
# Données de coexpression à partir des données complètes
#==========================================#

#------------------------------------------#
# STRINGdb detailed links → neo4j
#------------------------------------------#
# Chargement sous forme de tibble
links.detailed <- read_delim("downloads/511145.protein.links.detailed.v12.0.txt.gz", delim=" ", col_types = "ccnnnnnnnn")

# Génération des fichiers CSV pour l’import
links.detailed %>%
    filter(coexpression>0 & protein1 < protein2) %>%
    dplyr::select(protein1, protein2, coexpression) %>%
    mutate(organism=str_extract(protein1, '^\\d+'), id1=str_extract(protein1, 'b\\d+'), id2=str_extract(protein2, 'b\\d+')) %>%
    dplyr::select(organism:id2,coexpression) %>%
    write_csv("neo4j.import/string.coexpression.csv")

# Import dans neo4j
'
LOAD CSV WITH HEADERS FROM "file:///string.coexpression.csv" AS line
MATCH (g1:Gene),(g2:Gene)
 WHERE g1.id=line.id1 AND g2.id=line.id2
WITH g1,g2, toInteger(line.coexpression) AS value 
MERGE (g1)-[r:STRINGdb]-(g2)
 SET r.coexpression=value
'  %>% cq

# Vérification
'MATCH ()-[r:STRINGdb]-() RETURN count(r)' %>% cq %>% unlist

# Extraction des liens de coexpression d’au moins 0.4
coex <- "MATCH (g1:Gene)-[r:STRINGdb]-(g2:Gene) WHERE r.coexpression>=995 RETURN g1.id, g2.id, r.coexpression" %>% cq 
tibble(id1=coex$g1.id, id2=coex$g2.id, coexpression=coex$r.coexpression)

# Extraction d’un sous-graphe : on utilise ici neo4r et la fonction fournie (call_neo4j) avec le paramètre type="graph" afin de récupérer un sous-graphe.
g.coexpr <- 'MATCH p = ()-[r:STRINGdb]-() WHERE r.coexpression>=995 RETURN p'  %>% cq(result="graph")

# Manipulation pour igraph pour récupérer les propriétés des sommets sous forme de tibble
g.coexpr$nodes <-  g.coexpr$nodes |> 
    map(\(node) c(neo_id=node$id, node$properties)) %>% 
    bind_rows()
dtcoexprnodes <- g.coexpr$nodes %>% dted()
savePath <- paste0("generated.data/",pwd,"/proteins_coexpression_nodes_table.html")
saveWidget(dtcoexprnodes, savePath, selfcontained = TRUE)

# De même pour les arcs/arêtes du graphe à l’aide de la fonction unnest_relationships pour le passage dans igraph : besoin de réordonner les liens
g.coexpr$relationships <- g.coexpr$relationships |> 
    map(\(edge) c(startNode=edge$startNode, endNode=edge$endNode, type=edge$type, edge$properties)) %>% 
    bind_rows
dtcoexprrelation <- g.coexpr$relationships %>% dted()
savePath <- paste0("generated.data/",pwd,"/proteins_coexpression_relations_table.html")
saveWidget(dtcoexprrelation, savePath, selfcontained = TRUE)


#==========================================#
# Graphe
#==========================================#
# Chargement de la librairies et positionnement des valeurs par défaut pour certains paramètres
if (!require("igraph", quietly = TRUE)){install.packages("igraph", repos='https://mirror.ibcp.fr/pub/CRAN/')}
library(igraph)

igraph.options(vertex.color=NA)
igraph.options(vertex.label.cex=.6) # font size
igraph.options(vertex.label.family='sans')
igraph.options(vertex.size=2)
igraph.options(edge.label.cex=.6)
igraph.options(edge.label.family='sans')

# Création du graphe dans igraph à partir des sommets et des relations
# remarque : la première colonne du df passé comme vertices est considérée comme l’identifiant des sommes, donc neo_id ici
g.coexpr = graph_from_data_frame(d=g.coexpr$relationships, directed=FALSE, vertices = g.coexpr$nodes)

# Plot
savePath <- paste0("generated.data/",pwd,"/coexpression.png")
png(savePath, width = 1800, height = 1600, res = 200)
plot(g.coexpr, vertex.label=NA, main='coexpression')
dev.off()


#==========================================#
# Autres scores StringDB à intégrer
#==========================================#

#------------------------------------------#
# experimental
#------------------------------------------#
# Génération des fichiers CSV pour l’import
links.detailed %>%
    filter(experimental > 0 & protein1 < protein2) %>%
    dplyr::select(protein1, protein2, experimental) %>%
    mutate(
        organism = stringr::str_extract(protein1, '^\\d+'),
        id1 = stringr::str_extract(protein1, 'b\\d+'),
        id2 = stringr::str_extract(protein2, 'b\\d+')
    ) %>%
    dplyr::select(organism:id2, experimental) %>%
    write_csv("neo4j.import/string.experimental.csv")

# Import dans neo4j
'
LOAD CSV WITH HEADERS FROM "file:///string.experimental.csv" AS line
MATCH (g1:Gene),(g2:Gene)
    WHERE g1.id = line.id1 AND g2.id = line.id2
WITH g1,g2, toInteger(line.experimental) AS value 
MERGE (g1)-[r:STRINGdb_experimental]-(g2)
    SET r.experimental = value
' %>% cq

# Vérification
'MATCH ()-[r:STRINGdb_experimental]-() RETURN count(r)' %>% cq %>% unlist

#------------------------------------------#
# neighborhood
#------------------------------------------#
# Génération des fichiers CSV pour l’import
links.detailed %>%
    filter(neighborhood > 0 & protein1 < protein2) %>%
    dplyr::select(protein1, protein2, neighborhood) %>%
    mutate(
        organism = stringr::str_extract(protein1, '^\\d+'),
        id1 = stringr::str_extract(protein1, 'b\\d+'),
        id2 = stringr::str_extract(protein2, 'b\\d+')
    ) %>%
    dplyr::select(organism:id2, neighborhood) %>%
    write_csv("neo4j.import/string.neighborhood.csv")

# Import dans neo4j
'
LOAD CSV WITH HEADERS FROM "file:///string.neighborhood.csv" AS line
MATCH (g1:Gene),(g2:Gene)
    WHERE g1.id = line.id1 AND g2.id = line.id2
WITH g1,g2, toInteger(line.neighborhood) AS value 
MERGE (g1)-[r:STRINGdb_neighborhood]-(g2)
    SET r.neighborhood = value
' %>% cq

# Vérification
'MATCH ()-[r:STRINGdb_neighborhood]-() RETURN count(r)' %>% cq %>% unlist

#------------------------------------------#
# textmining
#------------------------------------------#
# Génération des fichiers CSV pour l’import
links.detailed %>%
    filter(textmining > 0 & protein1 < protein2) %>%
    dplyr::select(protein1, protein2, textmining) %>%
    mutate(
        organism = stringr::str_extract(protein1, '^\\d+'),
        id1 = stringr::str_extract(protein1, 'b\\d+'),
        id2 = stringr::str_extract(protein2, 'b\\d+')
    ) %>%
    dplyr::select(organism:id2, textmining) %>%
    write_csv("neo4j.import/string.textmining.csv")

# Import dans neo4j
'
LOAD CSV WITH HEADERS FROM "file:///string.textmining.csv" AS line
MATCH (g1:Gene),(g2:Gene)
    WHERE g1.id = line.id1 AND g2.id = line.id2
WITH g1,g2, toInteger(line.textmining) AS value 
MERGE (g1)-[r:STRINGdb_textmining]-(g2)
    SET r.textmining = value
' %>% cq

# Vérification
'MATCH ()-[r:STRINGdb_textmining]-() RETURN count(r)' %>% cq %>% unlist

#------------------------------------------#
# database
#------------------------------------------#
# Génération des fichiers CSV pour l’import
links.detailed %>%
    filter(database > 0 & protein1 < protein2) %>%
    dplyr::select(protein1, protein2, database) %>%
    mutate(
        organism = stringr::str_extract(protein1, '^\\d+'),
        id1 = stringr::str_extract(protein1, 'b\\d+'),
        id2 = stringr::str_extract(protein2, 'b\\d+')
    ) %>%
    dplyr::select(organism:id2, database) %>%
    write_csv("neo4j.import/string.database.csv")

# Import dans neo4j
'
LOAD CSV WITH HEADERS FROM "file:///string.database.csv" AS line
MATCH (g1:Gene),(g2:Gene)
    WHERE g1.id = line.id1 AND g2.id = line.id2
WITH g1,g2, toInteger(line.database) AS value 
MERGE (g1)-[r:STRINGdb_database]-(g2)
    SET r.database = value
' %>% cq

# Vérification
'MATCH ()-[r:STRINGdb_database]-() RETURN count(r)' %>% cq %>% unlist

#------------------------------------------#
# combined_score
#------------------------------------------#
# Génération des fichiers CSV pour l’import
links.detailed %>%
    filter(combined_score > 0 & protein1 < protein2) %>%
    dplyr::select(protein1, protein2, combined_score) %>%
    mutate(
        organism = stringr::str_extract(protein1, '^\\d+'),
        id1 = stringr::str_extract(protein1, 'b\\d+'),
        id2 = stringr::str_extract(protein2, 'b\\d+')
    ) %>%
    dplyr::select(organism:id2, combined_score) %>%
    write_csv("neo4j.import/string.combined_score.csv")

# Import dans neo4j
'
LOAD CSV WITH HEADERS FROM "file:///string.combined_score.csv" AS line
MATCH (g1:Gene),(g2:Gene)
    WHERE g1.id = line.id1 AND g2.id = line.id2
WITH g1,g2, toInteger(line.combined_score) AS value 
MERGE (g1)-[r:STRINGdb_combined]-(g2)
    SET r.combined_score = value
' %>% cq

# Vérification
'MATCH ()-[r:STRINGdb_combined]-() RETURN count(r)' %>% cq %>% unlist


#==========================================#
# Combinaison de matrices de similarité
#==========================================#
# Récupération du graphe g avec les liens de coexpression, neighborhood et experiment à partir de neo4j
# STRINGdb : r.coexpression
# STRINGdb_experimental : r.experimental (noté ppi)
# STRINGdb_neighborhood : r.neighborhood
# STRINGdb_textmining : r.textmining
# STRINGdb_database : r.database
# STRINGdb_combined : r.combined_score

# On fusionne tout ça par paire de gènes (id1,id2)
edges_raw <- "
MATCH (g1:Gene)-[r:STRINGdb]-(g2:Gene)
OPTIONAL MATCH (g1)-[re:STRINGdb_experimental]-(g2)
OPTIONAL MATCH (g1)-[rn:STRINGdb_neighborhood]-(g2)
OPTIONAL MATCH (g1)-[rt:STRINGdb_textmining]-(g2)
OPTIONAL MATCH (g1)-[rd:STRINGdb_database]-(g2)
OPTIONAL MATCH (g1)-[rc:STRINGdb_combined_score]-(g2)
WHERE g1.id < g2.id   // pour ne garder chaque paire qu'une seule fois
RETURN g1.id AS id1,
       g2.id AS id2,
       r.coexpression AS coexpression,
       re.experimental AS ppi,
       rn.neighborhood AS neighborhood,
       rt.textmining AS textmining,
       rd.database AS database,
       rc.combined_score AS combined_score
" %>% cq
# remplissage des NA par 0
edges_raw[is.na(edges_raw)] <- 0

# Mise en forme en tibble
edges_df <- tibble(
    id1 = edges_raw$id1,
    id2 = edges_raw$id2,
    coexpression = as.numeric(edges_raw$coexpression),
    ppi = as.numeric(edges_raw$ppi),
    neighborhood = as.numeric(edges_raw$neighborhood),
    textmining = as.numeric(edges_raw$textmining),
    database = as.numeric(edges_raw$database),
    combined_score = as.numeric(edges_raw$combined_score)
)

# Création du graphe non orienté dans igraph
# Les colonnes id1 / id2 servent d’extrémités des arêtes,
# et les autres colonnes deviennent des attributs d’arêtes.
g <- graph_from_data_frame(
    d = edges_df,
    directed = FALSE,
    vertices = NULL # les sommets sont déduits de id1 / id2
)

# Calcul du score combiné par STRINGdb (source: https://string-db.org/help/faq/#how-are-the-scores-computed)
prior <- 0.041
no_prior <- function(x, prior = 0.041) (ifelse(is.na(x), 0, x) / 1000 - prior) / (1-prior)
s_coexp_nop <- no_prior(E(g)$coexpression)
s_ppi_nop <- no_prior(E(g)$ppi)
s_neighborhood_nop <- no_prior(E(g)$neighborhood)
s_tot_nop <- 1 - (1 - s_coexp_nop) * (1 - s_ppi_nop) * (1 - s_neighborhood_nop)
E(g)$combined_score <- round(1000 * (s_tot_nop + prior *(1 - s_tot_nop)))


#==========================================#
# Prioritisation de gènes
#==========================================#
# Gènes “training”.
# Choix du cycle des citrates (pathway TCA).
# Récupération des gènes impliqués dans le pathway
training.genes <- "MATCH (Pathway {id: 'TCA'})-[:requires]->(g: Gene {organism: 511145}) return g.bnumber " %>% cq %>% .$g %>% unlist

# Gènes “candidats” (ensemble du génome)
candidates <- "MATCH (g:Gene {organism: 511145}) RETURN g.bnumber" %>% cq %>% .$g 

# Distance d’un gène à l’ensemble de référence
score <- function(gene.ids, ref.genes, datasource) {
    # mapping datasource -> type de relation dans Neo4j
    rel_type <- switch(
        datasource,
        coexpression = "STRINGdb",
        experimental = "STRINGdb_experimental",
        neighborhood = "STRINGdb_neighborhood",
        textmining   = "STRINGdb_textmining",
        database     = "STRINGdb_database",
        stop(paste("Datasource inconnue :", datasource))
    )
    ref_str <- paste("'", ref.genes, "'", sep = "", collapse = ",")
    sapply(gene.ids, function(gene.id) {
        query <- paste0(
        "MATCH (g1:Gene {id: '", gene.id, "', organism: 511145})",
        "-[r:", rel_type, "]-",
        "(g2:Gene {organism: 511145}) ",
        "WHERE g2.id IN [", ref_str, "] ",
        "RETURN r.", datasource, " AS score"
        )
        res <- query %>% cq %>% .$score
        ifelse(is.null(res), 0, mean(replace_na(res, 0)))
    })
}

# test
gene.ids <- 'b0001'
gene.ids <- 'b0116'
ref.genes <- training.genes
datasource <- 'coexpression'
score(gene.ids, training.genes, datasource)

# Application à l’ensemble des gènes
scores <- tibble(candidate=candidates) %>%
  mutate(
    coexpression = score(gene.id=candidate, ref.genes=training.genes, datasource='coexpression'),
    experimental = score(gene.id=candidate, ref.genes=training.genes, datasource='experimental'),
    neighborhood = score(gene.id=candidate, ref.genes=training.genes, datasource='neighborhood'),
    textmining = score(gene.id=candidate, ref.genes=training.genes, datasource='textmining'),
    database = score(gene.id=candidate, ref.genes=training.genes, datasource='database')
    )
scores <- scores %>% replace(is.na(.), 0)

# Plots
savePath <- paste0("generated.data/",pwd,"/priorisation_plot.png")
png(savePath, width = 1800, height = 1600, res = 200)
scores %>%
    ggplot(aes(x=coexpression, y=neighborhood)) +
    geom_point(alpha=0.2, aes(color=candidate %in% training.genes, shape=candidate %in% training.genes))
dev.off()

# ACP pour l’affichage
pca = scores %>%
    dplyr::select(-candidate) %>%
    prcomp
pca %>% summary

# Visualisation
library(factoextra)
savePath <- paste0("generated.data/",pwd,"/priorisation_plot-pca.png")
png(savePath, width = 1800, height = 1600, res = 200)
pca %>% fviz_pca_ind(col.ind = scores$candidate %in% training.genes, addEllipses = T, alpha.ind = .2, geom='point')
dev.off()

# Contribution des sources de données
savePath <- paste0("generated.data/",pwd,"/priorisation_plot-pca-contribution.png")
png(savePath, width = 1800, height = 1600, res = 200)
pca %>% fviz_pca_biplot(col.ind = scores$candidate %in% training.genes)
dev.off()

# Analyse discriminante linéaire
if (!require("MASS", quietly = TRUE)){install.packages("MASS", repos='https://mirror.ibcp.fr/pub/CRAN/')}
library(MASS)
mat <- scores %>% dplyr::select(-candidate) %>% as.matrix
model <- lda(x=mat, grouping = scores$candidate %in% training.genes)

# Distribution des scores sous forme de boîtes à moustache pour les 2 classes
savePath <- paste0("generated.data/",pwd,"/priorisation_plot-LDA-boxplot.png")
png(savePath, width = 1800, height = 1600, res = 200)
mat %>% 
    scale( center=T, scale=F ) %*% model$scaling %>% 
    as_tibble %>%
    mutate(gene.id=scores$candidate) %>%
    ggplot(aes(x=LD1, y=gene.id %in% training.genes, color=gene.id %in% training.genes, shape=gene.id %in% training.genes)) +
    geom_violin() +
    geom_boxplot(varwidth = T) +
    geom_jitter(height = 0.2, alpha=0.1, color='grey') +
    theme_light()
dev.off()

# Density
savePath <- paste0("generated.data/",pwd,"/priorisation_plot-density.png")
png(savePath, width = 1800, height = 1600, res = 200)
mat %>% 
    scale( center=T, scale=F ) %*% model$scaling %>% 
    as_tibble %>%
    mutate(gene.id=scores$candidate) %>%
    ggplot(aes(LD1)) + #, color=gene.id %in% training.genes, shape=gene.id %in% training.genes)) +
    geom_density()
dev.off()

#------------------------------------------#
# Tentative d’optimisation
#------------------------------------------#

# Récupération des “scores” entre les candidats et les gènes de référence en une requête
train.str <- paste("'", training.genes,"'", sep = '', collapse = ',')
datasources <- c('coexpression', 'experimental', 'neighborhood', 'textmining', 'database' )
ds.str <- paste(" r.", datasources, sep='', collapse=',')

tb <- paste0("
MATCH (candidate:Gene {organism: 511145}), (training:Gene {organism: 511145})
WHERE training.id IN [", train.str, "] AND candidate.id <> training.id

OPTIONAL MATCH (candidate)-[rc:STRINGdb]-(training)
OPTIONAL MATCH (candidate)-[re:STRINGdb_experimental]-(training)
OPTIONAL MATCH (candidate)-[rn:STRINGdb_neighborhood]-(training)
OPTIONAL MATCH (candidate)-[rt:STRINGdb_textmining]-(training)
OPTIONAL MATCH (candidate)-[rd:STRINGdb_database]-(training)

RETURN DISTINCT
    candidate.id AS candidate_id,
    training.id  AS training_id,
    rc.coexpression  AS coexpression,
    re.experimental  AS experimental,
    rn.neighborhood  AS neighborhood,
    rt.textmining    AS textmining,
    rd.database      AS database
") %>% cq

# Replace NAs with 0s and compute mean (or other) by candidate
scores <- tb %>%
    dplyr::select(-training_id) %>%
    replace(is.na(.), 0) %>%
    group_by(candidate_id) %>%
    summarise(across(all_of(datasources), mean), .groups = "drop")

# Suppression des colonnes entièrement nulles (aucune info)
scores <- scores %>%
    dplyr::select(candidate_id, where(~ any(. != 0)))

# Analyse discriminante linéaire
mat <- scores %>% dplyr::select(-candidate_id) %>% as.matrix
model <- lda(x=mat, grouping = scores$candidate_id %in% training.genes)

# Distribution des scores sous forme de boîtes à moustache pour les 2 classes
savePath <- paste0("generated.data/",pwd,"/priorisation_plot-LDA-boxplot-optimisation.png")
png(savePath, width = 1800, height = 1600, res = 200)
mat %>% 
    scale( center=T, scale=F ) %*% model$scaling %>% 
    as_tibble %>%
    mutate(gene.id=scores$candidate_id, training=gene.id %in% training.genes) %>%
    ggplot(aes(x=LD1, y=training, color=training, shape=training)) +
    geom_violin() +
    geom_jitter(alpha=.5, height = .2)
dev.off()


#==========================================#
# END SCRIPT9
#==========================================#
