#==========================================#
# SCRIPT2
#==========================================#
# Script pour générer les fichiers de référence JSON


setwd("~/GitHub/fr.utoulouse.GestionDonneesNonStructurees.Neo4J.Ecoli")

#------------------------------------------#
# Chargement des librairies
#------------------------------------------#
if (!require("neo2R", quietly = TRUE)){install.packages("neo2R", repos='https://mirror.ibcp.fr/pub/CRAN/')}
library(neo2R)
library(tidyverse)
library(reticulate)
if (!require("jsonlite", quietly = TRUE)){install.packages("jsonlite", repos='https://mirror.ibcp.fr/pub/CRAN/')}
library(jsonlite)

#-----------------------------------------#
# Chargement du fichier et son contenu
#------------------------------------------#
uniprot <- read_tsv("downloads/uniprotkb_proteome_UP000000625.tsv.gz", show_col_types = F)


#------------------------------------------#
# Extraction des bnumbers
#------------------------------------------#
mapping <- uniprot %>% 
	select(Entry, names = `Gene Names (ordered locus)`) %>%
	mutate(bnumber=str_extract(names, 'b\\d+')) %>%
	select(bnumber, uniprotID=Entry) %>%
	filter(!is.na(bnumber)) %>% # 2023 → P0DQD7 and P0A6D5 are lost (no bnumber)
	arrange(bnumber)
bnumber.uniprot <- mapping

#------------------------------------------#
# Chargement de la liste InterPro (format court)
#------------------------------------------#
interpro_list <- read_tsv("downloads/interpro_entry_list_short.tsv", show_col_types = FALSE)
# on ne garde que l'identifiant et la description courte
interpro_desc <- interpro_list %>% transmute(id = ENTRY_AC, desc = ENTRY_NAME)

#------------------------------------------#
# Chargement des fichiers PMID
#------------------------------------------#
# table bnumber <-> PMID
bnum_pmid <- read_tsv("generated.data/bnumber.PMID.tsv", show_col_types = FALSE) %>%
    mutate(
        bnumber = as.character(bnumber),
        PMID = as.character(PMID)
    )

# table de mapping bnumbers valides provenant de NCBI
mapping_ncbi <- read_tsv("generated.data/mapping.bnumber_ncbi.tsv", show_col_types = FALSE) %>%
    mutate(bnumber = as.character(bnumber)) %>%
    select(bnumber) %>%
    distinct()

# on ne garde que les bnumbers présents dans le mapping NCBI
bnum_pmid <- bnum_pmid %>%
    inner_join(mapping_ncbi, by = "bnumber")

#==========================================#
# REFORMATAGE
#==========================================#

#------------------------------------------#
# Reformatage des données : KEYWORDS
#------------------------------------------#
keywords <- uniprot %>% 
	select(uniprotID=Entry, keyword=Keywords) %>%
	right_join(mapping, by="uniprotID") %>% # right join to remove those without bnumber
	separate_rows(keyword, sep=';') %>%
	select(bnumber, keyword) %>%
	arrange(bnumber)

# Sauvegarde des fichiers CSV pour import dans Neo4j
keywords %>% 
    select(keyword) %>%
    unique %>%
    write_csv("neo4j.import/uniprot.keywords.csv")
# Liens Keyword → Gene
keywords %>% write_csv("neo4j.import/uniprot.keywords.genes.csv")

# Regroupement des bnumbers par keywords
ref_sets <- keywords %>% 
	group_by(keyword) %>%
	summarise(count=n(), elements = list(bnumber)) %>%
	ungroup %>%
	filter(count>1)  %>%
	select(id=keyword, desc=count, elements)

# Sauvegarde au format JSON
ref_sets %>% 
	toJSON(pretty = TRUE) %>% 
	write("reference.sets/uniprot.keywords.sets.json")


#------------------------------------------#
# Reformatage des données : INTERPRO
#------------------------------------------#
interpro_annot <- uniprot %>%
	select(uniprotID = Entry,
        interpro = InterPro) %>%
  # on ne garde que les protéines ayant un bnumber
	right_join(mapping, by = "uniprotID") %>%
  # séparer les domaines multiples (séparés par ;)
	separate_rows(interpro, sep = ";") %>%
	mutate(interpro = str_trim(interpro)) %>%
	filter(!is.na(interpro), interpro != "") %>%
	select(bnumber, interpro)

# Regroupement des bnumbers par domaine InterPro
interpro_sets <- interpro_annot %>%
    group_by(interpro) %>%
    summarise(
        elements = list(sort(unique(bnumber))),
        count    = length(unlist(elements)),
        .groups  = "drop"
    ) %>%
    # on ne garde que les domaines présents sur au moins 2 protéines
    filter(count > 1)

# Ajout des descriptions InterPro
ref_sets <- interpro_sets %>%
    left_join(interpro_desc, by = c("interpro" = "id")) %>%
    # si la description est manquante, on met l'ID comme desc par défaut
    mutate(desc = if_else(is.na(desc), interpro, desc)) %>%
    select(id = interpro, desc, elements)

# Sauvegarde au format JSON
ref_sets %>%
    toJSON(pretty = TRUE) %>%
    write("reference.sets/uniprot.InterPro.sets.json")


#------------------------------------------#
# Reformatage des données : GOTerms
#------------------------------------------#
go_annot <- uniprot %>%
    select(uniprotID = Entry,
            go_ids = `Gene Ontology IDs`) %>%
    # ne garder que les protéines ayant un bnumber
    right_join(mapping, by = "uniprotID") %>%
    # séparer les GO multiples (séparés par ;)
    separate_rows(go_ids, sep = ";") %>%
    mutate(go_ids = str_trim(go_ids)) %>%
    filter(!is.na(go_ids), go_ids != "") %>%
    select(bnumber, go_id = go_ids)

# Regroupement des bnumbers par GO term
go_sets <- go_annot %>%
    group_by(go_id) %>%
    summarise(
        elements = list(sort(unique(bnumber))),
        count    = length(unlist(elements)),
        .groups  = "drop"
    ) %>%
    # on ne garde que les GO présents sur au moins 2 protéines (comme pour keywords)
    filter(count > 1)

# Ajout des descriptions en générant les liens AMIGO
descriptGO <- "https://amigo.geneontology.org/amigo/term/"
go_sets <- go_sets %>%
    mutate(desc = paste0(descriptGO, go_id))

# Construction de la table finale au format attendu
ref_sets <- go_sets %>%
    select(id = go_id, desc, elements)

# Sauvegarde au format JSON
ref_sets %>%
	toJSON(pretty = TRUE) %>%
	write("reference.sets/uniprot.GOTerms.sets.json")

GOTerms <- uniprot %>% 
  select(uniprotID=Entry, GOTerm=`Gene Ontology IDs`) %>%
  right_join(mapping, by = join_by(uniprotID)) %>% # right join to remove those without bnumber
  separate_rows(GOTerm, sep='; ') %>%
  select(bnumber, GOTerm) %>%
  arrange(bnumber)
# Save CSV file for Neo4J
GOTerms %>% write_csv("neo4j.import/uniprot.GOTerm.bnumber.csv")

#------------------------------------------#
# Chargement des fichiers EcoCyc
#------------------------------------------#
ecocyc_genes <- read_tsv("downloads/All-instances-of-Genes-in-Escherichia-coli-K-12-substr.-MG1655.tsv", show_col_types = FALSE)
ecocyc_pwys <- read_tsv("downloads/All-instances-of-Pathways-in-Escherichia-coli-K-12-substr.-MG1655.tsv", show_col_types = FALSE)
ecocyc_tus <- read_tsv("downloads/all-transcription-units.tsv", show_col_types = FALSE)

#------------------------------------------#
# Reformatage des données : PATHWAYS
#------------------------------------------#
# mapping EcoCyc gene_id -> bnumber via la table "Genes"
ecocyc_mapping <- ecocyc_genes %>%
    mutate(bnumber = str_extract(Names, "b\\d+")) %>%
    select(gene_id = Genes, bnumber) %>%
    filter(!is.na(bnumber))

pathways_annot <- ecocyc_pwys %>%
    select(pathway_id = `Object ID`,
            pathway_name = `Common-Name`,
            genes = `Genes of pathway`) %>%
    # certains pathways ont plusieurs gènes séparés par " // "
    separate_rows(genes, sep = " // ") %>%
    mutate(genes = str_trim(genes)) %>%
    filter(!is.na(genes), genes != "") %>%
    # rattacher les bnumbers
    left_join(ecocyc_mapping, by = c("genes" = "gene_id")) %>%
    filter(!is.na(bnumber)) %>%
    select(bnumber, pathway_id, pathway_name)

# Regroupement des bnumbers par Pathway
pathway_sets <- pathways_annot %>%
    group_by(pathway_id, pathway_name) %>%
    summarise(
        elements = list(sort(unique(bnumber))),
        count = length(unlist(elements)),
        .groups = "drop"
    ) %>%
    # on ne garde que les pathways avec au moins 2 gènes
    filter(count > 1)

# Construction de la table finale au format attendu
ref_sets <- pathway_sets %>%
    select(id = pathway_id, desc = pathway_name, elements)

# Sauvegarde au format JSON
ref_sets %>%
    toJSON(pretty = TRUE) %>%
    write("reference.sets/ecocyc.Pathways.sets.json")

#------------------------------------------#
# Reformatage des données : TU (transcription units)
#------------------------------------------#
# on réutilise ecocyc_mapping défini plus haut (gene_id -> bnumber)
tu_annot <- ecocyc_tus %>%
    select(TU_id = `Object ID`,
            strand = `Strand`,
            left = `Left`,
            right = `Right`,
            genes = `Genes of transcription unit`) %>%
    separate_rows(genes, sep = " // ") %>%
    mutate(genes = str_trim(genes)) %>%
    filter(!is.na(genes), genes != "") %>%
    left_join(ecocyc_mapping, by = c("genes" = "gene_id")) %>%
    filter(!is.na(bnumber)) %>%
    select(bnumber, TU_id, strand, left, right)

# Regroupement des bnumbers par TU
tu_sets <- tu_annot %>%
    group_by(TU_id, strand, left, right) %>%
    summarise(
        elements = list(sort(unique(bnumber))),
        count    = length(unlist(elements)),
        .groups  = "drop"
    ) %>%
    # on ne garde que les TUs avec au moins 2 gènes
    filter(count > 1)

# Ajout de la description


# Construction de la table finale au format attendu
ref_sets <- tu_sets %>%
    mutate(desc = paste0(strand, ":", left, ":", right)) %>%  # ex: "-:1825688:1825625"
    select(id = TU_id, desc, elements)

# Sauvegarde au format JSON
ref_sets %>%
    toJSON(pretty = TRUE) %>%
    write("reference.sets/ecocyc.TUs.sets.json")

#------------------------------------------#
# Reformatage des données : PMID
#------------------------------------------#
pmid_sets <- bnum_pmid %>%
    group_by(PMID) %>%
    summarise(
        elements = list(sort(unique(bnumber))),
        count    = length(unlist(elements)),
        .groups  = "drop"
    ) %>%
    # on ne garde que les PMIDs associés à au moins 2 gènes (comme pour keywords)
    filter(count > 1)

# Ajout des descriptions en générant les liens AMIGO
descriptPMID <- "https://pubmed.ncbi.nlm.nih.gov/"
pmid_sets <- pmid_sets %>%
    mutate(desc = paste0(descriptPMID, PMID))

# Construction de la table finale au format attendu
ref_sets <- pmid_sets %>%
    select(id = PMID, desc, elements)

# Sauvegarde au format JSON
ref_sets %>%
    toJSON(pretty = TRUE) %>%
    write("reference.sets/ncbi.PMIDs.sets.json")


#==========================================#
# END SCRIPT2
#==========================================#
