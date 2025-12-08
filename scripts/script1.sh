#!/bin/bash

#==========================================#
# SCRIPT1
#==========================================#

#------------------------------------------#
# Création ou chargement de l'environnement
#------------------------------------------#
ENV_NAME="gdnsap"

eval "$(conda shell.bash hook)"

if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    echo "L'environnement '$ENV_NAME' existe déjà, activation..."
    conda activate "$ENV_NAME"
else
    echo "Création de l'environnement '$ENV_NAME'..."
    mamba create -n "$ENV_NAME" r-base \
                                r-tidyverse \
                                r-reticulate \
                                r-codetools \
                                r-dt \
                                r-kableextra \
                                bioconductor-stringdb \
                                r-rmdformats \
                                r-factoextra \
                                scipy \
                                igraph \
                                py2neo \
                                monotonic \
                                packaging \
                                numpy \
                                pandas \
                                scikit-learn \
                                biopython \
                                neo4j-python-driver \
                                r-markdown -y
    echo "Activation de l'environnement '$ENV_NAME'..."
    conda activate "$ENV_NAME"
fi

pip -m pip install neo4j

#------------------------------------------#
# Création de la structure des répertoires
#------------------------------------------#
mkdir -p downloads                         # données d'origine téléchargées
mkdir -p reference.sets                    # ensembles de référence générés pour la recherche d'enrichissement
mkdir -p query.sets                        # ensembles requête
mkdir -p neo4j.import                      # répertoire pour les fichiers de données à importer dans neo4j
mkdir -p neo4j.data                        # répertoire où neo4j stocke ses données
mkdir -p scripts                           # scripts de recherche d'enrichissement et de prioritization
mkdir -p tests                             # résultats des tests
mkdir -p generated.data                    # données générées par les scripts
mkdir -p generated.data/set.01             # données générées pour le set.01
mkdir -p generated.data/set.02             # données générées pour le set.02
mkdir -p generated.data/set.03             # données générées pour le set.03
mkdir -p generated.data/set.M2.8           # données générées pour le set.M2.8
mkdir -p tmp                               # fichiers temporaires

#------------------------------------------#
# Téléchargements des données du projet
#------------------------------------------#

# Téléchargement du protéome d'Escherichia coli K12 MG1655 depuis UniProt
curl -o downloads/uniprotkb_proteome_UP000000625.tsv.gz "https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Cid%2Cgene_oln%2Cgo_id%2Cxref_interpro%2Ckeyword&format=tsv&query=%28proteome%3AUP000000625%29"

# Visualisation
zcat downloads/uniprotkb_proteome_UP000000625.tsv.gz | head | column -ts $'\t'

# Téléchargement de la liste InterPro au format court
curl -o downloads/interpro_entry_list_short.tsv "https://ftp.ebi.ac.uk/pub/databases/interpro/releases/latest/entry.list"

# Téléchargement des données depuis EcoCyc
#curl -o downloads/All-instances-of-Genes-in-Escherichia-coli-K-12-substr.-MG1655.tsv "https://www.ecocyc.org/groups/export?id=biocyc13-55140-3842501533&tsv-type=FRAMES"
#curl -o downloads/All-instances-of-Pathways-in-Escherichia-coli-K-12-substr.-MG1655.tsv "https://www.ecocyc.org/groups/export?id=biocyc17-55140-3842483872&tsv-type=FRAMES"
#curl -o downloads/all-transcription-units.tsv "https://www.ecocyc.org/groups/export?id=biocyc17-55140-3842483291&tsv-type=FRAMES"
# Ne s'automatise pas facilement à cause de la gestion des sessions sur EcoCyc
#curl: (56) Recv failure: Connection reset by peer
# Alternative : téléchargement manuel depuis le depot Gitlab
curl -o downloads/All-instances-of-Genes-in-Escherichia-coli-K-12-substr.-MG1655.tsv "https://gitlab.com/rbarriot/data/-/raw/main/M2.GDNS-AP/All-instances-of-Genes-in-Escherichia-coli-K-12-substr.-MG1655.txt"
curl -o downloads/All-instances-of-Pathways-in-Escherichia-coli-K-12-substr.-MG1655.tsv "https://gitlab.com/rbarriot/data/-/raw/main/M2.GDNS-AP/All-instances-of-Pathways-in-Escherichia-coli-K-12-substr.-MG1655.txt"
curl -o downloads/all-transcription-units.tsv "https://gitlab.com/rbarriot/data/-/raw/main/M2.GDNS-AP/all-transcription-units.txt"


# Téléchargement des données du ficher d’annotation du génome depuis NCBI
curl -o downloads/GCF_000005845.2.zip "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_000005845.2/download?include_annotation_type=GENOME_GBFF&hydrated=FULLY_HYDRATED"
unzip -p downloads/GCF_000005845.2.zip ncbi_dataset/data/GCF_000005845.2/genomic.gbff | gzip > downloads/GCF_000005845.2.gbff.gz

# Téléchargement des données de neo4j
curl https://gitlab.com/rbarriot/data/-/raw/main/M2.GDNS-AP/pmid.selected.tsv > neo4j.import/pmid.selected.tsv

# Récupération de la dernière version au format OBO de la Gene Ontology
wget http://current.geneontology.org/ontology/go-basic.obo
mv go-basic.obo downloads/go-basic.obo

# Téléchargement du fichier pour STRINGdb detailed links → neo4j (Ecoli K12 MG1655 taxon id : 511145)
curl -o downloads/511145.protein.links.detailed.v12.0.txt.gz "https://stringdb-downloads.org/download/protein.links.detailed.v12.0/511145.protein.links.detailed.v12.0.txt.gz"

#------------------------------------------#
# Prétraitement
#------------------------------------------#
cp scripts/script1-bis.py scripts/gbff.to.bnumber.GeneID.py
cp scripts/script1-ter.py scripts/gbff.bnumber.location.py
cp scripts/script3.py scripts/blastsets.json.py
cp scripts/script7.py scripts/blastsets.neo4j.py


python3 scripts/gbff.to.bnumber.GeneID.py downloads/GCF_000005845.2.gbff.gz > generated.data/mapping.bnumber_ncbi.tsv
python3 scripts/gbff.bnumber.location.py downloads/GCF_000005845.2.gbff.gz > generated.data/ncbi.bnumber.location.tsv

# La commande suivante remplace le script optionel script-optional.py
curl -s https://gitlab.com/rbarriot/data/-/raw/main/M2.GDNS-AP/bnumber.PMID.tsv > generated.data/bnumber.PMID.tsv

# Réutilisation du code GeneOntology.py pour générer les fichiers correspondants aux sommets et arcs à charger dans neo4j
curl -s https://gitlab.com/rbarriot/data/-/raw/main/M2.GDNS-AP/GeneOntology.py > scripts/GeneOntology.py
python3 scripts/GeneOntology.py downloads/go-basic.obo nodes > neo4j.import/go.nodes.tsv
python3 scripts/GeneOntology.py downloads/go-basic.obo edges is_a > neo4j.import/go.is_a.edges.tsv
python3 scripts/GeneOntology.py downloads/go-basic.obo edges part_of > neo4j.import/go.part_of.edges.tsv

# Récupération des fichiers de sets de requête
fichier="query.sets/liste.sets.txt"
while read set;
do \
    curl -s https://gitlab.com/rbarriot/data/-/raw/main/M2.GDNS-AP/$set.txt > query.sets/$set.txt ;
done < $fichier
# curl -s https://gitlab.com/rbarriot/data/-/raw/main/M2.GDNS-AP/set.01.txt > query.sets/set.01.txt
# curl -s https://gitlab.com/rbarriot/data/-/raw/main/M2.GDNS-AP/set.02.txt > query.sets/set.02.txt
# curl -s https://gitlab.com/rbarriot/data/-/raw/main/M2.GDNS-AP/set.03.txt > query.sets/set.03.txt
# curl -s https://gitlab.com/rbarriot/data/-/raw/main/M2.GDNS-AP/set.M2.8.txt > query.sets/set.M2.8.txt

#==========================================#
# END SCRIPT1
#==========================================#
