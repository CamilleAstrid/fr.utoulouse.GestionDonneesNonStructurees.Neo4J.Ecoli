#!/bin/bash

#==========================================#
# SCRIPT5
#==========================================#
# Script pour lancer Neo4J

#------------------------------------------#
# Prétraitement des fichiers à importer dans Neo4J
#------------------------------------------#
cp generated.data/ncbi.bnumber.location.tsv neo4j.import/


#------------------------------------------#
# Lancement de Neo4J
#------------------------------------------#
podman run --rm -d -it --name neo4j \
   -p 7474:7474 \
   -p 7687:7687 \
   -v ./neo4j.data:/data \
   -v ./neo4j.import:/import \
   -e NEO4J_AUTH=neo4j/omybioinfo \
   docker.io/neo4j:5.26.12-community
# -it : mode interactif a été remplacé par -d pour détacher le conteneur et le laisser tourner en arrière-plan

# récupération des logs
#podman logs -f neo4j

# arrêt du conteneur
#podman stop neo4j

#sudo apt-get install firefox
function openUrl() {
   local URL="$1";
   firefox $URL || xdg-open $URL || sensible-browser $URL || x-www-browser $URL || gnome-open $URL;
}
openUrl http://localhost:7474/ & # authentification : neo4j / omybioinfo

#Ouverture d’un shell cq sur le conteneur en cours (depuis le shell)
#podman exec -it neo4j ./bin/cypher-shell

#==========================================#
# END SCRIPT5
#==========================================#
