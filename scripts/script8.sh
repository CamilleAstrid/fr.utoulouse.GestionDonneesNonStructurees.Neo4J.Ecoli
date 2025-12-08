#!/bin/bash

#==========================================#
# SCRIPT8
#==========================================#
# Script pour exécuter les tests d'enrichissement

fichier="query.sets/liste.sets.txt"

#------------------------------------------#
#  Fonction barre de progression
#------------------------------------------#
progress_bar() {
    local current=$1
    local total=$2
    local extra="$3" # nom du set
    local width=40 # largeur de la barre

    # pourcentage
    local percent=$(( current * 100 / total ))
    # nombre de blocs "plein"
    local filled=$(( current * width / total ))
    # nombre de blocs "vides"
    local empty=$(( width - filled ))

    # construire la barre
    local bar_filled
    local bar_empty
    bar_filled=$(printf "%${filled}s" | tr ' ' '#')
    bar_empty=$(printf "%${empty}s" | tr ' ' '-')

    # \r pour réécrire sur la même ligne, pas de \n à la fin
    printf "\r[%s%s] %3d%% (%d/%d) %s" "$bar_filled" "$bar_empty" "$percent" "$current" "$total" "$extra"
}

# Nombre total de sets
total_sets=$(wc -l < "$fichier")
# Nombre de tests d'enrichissement à exécuter
total_tests=$((total_sets * 3))
# Compteur de tests en cours
current_test=0

#------------------------------------------#
# Exécution des tests d'enrichissement
#------------------------------------------#
while read set;
do \
    # mettre à jour la barre de progression
    current_set=$((current_set+1))
    progress_bar "$current_set" "$total_tests" "$set"
    python3 scripts/blastsets.neo4j.py -q query.sets/$set.txt -t Pathway -c | column -ts $'\t' > tests/$set.neo4j.enrichments.tsv
    
    current_set=$((current_set+1))
    progress_bar "$current_set" "$total_tests" "$set"
    python3 scripts/blastsets.neo4j.py -q query.sets/$set.txt -t GOTerm --species 511145 -c | column -ts $'\t' > tests/$set.neo4j.GOTerm.enrichments.tsv
    
    current_set=$((current_set+1))
    progress_bar "$current_set" "$total_tests" "$set"
    python3 scripts/blastsets.neo4j.py -q query.sets/$set.txt -t GOTerm  -c > tmp/$set.out ;
done < $fichier

#==========================================#
# END SCRIPT8
#==========================================#