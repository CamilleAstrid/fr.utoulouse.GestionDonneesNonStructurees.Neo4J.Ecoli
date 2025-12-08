#!/bin/bash

#==========================================#
# SCRIPT10
#==========================================#
# Script pour lancer le script9.R dans une boucle

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
total_tests=$(wc -l < "$fichier")
# Compteur de tests en cours
current_test=0

#------------------------------------------#
# Calculs de similarité et priorisation des gènes candidats
#------------------------------------------#
while read set;
do \
	current_test=$((current_test+1))
	progress_bar "$current_test" "$total_tests" "$set"
	R -q --vanilla < scripts/script9-auto.R $set > step7.log 2>&1 ;
done < $fichier

#==========================================#
# END SCRIPT10
#==========================================#
