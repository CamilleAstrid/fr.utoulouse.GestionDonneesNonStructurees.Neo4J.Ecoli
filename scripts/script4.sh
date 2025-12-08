#!/bin/bash

#==========================================#
# SCRIPT4
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
total_tests=$((total_sets * 7))
# Compteur de tests en cours
current_test=0

#------------------------------------------#
# Exécution des tests d'enrichissement
#------------------------------------------#
while read set;
do \
    # mettre à jour la barre de progression
    current_test=$((current_test+1))
    progress_bar "$current_test" "$total_tests" "$set"

    # résultats attendus pour le nombre de résultats significatifs obtenus pour le set01
    # sans ajustement FDR
    echo '# Enrichissements significatifs pour le' $set > tests/$set.uniprot.keywords.enrichments.count.txt
    python3 scripts/blastsets.json.py -q query.sets/$set.txt -t reference.sets/uniprot.keywords.sets.json | wc -l >> tests/$set.uniprot.keywords.enrichments.count.txt
    # avec ajustement FDR
    echo '# Enrichissements significatifs pour le' $set 'avec ajustement FDR' >> tests/$set.uniprot.keywords.enrichments.count.txt
    python3 scripts/blastsets.json.py --adjust -q query.sets/$set.txt -t reference.sets/uniprot.keywords.sets.json | wc -l >> tests/$set.uniprot.keywords.enrichments.count.txt

    
    # résultats attendus pour l'enrichissement du set01

    current_test=$((current_test+1))
    progress_bar "$current_test" "$total_tests" "$set"
    python3 scripts/blastsets.json.py -q query.sets/$set.txt -t reference.sets/uniprot.keywords.sets.json | column -ts $'\t' > tests/$set.uniprot.keywords.enrichments.tsv

    current_test=$((current_test+1))
    progress_bar "$current_test" "$total_tests" "$set"
    python3 scripts/blastsets.json.py -q query.sets/$set.txt -t reference.sets/uniprot.InterPro.sets.json | column -ts $'\t' > tests/$set.uniprot.InterPro.enrichments.tsv

    current_test=$((current_test+1))
    progress_bar "$current_test" "$total_tests" "$set"
    python3 scripts/blastsets.json.py -q query.sets/$set.txt -t reference.sets/uniprot.GOTerms.sets.json | column -ts $'\t' > tests/$set.uniprot.GOTerms.enrichments.tsv

    current_test=$((current_test+1))
    progress_bar "$current_test" "$total_tests" "$set"
    python3 scripts/blastsets.json.py -q query.sets/$set.txt -t reference.sets/ecocyc.pathways.sets.json | column -ts $'\t' > tests/$set.ecocyc.pathways.enrichments.tsv

    current_test=$((current_test+1))
    progress_bar "$current_test" "$total_tests" "$set"
    python3 scripts/blastsets.json.py -q query.sets/$set.txt -t reference.sets/ecocyc.TUs.sets.json | column -ts $'\t' > tests/$set.ecocyc.TUs.enrichments.tsv

    current_test=$((current_test+1))
    progress_bar "$current_test" "$total_tests" "$set"
    python3 scripts/blastsets.json.py -q query.sets/$set.txt -t reference.sets/ncbi.PMIDs.sets.json -c | column -ts $'\t' > tests/$set.ncbi.PMIDs.enrichments.tsv ;
done < $fichier


#==========================================#
# END SCRIPT4
#==========================================#