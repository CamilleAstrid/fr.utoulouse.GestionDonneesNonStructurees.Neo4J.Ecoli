#!/bin/bash

#==========================================#
# MAIN SCRIPT
#==========================================#
# WARNING ! Avant exécution, il faut générer le fichier query.sets/liste.sets.txt dans lequel se trouve la liste des sets de requête à utiliser
# WARNING ! Avant exécution, il faut télécharger les données depuis EcoCyc
# https://www.ecocyc.org/groups/export?id=biocyc13-55140-3842501533&tsv-type=FRAMES
# https://www.ecocyc.org/groups/export?id=biocyc17-55140-3842483872&tsv-type=FRAMES
# https://www.ecocyc.org/groups/export?id=biocyc17-55140-3842483291&tsv-type=FRAMES


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

    printf "[%s%s] %3d%% (%d/%d) %s\n" "$bar_filled" "$bar_empty" "$percent" "$current" "$total" "$extra"
}

# Nombre total de sets
total_tests=7
# Compteur de tests en cours
current_test=0

#------------------------------------------#
# Exécution des différents scripts dans l'ordre
#------------------------------------------#
echo "Lancement du pipeline d'analyse..."

# Étape 1 : Création de l'environnement, de la structure des répertoires et téléchargement des données
current_test=$((current_test+1))
progress_bar "$current_test" "$total_tests" "Étape 1/7 : Création de l'environnement, de la structure des répertoires et téléchargement des données"
./scripts/script1.sh > step1.log 2>&1

# Étape 2 : Préparation des ensembles de référence
current_test=$((current_test+1))
progress_bar "$current_test" "$total_tests" "Étape 2/7 : Préparation des ensembles de référence"
R -q --vanilla < scripts/script2.R > step2.log 2>&1

# Étape 3 : Préparation des ensembles de requête et lancement de l'enrichissement
current_test=$((current_test+1))
progress_bar "$current_test" "$total_tests" "Étape 3/7 : Préparation des ensembles de requête et lancement de l'enrichissement"
./scripts/script4.sh > step3.log 2>&1

# Étape 4 : Lancement de Neo4J
current_test=$((current_test+1))
progress_bar "$current_test" "$total_tests" "Étape 4/7 : Lancement de Neo4J"
./scripts/script5.sh > step4.log 2>&1

# Étape 5 : Importation des données dans Neo4J
current_test=$((current_test+1))
progress_bar "$current_test" "$total_tests" "Étape 5/7 : Importation des données dans Neo4J"
R -q --vanilla < scripts/script6.R > step5.log 2>&1

# Étape 6 : Exécution des tests d'enrichissement avec Neo4J
current_test=$((current_test+1))
progress_bar "$current_test" "$total_tests" "Étape 6/7 : Exécution des tests d'enrichissement avec Neo4J"
./scripts/script8.sh > step6.log 2>&1

# Étape 7 : Calculs de similarité et priorisation des gènes candidats
current_test=$((current_test+1))
progress_bar "$current_test" "$total_tests" "Étape 7/7 : Calculs de similarité et priorisation des gènes candidats"
R -q --vanilla < scripts/script9.R > step7.log 2>&1

#==========================================#
# END MAIN SCRIPT
#==========================================#