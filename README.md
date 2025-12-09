# fr.utoulouse.GestionDonneesNonStructurees.Neo4J.Ecoli

## Description

Ce projet présente une analyse des données de *Escherichia coli K-12 MG1655* avec l'aide de base de données graphes.
Le travail s’inscrit dans le cadre du **Master Bioinformatique et Biologie des Systèmes** de l’Université de Toulouse, unité d'enseignement **Gestion des Données Non Structurées**.

## Structure

### Dossiers et fichiers
* `README.md` : Fichier de présentation du projet (vous y êtes !).
* `LICENSE` : Licence d’utilisation.
* `.gitignore` : Liste des fichiers et/ou dossiers à ignorer pour le git.
* `Rapport.md` : Consignes du rapport de projet.
* `Rapport_RODRIGUES.md` : Version finale du rapport rédigé pour l’UE en markdown.
* `Rapport_RODRIGUES.pdf` : Version finale du rapport rédigé pour l’UE au format pdf.
* `M2BBS_GDNS-APG_schema_base_neo4j.png` : Schéma de la base Neo4j pour le projet.
* **generated.data/** : dossier regroupant les fichiers de données dérivées (générées par les scripts à partir des sources brutes).
* **neo4j.data/** : dossier contenant les fichiers internes de Neo4j (base, transactions, configuration).
* **neo4j.import/** : ensemble des fichiers prêts à être importés dans Neo4j.
* **query.sets/** : définition des ensembles de gènes à analyser.
  * `liste.sets.txt` : liste des sets à traiter.
* **reference.sets/** : ensembles de référence utilisés pour les tests d’enrichissement.
* **scripts/** : dossier regroupant l’ensemble des scripts utilisés pour télécharger, transformer, intégrer et analyser les données.
* **tests/** : résultats textuels des tests d’enrichissement et de cohérence.
* **tmp/** : fichiers temporaires générés lors des analyses.
* **downloads/** : données brutes téléchargées depuis les bases externes.
* **analyses/** : figures finales utilisées pour le rapport et l’interprétation des résultats.

Cette organisation vise à séparer clairement **les données brutes**, **les données dérivées**, **les scripts** et **les résultats**, afin de faciliter la compréhension du pipeline et la reproductibilité des analyses.

## Outils utilisés
* podman v5.5.2
* Neo4j via une image : neo4j:5.26.12-community

## Prérequis
![mamba](https://img.shields.io/badge/mamba-2.1.1-green)
![Python](https://img.shields.io/badge/Python-3.12.11-blue)
![R](https://img.shields.io/badge/R-4.4.3-darkred)

### Langages
* bash
* cypher
* python
* R

### Packages

![Python-package](https://img.shields.io/badge/Python-sys-lightblue)
![Python-package](https://img.shields.io/badge/Python-re-lightblue)
![Python-package](https://img.shields.io/badge/Python-Bio_SeqIO-lightblue)
![Python-package](https://img.shields.io/badge/Python-Bio_Entrez-lightblue)
![Python-package](https://img.shields.io/badge/Python-gzip-lightblue)
![Python-package](https://img.shields.io/badge/Python-argparse-lightblue)
![Python-package](https://img.shields.io/badge/Python-os.path_isfile-lightblue)
![Python-package](https://img.shields.io/badge/Python-json-lightblue)
![Python-package](https://img.shields.io/badge/Python-scipy.stats_binom-lightblue)
![Python-package](https://img.shields.io/badge/Python-scipy.stats_hypergeom-lightblue)
![Python-package](https://img.shields.io/badge/Python-neo4j_GraphDatabase-lightblue)

![R-package](https://img.shields.io/badge/R-neo2R-red)
![R-package](https://img.shields.io/badge/R-tidyverse-red)
![R-package](https://img.shields.io/badge/R-reticulate-red)
![R-package](https://img.shields.io/badge/R-jsonlite-red)
![R-package](https://img.shields.io/badge/R-BiocManager-red)
![R-package](https://img.shields.io/badge/R-rrvgo-red)
![R-package](https://img.shields.io/badge/R-ontologyIndex-red)
![R-package](https://img.shields.io/badge/R-GOSemSim-red)
![R-package](https://img.shields.io/badge/R-org.Hs.eg.db-red)
![R-package](https://img.shields.io/badge/R-STRINGdb-red)
![R-package](https://img.shields.io/badge/R-kableExtra-red)
![R-package](https://img.shields.io/badge/R-htmlwidgets-red)
![R-package](https://img.shields.io/badge/R-DT-red)
![R-package](https://img.shields.io/badge/R-igraph-red)
![R-package](https://img.shields.io/badge/R-factoextra-red)
![R-package](https://img.shields.io/badge/R-MASS-red)

## Installation

```bash
git clone https://https://github.com/CamilleAstrid/fr.utoulouse.GestionDonneesNonStructurees.Neo4J.Ecoli.git
cd fr.utoulouse.GestionDonneesNonStructurees.Neo4J.Ecoli
```

## Usage

>[!WARNING]
>Avant utilisation des scripts, générer un fichier `liste.sets.txt` dans un dossier `query.sets/`.
>L'ensemble des téléchargements et des analyses vont reposer sur ce fichier. Il en existe un exemple dans ce dépôt mais il peut être modifié pour d'autres analyses.

```bash
# Étape 1 : Création de l'environnement, de la structure des répertoires et téléchargement des données
./scripts/script1.sh
mamba activate gdnsap

# Étape 2 : Préparation des ensembles de référence
R -q --vanilla < scripts/script2.R

# Étape 3 : Préparation des ensembles de requête et lancement de l'enrichissement
./scripts/script4.sh

# Étape 4 : Lancement de Neo4J
./scripts/script5.sh

# Étape 5 : Importation des données dans Neo4J
#R -q --vanilla < scripts/script6.R
# A exécuter depuis Rstudio pas à pas (délai de traitement de Neo4J à optimiser)
```

>[!WARNING]
>L'étape 5 est a effectué depuis RStudio (ou autre IDE) avant que la version finale de l'automatisation des scripts ne soit aboutie. Il faut lancer les sections pas à pas à cause du délai de traitement de Neo4J via R via le conteneur (optimisation ultérieure).

```bash
# Étape 6 : Exécution des tests d'enrichissement avec Neo4J
./scripts/script8.sh

# Étape 7 : Calculs de similarité et priorisation des gènes candidats
./scripts/script10.sh
```

>[!NOTE]
>A terme, le script SCRIPT.sh permettra d'exécuter l'ensemble des analyses automatiquement avec redirection des logs.

## Licence
Ce projet et donc l'ensemble des éléments de ce répertoire est sous [licence GPL-v3](https://github.com/CamilleAstrid/fr.utoulouse.GestionDonneesNonStructurees.Neo4J.Ecoli/blob/main/LICENSE) (sauf cas précisé).


## Auteurs

Les données et les codes initiaux sont issus des enseignements de Roland Barriot et issus du dépôt GitLab : https://gitlab.com/rbarriot/data/-/tree/main/M2.GDNS-AP?ref_type=heads.

Les adaptations des scripts, leur automatisation et le rapport sont la propriété intellectuelle de Camille-Astrid Rodrigues.

>[!NOTE]
>Pour toute question, veuillez me contacter par mail : [Camille-Astrid Rodrigues](mailto:camilleastrid.cr@gmail.com)   
>Si des ajustements ou des ajouts sont nécessaires, n'hésitez pas à me le signaler !
