# fr.utoulouse.GestionDonneesNonStructurees.Neo4J.Ecoli

## Usage
```bash
git clone
cd

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

# Étape 6 : Exécution des tests d'enrichissement avec Neo4J
./scripts/script8.sh

# Étape 7 : Calculs de similarité et priorisation des gènes candidats
./scripts/script10.sh

```
