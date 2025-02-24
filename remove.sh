#!/bin/bash

# Vérifie si le nombre d'arguments est correct
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <A> <B>"
    exit 1
fi

# Récupère les arguments
A=$1
B=$2

# Boucle pour supprimer les dossiers de resultatsA à resultatsB
for ((i=A; i<=B; i++)); do
    folder="result${i}"
    if [ -d "$folder" ]; then
        echo "Suppression du dossier $folder"
        rm -r "$folder"
    else
        echo "Le dossier $folder n'existe pas."
    fi
done

echo "Suppression terminée."

