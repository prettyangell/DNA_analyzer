# Analyseur ADN

## Description

L'Analyseur ADN est une application Python développée avec la bibliothèque **PyQt5** qui propose diverses fonctionnalités pour l'analyse des séquences ADN. Il permet aux utilisateurs d'effectuer des opérations sur une ou plusieurs séquence d'ADN.

## Fonctionnalités

- Importer des séquences ADN depuis un fichier FASTA ou générer une séquence aléatoire.
- Valider l'intégrité des séquences ADN.
- Calculer les fréquences nucléotidiques dans une séquence ADN.
- Transcrire les séquences ADN en ARN.
- Traduire les séquences ARN en protéines à l'aide d'une table de codons prédéfinie.
- Calculer le complément inverse d'une séquence ADN.
- Calculer le pourcentage de contenu GC d'une séquence ADN.
- Analyser les fréquences de codons dans une séquence ADN.
- Effectuer des mutations aléatoires sur une séquence ADN.
- Rechercher des motifs au sein d'une séquence ADN.
- Générer une séquence d'ADN consensus et une matrice de profil à partir d'un fichier FASTA.

## Les videos de tests

- **test 1:** génération de chaine d'ADN aléatoire et test de toutes les fonctionnalitées (sauf la matrice profil et la chaine consensus)
  [Watch the Video](https://drive.google.com/file/d/10ppMwEa_9qwS9sSiwgMr4rCdRT3cZOHg/view?usp=drive_link)
- **test 2:** lecture de plusieurs chaines d'ADN d'un fichier fasta et test de toutes les fonctionnalitées.
  [Watch the Video](https://drive.google.com/file/d/18Sk7ESid24_LbkWYTmAUlRcetAZLPfXv/view?usp=drive_link)
- **test 3:** tests des cas particuliers:[Watch the Video](https://drive.google.com/file/d/1ypELhE7Hktf3W9guPoUob6_Q7vZcs-z3/view?usp=drive_link)
  - le fichier fasta contien une seule chaine
  - le motif et en minuscule
  - le motif n'est pas valid (trop long ou contient d'autres caractéres)
- **test 4:**d'autres tests particuliers:[Watch the Video](https://drive.google.com/file/d/1XA-IXjM6GX2mOf4RJyt7TPRmIpmbZ50t/view?usp=drive_link)
  - le fichier n'est pas de type FASTA
  - le fichier fasta contien une chaine fausse(invalide)
  - le fichier fasta contien des chianes de differentes tailles(on peut pas avoir la matrice profil et la chaine consensus)

## Contributeurs

- Bouhraoua Yousra Hind

## Utilisation

Exécutez l'application en lançant le fichier DNA_analyzer.py :

```bash
python DNA_analyzer.py
```

1- Choisissez entre importer une séquence ADN depuis un **fichier FASTA**(vous pouvez utilisez le fichier **"input.txt"**) ou générer une séquence aléatoire.
2- Effectuez diverses opérations depuis le menu pour analyser la séquence ADN.
3- Visualisez les résultats et enregistrez-les dans un fichier texte si nécessaire(vous pouvez utilisez le fichier **"Results.txt"**).
