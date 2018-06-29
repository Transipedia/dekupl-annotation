# TODO

## Annotation

[x] Ajouter le fichier BED d'output des alignements
[x] Ajouter le fichier Loci d'output de l'annotation
[ ] Ajouter un graphique sur l'architecture de DEkupl-annot
[ ] Génerer le tableau des ontologies automatiquement dans le README
[ ] Ajouter les dernieres colonnes manquantes (voir fichier Drive)
[ ] Afficher les headers dans le fichier output (comme un VCF)
[ ] Faire des tests pour GSNAP et SWitche afin de vérifier que les dépendances sont bien installées
    [ ] GSNAP
    [ ] Samtools
    [ ] R & DESeq2
[ ] Logs de GSNAP dans le dossier de résultats et pas un fichier temporaire
[ ] Choix de l'annotations si plusieurs chevauchantes :
    1. exon over intron
    2. Taille du chevauchement
    3. Longeur du gène
[ ] Annotation des jonctions connus (dans le GFF)
[ ] Ontologie SNV: on conserve la règle pour rétro-compatibilité
[x] Ajout d'une colonne avec le CIGAR
[x] On renomme 5p et 4p en downstream/upstream
[ ] Ajout d'une colonne avec la liste des features GFF qui overlappent le contig
[ ] Ontologie UTR  : pas de modifications
[ ] Benchmark des procédure pour la recherche de fusion : GNSAP 2-step, STAR, Blast
[ ] Ajout du mapping BLAST pour les primers
[ ] Ajout du mapping BLAST pour les contigs non mappés par GSNAP
[ ] Ajouter des time sur les logs pour faire des stats

## Indexing

[ ] Un seul index pour dekupl-annotation et dekupl-run
[x] On vérifie à la constriction de l'index que les fichier sont compatibles : 
    1. Mêmes références dans le GFF et le fasta
    2. Présence de features "gene"
    3. Vérification que les ID sont uniques
    4. Warning sur les données manquantes utilisés par l'annotation (ex: Name dans le GFF)