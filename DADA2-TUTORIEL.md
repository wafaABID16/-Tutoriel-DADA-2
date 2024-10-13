Tutoriel DADA2
================

``` r
# Charger le package dada2 pour l'analyse des données de séquençage

library(dada2); packageVersion("dada2")
```

    ## Loading required package: Rcpp

    ## [1] '1.28.0'

``` r
# Définir le chemin vers le répertoire contenant les fichiers fastq

path <- "~/MiSeq_SOP" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)
```

    ##  [1] "F3D0_S188_L001_R1_001.fastq"   "F3D0_S188_L001_R2_001.fastq"  
    ##  [3] "F3D1_S189_L001_R1_001.fastq"   "F3D1_S189_L001_R2_001.fastq"  
    ##  [5] "F3D141_S207_L001_R1_001.fastq" "F3D141_S207_L001_R2_001.fastq"
    ##  [7] "F3D142_S208_L001_R1_001.fastq" "F3D142_S208_L001_R2_001.fastq"
    ##  [9] "F3D143_S209_L001_R1_001.fastq" "F3D143_S209_L001_R2_001.fastq"
    ## [11] "F3D144_S210_L001_R1_001.fastq" "F3D144_S210_L001_R2_001.fastq"
    ## [13] "F3D145_S211_L001_R1_001.fastq" "F3D145_S211_L001_R2_001.fastq"
    ## [15] "F3D146_S212_L001_R1_001.fastq" "F3D146_S212_L001_R2_001.fastq"
    ## [17] "F3D147_S213_L001_R1_001.fastq" "F3D147_S213_L001_R2_001.fastq"
    ## [19] "F3D148_S214_L001_R1_001.fastq" "F3D148_S214_L001_R2_001.fastq"
    ## [21] "F3D149_S215_L001_R1_001.fastq" "F3D149_S215_L001_R2_001.fastq"
    ## [23] "F3D150_S216_L001_R1_001.fastq" "F3D150_S216_L001_R2_001.fastq"
    ## [25] "F3D2_S190_L001_R1_001.fastq"   "F3D2_S190_L001_R2_001.fastq"  
    ## [27] "F3D3_S191_L001_R1_001.fastq"   "F3D3_S191_L001_R2_001.fastq"  
    ## [29] "F3D5_S193_L001_R1_001.fastq"   "F3D5_S193_L001_R2_001.fastq"  
    ## [31] "F3D6_S194_L001_R1_001.fastq"   "F3D6_S194_L001_R2_001.fastq"  
    ## [33] "F3D7_S195_L001_R1_001.fastq"   "F3D7_S195_L001_R2_001.fastq"  
    ## [35] "F3D8_S196_L001_R1_001.fastq"   "F3D8_S196_L001_R2_001.fastq"  
    ## [37] "F3D9_S197_L001_R1_001.fastq"   "F3D9_S197_L001_R2_001.fastq"  
    ## [39] "filtered"                      "HMP_MOCK.v35.fasta"           
    ## [41] "Mock_S280_L001_R1_001.fastq"   "Mock_S280_L001_R2_001.fastq"  
    ## [43] "mouse.dpw.metadata"            "mouse.time.design"            
    ## [45] "stability.batch"               "stability.files"

``` r
# Les noms des fichiers fastq pour les lectures avant et arrière suivent le format : SAMPLENAME_R1_001.fastq et SAMPLENAME_R2_001.fastq

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extraire les noms d'échantillons en supposant que les noms de fichiers ont le format : SAMPLENAME_XXX.fastq

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
```

``` r
# Visualiser le profil de qualité des deux premiers fichiers R1

plotQualityProfile(fnFs[1:2])
```

![](DADA2-TUTORIEL_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
# Visualiser le profil de qualité des deux premiers fichiers R2

plotQualityProfile(fnRs[1:2])
```

![](DADA2-TUTORIEL_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
# Définir les chemins pour les fichiers filtrés

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
```

``` r
# Assigner les noms d'échantillons aux fichiers filtrés pour un meilleur suivi

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
```

    ##                               reads.in reads.out
    ## F3D0_S188_L001_R1_001.fastq       7793      7113
    ## F3D1_S189_L001_R1_001.fastq       5869      5299
    ## F3D141_S207_L001_R1_001.fastq     5958      5463
    ## F3D142_S208_L001_R1_001.fastq     3183      2914
    ## F3D143_S209_L001_R1_001.fastq     3178      2941
    ## F3D144_S210_L001_R1_001.fastq     4827      4312

``` r
# Estimer les erreurs dans les lectures avant à partir des fichiers filtrés

errF <- learnErrors(filtFs, multithread=TRUE)
```

    ## 33514080 total bases in 139642 reads from 20 samples will be used for learning the error rates.

``` r
# Estimer les erreurs dans les lectures arrière à partir des fichiers filtrés

errR <- learnErrors(filtRs, multithread=TRUE)
```

    ## 22342720 total bases in 139642 reads from 20 samples will be used for learning the error rates.

``` r
# Visualiser les taux d'erreur estimés pour les lectures avant

plotErrors(errF, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis
    ## Transformation introduced infinite values in continuous y-axis

![](DADA2-TUTORIEL_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
# Appliquer le modèle d'erreur estimé aux lectures avant pour générer des    séquences amplicon

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
```

    ## Sample 1 - 7113 reads in 1979 unique sequences.
    ## Sample 2 - 5299 reads in 1639 unique sequences.
    ## Sample 3 - 5463 reads in 1477 unique sequences.
    ## Sample 4 - 2914 reads in 904 unique sequences.
    ## Sample 5 - 2941 reads in 939 unique sequences.
    ## Sample 6 - 4312 reads in 1267 unique sequences.
    ## Sample 7 - 6741 reads in 1756 unique sequences.
    ## Sample 8 - 4560 reads in 1438 unique sequences.
    ## Sample 9 - 15637 reads in 3590 unique sequences.
    ## Sample 10 - 11413 reads in 2762 unique sequences.
    ## Sample 11 - 12017 reads in 3021 unique sequences.
    ## Sample 12 - 5032 reads in 1566 unique sequences.
    ## Sample 13 - 18075 reads in 3707 unique sequences.
    ## Sample 14 - 6250 reads in 1479 unique sequences.
    ## Sample 15 - 4052 reads in 1195 unique sequences.
    ## Sample 16 - 7369 reads in 1832 unique sequences.
    ## Sample 17 - 4765 reads in 1183 unique sequences.
    ## Sample 18 - 4871 reads in 1382 unique sequences.
    ## Sample 19 - 6504 reads in 1709 unique sequences.
    ## Sample 20 - 4314 reads in 897 unique sequences.

``` r
# Appliquer le modèle d'erreur estimé aux lectures arrière pour générer des séquences amplicon

dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
```

    ## Sample 1 - 7113 reads in 1660 unique sequences.
    ## Sample 2 - 5299 reads in 1349 unique sequences.
    ## Sample 3 - 5463 reads in 1335 unique sequences.
    ## Sample 4 - 2914 reads in 853 unique sequences.
    ## Sample 5 - 2941 reads in 880 unique sequences.
    ## Sample 6 - 4312 reads in 1286 unique sequences.
    ## Sample 7 - 6741 reads in 1803 unique sequences.
    ## Sample 8 - 4560 reads in 1265 unique sequences.
    ## Sample 9 - 15637 reads in 3414 unique sequences.
    ## Sample 10 - 11413 reads in 2522 unique sequences.
    ## Sample 11 - 12017 reads in 2771 unique sequences.
    ## Sample 12 - 5032 reads in 1415 unique sequences.
    ## Sample 13 - 18075 reads in 3290 unique sequences.
    ## Sample 14 - 6250 reads in 1390 unique sequences.
    ## Sample 15 - 4052 reads in 1134 unique sequences.
    ## Sample 16 - 7369 reads in 1635 unique sequences.
    ## Sample 17 - 4765 reads in 1084 unique sequences.
    ## Sample 18 - 4871 reads in 1161 unique sequences.
    ## Sample 19 - 6504 reads in 1502 unique sequences.
    ## Sample 20 - 4314 reads in 732 unique sequences.

``` r
# Affiche le résultat du processus DADA appliqué aux lectures avant pour le premier échantillon

dadaFs[[1]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 128 sequence variants were inferred from 1979 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

``` r
# Fusionner les résultats des lectures avant et arrière pour créer des séquences d'amplicon assemblées 

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
```

    ## 6540 paired-reads (in 107 unique pairings) successfully merged out of 6891 (in 197 pairings) input.

    ## 5028 paired-reads (in 101 unique pairings) successfully merged out of 5190 (in 157 pairings) input.

    ## 4986 paired-reads (in 81 unique pairings) successfully merged out of 5267 (in 166 pairings) input.

    ## 2595 paired-reads (in 52 unique pairings) successfully merged out of 2754 (in 108 pairings) input.

    ## 2553 paired-reads (in 60 unique pairings) successfully merged out of 2785 (in 119 pairings) input.

    ## 3646 paired-reads (in 55 unique pairings) successfully merged out of 4109 (in 157 pairings) input.

    ## 6079 paired-reads (in 81 unique pairings) successfully merged out of 6514 (in 198 pairings) input.

    ## 3968 paired-reads (in 91 unique pairings) successfully merged out of 4388 (in 187 pairings) input.

    ## 14233 paired-reads (in 143 unique pairings) successfully merged out of 15355 (in 352 pairings) input.

    ## 10528 paired-reads (in 120 unique pairings) successfully merged out of 11165 (in 278 pairings) input.

    ## 11154 paired-reads (in 137 unique pairings) successfully merged out of 11797 (in 298 pairings) input.

    ## 4349 paired-reads (in 85 unique pairings) successfully merged out of 4802 (in 179 pairings) input.

    ## 17431 paired-reads (in 153 unique pairings) successfully merged out of 17812 (in 272 pairings) input.

    ## 5850 paired-reads (in 81 unique pairings) successfully merged out of 6095 (in 159 pairings) input.

    ## 3716 paired-reads (in 86 unique pairings) successfully merged out of 3894 (in 147 pairings) input.

    ## 6865 paired-reads (in 99 unique pairings) successfully merged out of 7191 (in 187 pairings) input.

    ## 4426 paired-reads (in 67 unique pairings) successfully merged out of 4603 (in 127 pairings) input.

    ## 4576 paired-reads (in 101 unique pairings) successfully merged out of 4739 (in 174 pairings) input.

    ## 6092 paired-reads (in 109 unique pairings) successfully merged out of 6315 (in 173 pairings) input.

    ## 4269 paired-reads (in 20 unique pairings) successfully merged out of 4281 (in 28 pairings) input.

``` r
# Affiche les premières lignes du dataframe de fusion pour le premier échantillon

head(mergers[[1]])
```

    ##                                                                                                                                                                                                                                                       sequence
    ## 1 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCAAGTCAGCGGTAAAATTGAGAGGCTCAACCTCTTCGAGCCGTTGAAACTGGTTTTCTTGAGTGAGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCTCAACTGACGCTCATGCACGAAAGTGTGGGTATCGAACAGG
    ## 2 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGCCAAGTCAGCGGTAAAATTGCGGGGCTCAACCCCGTACAGCCGTTGAAACTGCCGGGCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTGAGGCACGAAAGTGCGGGGATCAAACAGG
    ## 3 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGCTGTTAAGTCAGCGGTCAAATGTCGGGGCTCAACCCCGGCCTGCCGTTGAAACTGGCGGCCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCCCGACTGACGCTGAGGCACGAAAGCGTGGGTATCGAACAGG
    ## 4 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGCTTTTAAGTCAGCGGTAAAAATTCGGGGCTCAACCCCGTCCGGCCGTTGAAACTGGGGGCCTTGAGTGGGCGAGAAGAAGGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCCTTCCGGCGCCCTACTGACGCTGAGGCACGAAAGTGCGGGGATCGAACAGG
    ## 5 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGACTCTCAAGTCAGCGGTCAAATCGCGGGGCTCAACCCCGTTCCGCCGTTGAAACTGGGAGCCTTGAGTGCGCGAGAAGTAGGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCCTACCGGCGCGCAACTGACGCTCATGCACGAAAGCGTGGGTATCGAACAGG
    ## 6 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGATGCCAAGTCAGCGGTAAAAAAGCGGTGCTCAACGCCGTCGAGCCGTTGAAACTGGCGTTCTTGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTGAGGCACGAAAGCGTGGGTATCGAACAGG
    ##   abundance forward reverse nmatch nmismatch nindel prefer accept
    ## 1       579       1       1    148         0      0      1   TRUE
    ## 2       470       2       2    148         0      0      2   TRUE
    ## 3       449       3       4    148         0      0      1   TRUE
    ## 4       430       4       3    148         0      0      2   TRUE
    ## 5       345       5       6    148         0      0      1   TRUE
    ## 6       282       6       5    148         0      0      2   TRUE

``` r
#  créer une table de séquence (ASV) à partir des séquences fusionnées

seqtab <- makeSequenceTable(mergers)

# Affiche les dimensions de la table de séquences

dim(seqtab)
```

    ## [1]  20 293

``` r
#  Compte le nombre de séquences pour chaque longueur et affiche la distribution

table(nchar(getSequences(seqtab)))
```

    ## 
    ## 251 252 253 254 255 
    ##   1  88 196   6   2

``` r
# Identifier et retirer les séquences chimériques (artificielles) basées sur une méthode de consensus

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
```

    ## Identified 61 bimeras out of 293 input sequences.

``` r
# # Afficher les dimensions de la table de séquences après suppression des chimériques

dim(seqtab.nochim)
```

    ## [1]  20 232

``` r
# Calculer la proportion de séquences non chimériques par rapport à la table de séquences originale

sum(seqtab.nochim)/sum(seqtab)
```

    ## [1] 0.9640374

``` r
# Extraire les séquences uniques à l'aide de la fonction getUniques() et les compter avec la fonction sum()

getN <- function(x) sum(getUniques(x))

## Créer une matrice de suivi avec les statistiques de chaque étape du traitement

track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

# Renommer les colonnes de la matrice 'track' avec des descriptions appropriées

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")

# Assigner les noms d'échantillons comme noms de lignes pour la matrice 'track'

rownames(track) <- sample.names

# Affiche les premières lignes de la matrice 'track' pour inspection 

head(track)
```

    ##        input filtered denoisedF denoisedR merged nonchim
    ## F3D0    7793     7113      6976      6979   6540    6528
    ## F3D1    5869     5299      5227      5239   5028    5017
    ## F3D141  5958     5463      5331      5357   4986    4863
    ## F3D142  3183     2914      2799      2830   2595    2521
    ## F3D143  3178     2941      2822      2868   2553    2519
    ## F3D144  4827     4312      4151      4228   3646    3507

``` r
# Créer le répertoire "tax" dans ton répertoire personnel (home) s'il n'existe pas

if (!dir.exists("~/tax")) dir.create("~/tax")

# Télécharger le fichier Silva v132 dans le répertoire "tax"

download.file("https://zenodo.org/record/1172783/files/silva_nr_v132_train_set.fa.gz",
              destfile = "~/tax/silva_nr_v132_train_set.fa.gz", method = "auto")

# Vérifier que le fichier a été correctement téléchargé

file.exists("~/tax/silva_nr_v132_train_set.fa.gz")
```

    ## [1] TRUE

``` r
# Assigner la taxonomie aux séquences non chimériques à l'aide d'une base de données de formation SILVA

taxa <- assignTaxonomy(seqtab.nochim, "~/tax/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
```

``` r
# Copier les résultats de taxonomie dans une nouvelle variable 

taxa.print <- taxa 

#Retirer les noms de lignes pour l'affichage uniquement

rownames(taxa.print) <- NULL

# Afficher les premières lignes des résultats de taxonomie

head(taxa.print)
```

    ##      Kingdom    Phylum          Class         Order           Family          
    ## [1,] "Bacteria" "Bacteroidetes" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [2,] "Bacteria" "Bacteroidetes" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [3,] "Bacteria" "Bacteroidetes" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [4,] "Bacteria" "Bacteroidetes" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [5,] "Bacteria" "Bacteroidetes" "Bacteroidia" "Bacteroidales" "Bacteroidaceae"
    ## [6,] "Bacteria" "Bacteroidetes" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ##      Genus        
    ## [1,] NA           
    ## [2,] NA           
    ## [3,] NA           
    ## [4,] NA           
    ## [5,] "Bacteroides"
    ## [6,] NA

``` r
# Télecharger le package DECIPHER pour le traitement des séquences

library(DECIPHER); packageVersion("DECIPHER")
```

    ## Loading required package: Biostrings

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
    ##     table, tapply, union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## Loading required package: XVector

    ## Loading required package: GenomeInfoDb

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

    ## Loading required package: RSQLite

    ## Loading required package: parallel

    ## [1] '2.28.0'

``` r
# Extraire les séquences uniques de l'échantillon "Mock"

unqs.mock <- seqtab.nochim["Mock",]

# Conserver uniquement les ASVs présentes, en les triant par ordre décroissant 
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) 

# Afficher le nombre d'ASVs présentes dans la communauté "Mock"

cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
```

    ## DADA2 inferred 20 sample sequences present in the Mock community.

``` r
# Télecharger les séquences de référence à partir d'un fichier FASTA contenant des séquences de la communauté "Mock"

mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))

#Compter le nombre d'ASVs détectées qui correspondent aux séquences de référence

match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))

# Afficher le nombre d'ASVs qui correspondent aux séquences de référence

cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")
```

    ## Of those, 20 were exact matches to the expected reference sequences.

``` r
# Télecharger le package phyloseq pour l'analyse de la diversité microbienne

library(phyloseq); packageVersion("phyloseq")
```

    ## 
    ## Attaching package: 'phyloseq'

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     distance

    ## [1] '1.44.0'

``` r
# Télécharger le package Biostrings pour la manipulation des séquences biologiques

library(Biostrings); packageVersion("Biostrings")
```

    ## [1] '2.68.1'

``` r
# Télecharger le package ggplot2 pour la visualisation des données

library(ggplot2); packageVersion("ggplot2")
```

    ## [1] '3.4.3'

``` r
# Définir le thème par défaut pour les graphiques ggplot2 à un thème en noir et blanc

theme_set(theme_bw())
```

``` r
#Extraire les noms des échantillons à partir de la table des séquences non chimériques

samples.out <- rownames(seqtab.nochim)

#Extraire le sujet à partir du nom de l'échantillon, en supposant que le format soit "SOMETHING_D1", "SOMETHING_D2", etc

subject <- sapply(strsplit(samples.out, "D"), `[`, 1)

# Extraire le sexe à partir du nom du sujet

gender <- substr(subject,1,1)

# Extraire l'identifiant du sujet en supprimant la première lettre

subject <- substr(subject,2,999)

#  Extraire le jour à partir du nom de l'échantillon

day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))

# Créer un DataFrame contenant les informations sur les sujets, le sexe et le jour

samdf <- data.frame(Subject=subject, Gender=gender, Day=day)

# Ajouter une colonne "When" pour indiquer si l'échantillon a été collecté tôt ou tard

samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"

# Définir les noms de lignes du DataFrame pour correspondre aux noms des échantillons

rownames(samdf) <- samples.out
```

``` r
## Créer un objet phyloseq à partir de la table d'abondance OTU, des données d'échantillon et des données de taxonomie

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

#Supprimer l'échantillon "Mock" de l'objet phyloseq

ps <- prune_samples(sample_names(ps) != "Mock", ps) 
```

``` r
# Créer un objet DNAStringSet à partir des noms des taxons dans l'objet phyloseq

dna <- Biostrings::DNAStringSet(taxa_names(ps))

#Nommer les séquences dans l'objet DNAStringSet avec les noms des taxons

names(dna) <- taxa_names(ps)

# Fusionner l'objet DNAStringSet avec l'objet phyloseq

ps <- merge_phyloseq(ps, dna)

# Renommer les ASVs dans l'objet phyloseq

taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 232 taxa and 19 samples ]
    ## sample_data() Sample Data:       [ 19 samples by 4 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 232 taxa by 6 taxonomic ranks ]
    ## refseq()      DNAStringSet:      [ 232 reference sequences ]

``` r
# Créer un graphique de richesse avec les indices de diversité Shannon et Simpson, en coloriant par le moment de collecte (Early/Late)

plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When")
```

    ## Warning in estimate_richness(physeq, split = TRUE, measures = measures): The data you have provided does not have
    ## any singletons. This is highly suspicious. Results of richness
    ## estimates (for example) are probably unreliable, or wrong, if you have already
    ## trimmed low-abundance taxa from the data.
    ## 
    ## We recommended that you find the un-trimmed data and retry.

![](DADA2-TUTORIEL_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

``` r
# # Transformer les données en proportions pour le calcul des distances de Bray-Curtis

ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))

# Calculer une ordination NMDS basée sur les distances de Bray-Curtis

ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
```

    ## Run 0 stress 0.08043117 
    ## Run 1 stress 0.08043116 
    ## ... New best solution
    ## ... Procrustes: rmse 1.081946e-06  max resid 2.420815e-06 
    ## ... Similar to previous best
    ## Run 2 stress 0.08076337 
    ## ... Procrustes: rmse 0.01049328  max resid 0.03229206 
    ## Run 3 stress 0.09477198 
    ## Run 4 stress 0.1432873 
    ## Run 5 stress 0.08076336 
    ## ... Procrustes: rmse 0.01048315  max resid 0.03225926 
    ## Run 6 stress 0.1010631 
    ## Run 7 stress 0.08076337 
    ## ... Procrustes: rmse 0.01050296  max resid 0.0323237 
    ## Run 8 stress 0.1358136 
    ## Run 9 stress 0.08616061 
    ## Run 10 stress 0.1010629 
    ## Run 11 stress 0.08616061 
    ## Run 12 stress 0.08616061 
    ## Run 13 stress 0.08616061 
    ## Run 14 stress 0.08043117 
    ## ... Procrustes: rmse 6.06345e-06  max resid 1.568702e-05 
    ## ... Similar to previous best
    ## Run 15 stress 0.08076341 
    ## ... Procrustes: rmse 0.01058296  max resid 0.03258477 
    ## Run 16 stress 0.08076337 
    ## ... Procrustes: rmse 0.01050059  max resid 0.03231604 
    ## Run 17 stress 0.1212044 
    ## Run 18 stress 0.08076338 
    ## ... Procrustes: rmse 0.01052621  max resid 0.0323993 
    ## Run 19 stress 0.08076337 
    ## ... Procrustes: rmse 0.01051393  max resid 0.03235946 
    ## Run 20 stress 0.09477199 
    ## *** Best solution repeated 2 times

``` r
# Visualiser l'ordination NMDS avec les distances de Bray-Curtis, colorée par "When"

plot_ordination(ps.prop, ord.nmds.bray, color="When", title="Bray NMDS")
```

![](DADA2-TUTORIEL_files/figure-gfm/unnamed-chunk-35-1.png)<!-- -->

``` r
# Extraire les noms des 20 OTUs les plus abondants dans l'objet phyloseq

top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]

#  Transformer les données d'abondance pour les échantillons en proportions

ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))

# Prendre uniquement les taxons qui sont dans le top 20

ps.top20 <- prune_taxa(top20, ps.top20)

# Créer un graphique en barres représentant les familles des 20 OTUs les plus abondants, en fonction du jour et du moment de collecte

plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~When, scales="free_x")
```

![](DADA2-TUTORIEL_files/figure-gfm/unnamed-chunk-36-1.png)<!-- -->
