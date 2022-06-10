# Structure du dossier :

# path 
#  | Directement dans path : les fichiers CCS fwd et rev mélangés 
#  |------- /forward_filtered : dossier où l'on va stocker les fwd "filtrés"
#  |------- /reverse_filtered : dossier où l'on va stocker les rev "filtrés"
#  |------- /output : dossier où stocker les fichiers de sorties
#              |---------> seqtab.rds (sequence table créée en cours de script)
#              |---------> taxonomy.rds (table avec la taxonomie)
#  |--------/tax : dossier avec la table de concordance pour l'assignement de la taxonomie
###############|---------> silva_nr99_v138.1_train_set.fa.gz (table pour la concordance avec la taxonomie)

# Normalement "forward_filtered" et "reverse_filtered" sont créés automatiquement
# par la fonction filterAndTrim(), d'après les chemins specifiés dans le script
# mais il m'est arrivé de devoir les créer manuellement pour que ça fonctionne.
# Ces dossiers vont se remplir automatiquement suite à l'application de filterAndTrim()

# Le dossier "output" est à créer "manuellement" avant de lancer le script.

# Enfin il faudrait créer le dossier "tax" manuellement avant de lancer le script
# et mettre à l'intérieur la table pour reconnaître les taxons, que l'on télécharger 
# à l'adresse suivante:
# ----(MAJ AVEC LA DERNIERE VERSION)------
# https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz?download=1


# Charger le package {dada2}
############################
library(dada2)

# Config de l'environnement
###########################

# Lien vers le repértoire avec les fichiers fastq
################# A MODIFIER ####################
path <- "Limagne/Data/E969-Reads-forward-reverse"
#################################################

# Création de deux vecteurs avec les noms de fichiers CCS 
# (un avec les CCS fwd, l'autre avec les CCS rev) 
fastqFs <- sort(list.files(path, pattern=".fwd.fastq.gz"))
fastqRs <- sort(list.files(path, pattern=".rev.fastq.gz"))

# Spécifier le chemin des dossiers où stocker 
# les CCS filtrés après application de filterAndTrim()
filtpathF <- file.path(path, "forward_filtered") 
filtpathR <- file.path(path, "reverse_filtered") 

# Extrait les noms des échantillons pour fwd et rev
sample.names <- sapply(strsplit(basename(fastqFs), "[.]"), `[`, 2) 
sample.namesR <- sapply(strsplit(basename(fastqRs), "[.]"), `[`, 2) 
# Les deux vecteurs avec les noms des échantillons doivent être identiques sinon stop
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")

# Nettoyage des échantillons
############################

# Filtre et trie les échantillons selon taille et qualité des reads
track.filt <- filterAndTrim(
  # Fwd reads
  fwd=file.path(path, fastqFs), filt=file.path(filtpathF, fastqFs),
  # Rev reads
  rev=file.path(path, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
  #minQ=0, maxEE=3,
  minLen=800, maxLen=1600, 
  #truncLen=c(1800,1800),
  maxN=0, rm.phix=FALSE, 
  compress=TRUE, verbose=TRUE, multithread=TRUE
)

# Détermination du taux d'erreur
################################

# Création de deux vecteurs avec les noms de fichiers CCS filtrés
# (un avec les CCS fwd, l'autre avec les CCS rev)
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)

# On vérifie que l'on a le même nombre de fichiers pour fwd et rev
if(length(filtFs) != length(filtRs)) stop("Forward and reverse files do not match.")

# On définit les noms des échantillons pour les vecteurs nouvellement créés
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Fixe la seed
set.seed(100)
# Détermination du taux d'erreur pour les fwd
errF <- learnErrors(filtFs, nbases=1e8, errorEstimationFunction=PacBioErrfun, BAND_SIZE=32, multithread=TRUE) 
# Détermination du taux d'erreur pour les rev
errR <- learnErrors(filtRs, nbases=1e8, errorEstimationFunction=PacBioErrfun, BAND_SIZE=32, multithread=TRUE)

# Sample inference and merger of paired-end reads
#################################################

# Création du vecteur pour stocker les fichiers finaux
mergers <- vector("list", length(sample.names))

# Définition des noms des échantillons dans ce vecteur
names(mergers) <- sample.names

# Correction des CCS d'après les taux d'erreur estimés à l'étape précédente
# et fusion des CCS des 2 sens
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}

# Construct sequence table 
seqtab <- makeSequenceTable(mergers)
# Enlever les chimères
seqtab <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
# Sauvegarder la table
saveRDS(seqtab, file.path(path, "output/seqtab.rds"))

# Assignement de la taxonomie
#############################
# -------- Nom de la dernière version de la base SILVA a été modifié ici ----------------
tax <- assignTaxonomy(seqtab, file.path(path, "tax/silva_nr99_v138.1_train_set.fa.gz"), multithread=TRUE, tryRC=TRUE)
# Sauvegarde de la taxonomie
saveRDS(tax, file.path(path, "output/taxonomy.rds")) 

# Track reads through the pipeline
##################################
# Méthode inspirée de :
# https://benjjneb.github.io/dada2/tutorial.html
getNreads <- function(x) sum(getUniques(x))

track<-cbind(track.filt,sapply(mergers, getNreads),rowSums(seqtab))
colnames(track)<-c("input","filtered","merged","nonchim")

write.csv(track,file.path(path, "output/tracks_read.csv"))

# Export error rate
###################
plotErrors(errF, nominalQ=TRUE)
ggplot2::ggsave(file.path(path, "output/error_rate.pdf"))
