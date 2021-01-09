library(ape)
library(seqinr)
library(QSutils)
library(GenomicAlignments)
library(Biostrings)
library(DECIPHER)
library(ips)
library(lmds)

args <- commandArgs(trailingOnly = TRUE)

mysample <- read.table(args[1])
parent <- args[1]
parent <- sub("_daughter.txt", "", parent)

#mysample <- read.table("asmatic_daughter.txt")
#parent <- "asthmatic"

for (id in as.vector(mysample$V1)){


  fa1 <- read.fasta(paste0(id, "-V2/A_PB1.fasta")) #daughter
  ref <- paste0(id, "-V2/A_PB1.fasta")


  fa2 <- read.fasta("../../fastq/A3073_S44-V2/A_PB1.fasta") #non-asthmatic parent
  #fa2 <- read.fasta("../../fastq/A6016/A_PB1.fasta") #asthmatic parent

  mygenelen1 <- getLength(fa1)

  mygene <- "A_PB1"

  mygenelen2 <- getLength(fa2)

  mybam1 <- BamFile(paste0(id,"-V2/A_PB1.bam")) #daughter samples

  #mybam2 <- BamFile("../../fastq/A6016_S25-V2/A_PB1.bam") #asthmatic parent
  mybam2 <- BamFile("../../fastq/A3073_S44-V2/A_PB1.bam") #non-athmatic

  region1 <- GRanges(mygene , IRanges(1, mygenelen1))
  region2 <- GRanges(mygene , IRanges(1, mygenelen2))

  hap1 <- stackStringsFromBam(mybam1, param = region1, Lpadding.letter = "-", Rpadding.letter = "-")
  hap1

  hap2 <- stackStringsFromBam(mybam2, param = region2, Lpadding.letter = "-", Rpadding.letter = "-")
  hap2

  lstCollapsed1 <- Collapse(hap1)
  rm(hap1, mybam1)

  lstCollapsed2 <- Collapse(hap2)
  rm(hap2, mybam2)

  cons1 <- ConsSeq(lstCollapsed1$hseqs)
  cons2 <- ConsSeq(lstCollapsed2$hseqs)

  lstCorrected1 <- CorrectGapsAndNs(lstCollapsed1$hseqs[1:length(lstCollapsed1$hseqs)], cons1)
  lstCorrected2 <- CorrectGapsAndNs(lstCollapsed2$hseqs[1:length(lstCollapsed2$hseqs)], cons2)

  rm(cons2)

  lstCorrected1
  lstCorrected2

  lstRecollapsed1 <- Recollapse(lstCorrected1, lstCollapsed1$nr)
  lstRecollapsed2 <- Recollapse(lstCorrected2, lstCollapsed2$nr)

  features1 <- c(sprintf("daughter%d", seq(1,length(lstRecollapsed1$seqs))))
  features2 <- c(sprintf("parent%d", seq(1,length(lstRecollapsed2$seqs))))

  (names(lstRecollapsed1$seqs) <- features1)
  (names(lstRecollapsed2$seqs) <- features2)

  lstMerged <- c(lstRecollapsed1$seqs, lstRecollapsed2$seqs)

  rm(lstRecollapsed2, lstRecollapsed1)

  saveRDS(lstMerged, file = paste0("lstMerged.RDS"))


  lstMerged <- readRDS("lstMerged.RDS")
  lstMerged <- chartr(".", "-", lstMerged)
  lstMerged <- RemoveGaps(lstMerged, "all", processors=NULL)
  writeXStringSet(lstMerged, 'fname.fa')
  seqs <- "fname.fa"

  exec <- "mafft"
  call.mafft <- paste(exec, "--auto --thread -1 --addfragments", seqs, ref, "> alnformds.fa")

  system(call.mafft, intern = FALSE, ignore.stdout = FALSE)

  aln <- "alnformds.fa"

  lstMerged <- readDNAStringSet(aln, format = "fasta")
  lstMerged <- lstMerged[-1,]
  ref <- ReadAmplSeqs(ref)
  lstMerged <- CorrectGapsAndNs(lstMerged[1:length(lstMerged)], cons1)

  dst <- DNA.dist(lstMerged, model = "raw")

  dst <- as.matrix(dst)

  rm(lstMerged)

  rds <- paste0(parent, "-", id, "-dst.RDS")

  saveRDS(dst, rds) # a matrix object

  mat <- readRDS(rds)

  dist_2lm <- select_landmarks(mat, num_landmarks = 12500)
  fit <- cmdscale_landmarks(dist_2lm, ndim = 2)
  x <- fit[,1]
  y <- fit[,2]
  df <- data.frame(x, y, row.names = row.names(fit))

  mydir <- "MDS_PDFs/"
  dir.create(mydir)

  saveRDS(df, paste0(mydir, parent,"-",id,"-dataframe.mds.RDS"))

  pdf(paste0(mydir,parent, id, ".MDS.pdf"))
  sp2 <- ggplot(df, aes(x, y)) +
    geom_point(colour = ifelse(grepl("daughter", row.names(df)), 'red', 'blue') ,
               size = 4,
               alpha = 0.2,
               shape = 12) +
    geom_polygon(data= hulls, alpha = 0.2) +
    #stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="black") +
    coord_cartesian(ylim=c(-0.43,-0.1)) + coord_cartesian(xlim=c(-0.15,0.2)) +
    labs(colour = "Level",
         x = "Component 1",
         y = "Component 2",
         tag = "A8323",
         title = "Landmark Multiple Dimension Scaling",
         subtitle = "non-asthmatic donor")


  sp2
  dev.off()
  rm(list = ls())
}
