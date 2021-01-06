require(QSutils, quietly=TRUE)
library(seqinr)
library(QSutils)
library(parallel)
library(unix)
#install.packages("seqinr")
#BiocManager::install("GenomicAlignments")
library(GenomicAlignments)
#rlimit_as(1e12)

main <- function(){
	args <- commandArgs(trailingOnly = TRUE)
	base <- args[1] 
	gene <- args[2]
	bampath <- args[3] #"A2989_S43-V2/A_PB1.bam"
	fapath <- args[4] #"A2989_S43-V2/A_PB1.fasta"
	fa1 <- read.fasta(fapath)
	mybam <- BamFile(bampath)
	mygenelen <- getLength(fa1)
	mygene <- gene
	region1 <- GRanges(mygene , IRanges(1, mygenelen))
	#?stackStringsFromBam
	hap <- stackStringsFromBam(mybam, param = region1, Lpadding.letter = "-", Rpadding.letter = "-")
	hap
	lstCollapsed <- Collapse(hap)

	rm(hap, mybam)
	cons <- ConsSeq(lstCollapsed$hseqs)

	getLength(cons)
	lstCorrected <- CorrectGapsAndNs(lstCollapsed$hseqs[1:length(lstCollapsed$hseqs)], cons)
	rm(cons)
	#lstCorrected
	lstRecollapsed <- Recollapse(lstCorrected, lstCollapsed$nr)

	thr <- 0.01
	sz <- sum(lstRecollapsed$nr)
	#sz

	nr.pop <- lstRecollapsed$nr
	#nr.pop

	(fl <- nr.pop>=sz*thr/100)

	nr.sz <- nr.pop[fl]
	#nr.sz

	shannon_after_thr <- Shannon(nr.sz)

	#sz2 <- 50000
	#nr.sz
	#fl <- DSFT(nr.sz, sz2)
	#fl

	#length(lstRecollapsed$nr)

	rm(lstCollapsed, lstCollapsed)
	dst <- DNA.dist(lstRecollapsed$seqs, model = "raw")
	
	ND <- NucleotideDiversity(dst)
	ND <- as.numeric(ND)	
	Mfe <- as.numeric(MutationFreq(dst))
	dataplus <- data.frame(base, gene, sz, shannon_after_thr, ND, Mfe)
	smpsize <- data.frame(base, sz)
	write.table(dataplus, file="./DNA_diversity.csv", append = T, sep = "\t", row.names = F, col.names = F) 
	write.table(smpsize, file=paste0("./",gene,".tsv"), append = T, sep = "\t", row.names = F, col.names = F)
	rm(list = ls())
}
main()
	
