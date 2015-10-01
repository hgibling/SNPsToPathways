##### Getting Pathway p-value using aSPUpath and HYST #####

library(aSPU)
library(gage)
library(biomaRt)

# SNP data needed from Plink:
	# ID (rs number)
	# Chromosome location
	# Base pair location
	# P value
	# Correlation data


# generate KEGG gene sets

get.gs.kegg <- kegg.gsets(species="hsa", id.type="entrez")
gs.kegg <- get.gs.kegg$kg.sets


# generate master list of genes in KEGG pathways

gene.list <- c()

for (i in 1:length(gs.kegg)) {
	genes <- as.numeric(gs.kegg[[i]])
	gene.list <- c(gene.list, genes)
}

unique.genes <- unique(gene.list)


# get chromosome location for the master list of genes

mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

all.gene.info <- getBM(attributes=c("entrezgene", "chromosome_name", "start_position", "end_position"), filters=c("chromosome_name", "entrezgene"), values=list(chromosome_name=c(1:22, "X", "Y"), entrezgene=c(unique.genes)), mart=mart)

## 108 genes are not in this list--are not mapped to a chromosome in biomart
## some are microRNAs, others unlocalized scaffolds
## some are alternate reference loci?


# extract gene info for one gene set

genes.in.gs <- function(gene.set, master.list) {
	set.genes <- as.numeric(master.list[[gene.set]])
	master.position <- all.gene.info[,1] %in% set.genes
	return(all.gene.info[master.position,])
}


# load sample dataset

assoc.data <- read.table("~/Desktop/gwas.assoc", header=T)


# extract snp info

snp.info <- data.frame(SNP=assoc.data[,2], Chrom=assoc.data[,1], Position=assoc.data[,3])


### Run aSPUpath for all gene sets

aSPUpath.results <- data.frame(Pathway=NA, Pval=NA)

for (i in 1:length(gs.kegg)) {
	gene.info <- genes.in.gs(i, gs.kegg)
	ldmatrix <- 
	results <- aSPUsPath(assoc.data$P, 	# P values of SNPs
		corrSNP=ldmatrix, 				# correlation of SNPs to controls
		pow=c(1, 2, 4, 8, Inf), 		# SNP gamma values
		pow2=c(1, 2, 4, 8), 			# gene gamma values
		snp.info=snp.info, 				# SNP location info
		gene.info=gene.info, 			# gene location info
		n.perm=1000, 					# 1000 permutations
		Ps=T)							# using P values instead of Z scores
	aSPUpath.results[i,1] <- names(gs.kegg[i])
	aSPUpath.results[i,2] <- results[21]
}


# order pathways/gene sets by significance

aSPUpath.sig <- aSPUpath.results[order(aSPUpath.results$Pval)]


### Run HYST for all gene sets

hyst.resuts <- data.frame(Pathway=NA, Pval=NA)

for (i in 1:length(gs.kegg)) {
	gene.info <- genes.in.gs(i, gs.kegg)
	ldmatrix <- 
	results <- Hyst(assoc.data$P, 		# P values of SNPs
		ldmatrix=ldmatrix, 				# correlation of SNPs to controls
		snp.info=snp.info, 				# SNP location info
		gene.info=gene.info) 			# gene location info
	hyst.resuts[i,1] <- names(gs.kegg[i])
	hyst.resuts[i,2] <- results[21]
}


# order pathways/gene sets by significance

hyst.sig <- hyst.results[order(hyst.results$Pval)]
