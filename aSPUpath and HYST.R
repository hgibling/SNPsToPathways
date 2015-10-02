##### Getting Pathway p-values using aSPUpath or HYST #####

library(aSPU)
library(gage)
library(biomaRt)
library(topGO)


# load sample dataset

assoc.data <- read.table("~/Desktop/gwas.assoc", header=T)


# SNP data needed from Plink:
	# ID (rs number)
	# Chromosome location
	# Base pair location
	# P value
	# Correlation data


# generate KEGG gene sets

get.gs.kegg <- kegg.gsets(species="hsa", id.type="entrez")
gs.kegg <- get.gs.kegg$kg.sets


# generate GO Biological Pathways gene sets

gs.gobp <- annFUN.org("BP", mapping="org.Hs.eg.db", ID="entrez")


# generate master list of genes in KEGG or GO BP gene sets

get.master.list <- function(gene.set.list) {
	gene.list <- c()
	for (i in 1:length(gene.set.list)) {
		genes <- as.numeric(gene.set.list[[i]])
		gene.list <- c(gene.list, genes)
		}
	unique.genes <- unique(gene.list)
	return(unique.genes)
}

kegg.master <- get.master.list(gs.kegg)
gobp.master <- get.master.list(gs.gobp)


# get chromosome locations for the KEGG or GO BP master list of genes

mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

get.gene.info <- function (master.list) {
	all.gene.info <- getBM(attributes=c("entrezgene", "chromosome_name", 	"start_position", "end_position"), filters=c("chromosome_name", 	"entrezgene"), values=list(chromosome_name=c(1:22, "X", "Y"), 				entrezgene=c(master.list)), mart=mart)
}

kegg.gene.info <- get.gene.info(kegg.master)
gobp.gene.info <- get.gene.info(gobp.master)

## some genes are not in these lists--are not mapped to a chromosome in biomart
## some are microRNAs, others unlocalized scaffolds
## some are alternate reference loci?


# Extract gene info for one gene set

genes.in.gs <- function(gene.set.position, gene.set.list, master.gene.info) {
	set.genes <- as.numeric(gene.set.list[[gene.set.position]])
	master.position <- master.gene.info[,1] %in% set.genes
	return(master.gene.info[master.position,])
}


# Extract SNP info from dataset

snp.info <- data.frame(SNP=assoc.data[,2], Chrom=assoc.data[,1], Position=assoc.data[,3])


# Generate linkage disequilibrium matrix (SNP correlation matrix)

ld.matrix <-


### Run aSPUpath

# For KEGG gene sets

aSPUpath.kegg <- data.frame(Pathway=NA, Pval=NA)

for (i in 1:length(gs.kegg)) {
	gene.info <- genes.in.gs(i, gs.kegg, kegg.gene.info)
	results.kegg.a <- aSPUsPath(assoc.data$P, 	# P values of SNPs
		corrSNP=ld.matrix, 				# correlation of SNPs to controls
		pow=c(1, 2, 4, 8, Inf), 		# SNP gamma values
		pow2=c(1, 2, 4, 8), 			# gene gamma values
		snp.info=snp.info, 				# SNP location info
		gene.info=gene.info, 			# gene location info
		n.perm=1000, 					# 1000 permutations
		Ps=T)							# using P values instead of Z scores
	aSPUpath.kegg[i,1] <- names(gs.kegg[i])
	aSPUpath.kegg[i,2] <- results.kegg.a[21]
}

aSPUpath.kegg.sig <- aSPUpath.kegg[order(aSPUpath.kegg$Pval)]


# For GO BP gene sets

aSPUpath.gobp <- data.frame(Pathway=NA, Pval=NA)

for (i in 1:length(gs.gobp)) {
	gene.info <- genes.in.gs(i, gs.gobp, gobp.gene.info)
	results.gobp.a <- aSPUsPath(assoc.data$P, 	# P values of SNPs
		corrSNP=ld.matrix, 				# correlation of SNPs to controls
		pow=c(1, 2, 4, 8, Inf), 		# SNP gamma values
		pow2=c(1, 2, 4, 8), 			# gene gamma values
		snp.info=snp.info, 				# SNP location info
		gene.info=gene.info, 			# gene location info
		n.perm=1000, 					# 1000 permutations
		Ps=T)							# using P values instead of Z scores
	aSPUpath.gobp[i,1] <- names(gs.gobp[i])
	aSPUpath.gobp[i,2] <- results.gobp.a[21]
}

aSPUpath.gobp.sig <- aSPUpath.gobp[order(aSPUpath.gobp$Pval)]


### Run HYST

# For KEGG gene sets

hyst.kegg <- data.frame(Pathway=NA, Pval=NA)

for (i in 1:length(gs.kegg)) {
	gene.info <- genes.in.gs(i, gs.kegg, kegg.gene.info)
	results.kegg.h <- Hyst(assoc.data$P, 	# P values of SNPs
		ldmatrix=ld.matrix, 			# correlation of SNPs to controls
		snp.info=snp.info, 				# SNP location info
		gene.info=gene.info) 			# gene location info
	hyst.kegg[i,1] <- names(gs.kegg[i])
	hyst.kegg[i,2] <- results.kegg.h[21]
}

hyst.kegg.sig <- hyst.kegg[order(hyst.kegg$Pval)]


# For GO BP gene sets

hyst.gobp <- data.frame(Pathway=NA, Pval=NA)

for (i in 1:length(gs.kegg)) {
	gene.info <- genes.in.gs(i, gs.gobp, gobp.gene.info)
	results.gobp.h <- Hyst(assoc.data$P, 	# P values of SNPs
		ldmatrix=ld.matrix, 			# correlation of SNPs to controls
		snp.info=snp.info, 				# SNP location info
		gene.info=gene.info) 			# gene location info
	hyst.gobp[i,1] <- names(gs.kegg[i])
	hyst.gobp[i,2] <- results.gobp.h[21]
}

hyst.gobp.sig <- hyst.gobp[order(hyst.gobp$Pval)]

