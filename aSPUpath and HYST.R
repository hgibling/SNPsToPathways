##### Getting Pathway p-values using aSPUpath or HYST on pre-determined p values #####

library(aSPU)
library(gage)
library(biomaRt)
library(GSA)
library(impute)
library(dplyr)


### Load sample dataset

assoc.data <- read.table("gwas.assoc", header=T)


# SNP data needed from Plink:
	# ID (rs number)
	# Chromosome location
	# Base pair location
	# P value
	# Control genotypes


### Generate gene set collections

# KEGG

get.kegg.gs <- kegg.gsets(species="hsa", id.type="entrez")
kegg.gs <- get.kegg.gs$kg.sets


# GO Biological Pathways

get.gobp.gs <- GSA.read.gmt("Homo_sapiens_GSEA_GO_sets_bp_ids_highquality_April_2015.gmt")
gobp.gs <- get.gobp.gs$genesets
names(gobp.gs) <- get.gobp.gs$geneset.names


# Import predefined gene sets from Bader Lab

get.bader.gs <- GSA.read.gmt("Human_AllPathways_January_28_2015_symbol.gmt")


# Remove blank entry at end of each Bader Lab gene set

for (i in 1:length(get.bader.gs$genesets)) {
	get.bader.gs$genesets[[i]] <- get.bader.gs$genesets[[i]][-length(get.bader.gs	$genesets[[i]])]
}

bader.gs <- get.bader.gs$genesets
names(bader.gs) <- get.bader.gs$geneset.names


### Generate master list of genes in each gene set collection

get.master.list <- function(gene.set.list) {
	gene.list <- c()
	for (i in 1:length(gene.set.list)) {
		genes <- gene.set.list[[i]]
		gene.list <- c(gene.list, genes)
		}
	unique.genes <- unique(gene.list)
	return(unique.genes)
}

kegg.master <- as.numeric(get.master.list(kegg.gs))
gobp.master <- as.numeric(get.master.list(gobp.gs))
bader.master <- get.master.list(bader.gs)	# gene symbols, so not numeric


### Get chromosome locations for the master list of genes for each collection

mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

get.gene.info <- function (master.list, type="ID") {
	if (type=="symbol") {
		all.gene.info <- getBM(attributes=c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
		filters=c("chromosome_name", "hgnc_symbol"),
		values=list(chromosome_name=c(1:22, "X", "Y"),
		hgnc_symbol=c(master.list)), mart=mart)
	} else if (type=="ID") {
		all.gene.info <- getBM(attributes=c("entrezgene", "chromosome_name", "start_position", "end_position"),
		filters=c("chromosome_name", "entrezgene"),
		values=list(chromosome_name=c(1:22, "X", "Y"),
		entrezgene=c(master.list)), mart=mart)
	} else {
		stop("Type of gene identification must be either 'ID' (default) or 'symbol'.")
	}
}

kegg.gene.info <- get.gene.info(kegg.master)
gobp.gene.info <- get.gene.info(gobp.master)
bader.gene.info <- get.gene.info(bader.master, type="symbol")

# some genes are not in these lists--are not mapped to a chromosome in biomart
# some are microRNAs, others unlocalized scaffolds
# some are alternate reference loci?


### Extract gene info for one gene set

genes.in.gs <- function(gene.set.position, gene.set.list, master.gene.info) {
	set.genes <- gene.set.list[[gene.set.position]]
	master.position <- master.gene.info[,1] %in% set.genes
	return(master.gene.info[master.position,])
}


### Extract SNP info from dataset (MAF>5% in controls)

rare.snps <- which(assoc.data$F_U < 0.05)
assoc.clean <- assoc.data[-rare.snps,]
all.snp.info <- data.frame(SNP=assoc.clean[,2], Chrom=assoc.clean[,1], Position=assoc.clean[,3])


### Find SNPs associated with genes in a gene set

snps.in.gs <- function(gene.info) {
	snp.info <- data.frame(NULL)
	for (i in 1:nrow(gene.info)) {
		snp.info.gene <- filter(all.snp.info, Chrom==gene.info[i,2]) %>%
		filter(Position > gene.info[i,3]-20000 & Position < gene.info[i,4]+20000)
		snp.info <- unique(rbind(snp.info, snp.info.gene))
	}
	snp.info <- snp.info[order(snp.info$SNP),]
	return(snp.info)
}

# gets SNPs positioned anywhere between 20kb downstream and 20kb upstream of a gene 
# can be adjusted if desired


### Prepare genotype data for generating linkage disequilibrium matrix (SNP correlation matrix)

control.data <- read.table("GO_Quad_DATA-clean-CEU.traw", stringsAsFactors=F, header=T)

rownames(control.data) <- control.data$SNP


# Remove uneccesary columns and rare SNPS

control.genotypes <- control.data[-rare.snps,-(1:6)]


# Impute missing values

control.imputed <- impute.knn(as.matrix(control.genotypes))
control.rounded <- round(control.imputed$data)


### Generate linkage disequilibrium matrix (SNP correlation matrix) for SNPs within a gene set

get.ld.matrix <- function(snp.info) {
	position <- which(rownames(control.rounded) %in% snp.info[,1]==T)
	control.snps <- control.rounded[position,]
	control.order <- control.snps[order(rownames(control.snps)),]
	control.trans <- t(control.order)
	ld.matrix <- cor(control.trans)
	if (anyNA(ld.matrix)==T) {
		pos <- which(is.na(ld.matrix[1,])==T)
		ld.matrix <- ld.matrix[-pos, -pos]		# remove NAs from LD matrix
	}
	return(ld.matrix)
}


### Get P values for only those SNPs within a gene set

get.p.values <- function(snp.info) {
	position <- which(assoc.clean$SNP %in% snp.info[,1]==T)
	snps <- assoc.clean[position,]
	snps <- snps[order(snps$SNP),]
	return(snps$P)
}


### Run aSPUpath or HYST for any of the three gene set colelctions

run.snp.gsa <- function(collection, method, min=10, max=300) {
	results.df <- data.frame(Pathway=NA, Pval=NA)
	if (collection=="kegg") {
		gs <- kegg.gs
		all.gene.info <- kegg.gene.info
	} else if (collection=="gobp") {
		gs <- gobp.gs
		all.gene.info <- gobp.gene.info
	} else if (collection=="bader") {
		gs <- bader.gs
		all.gene.info <- bader.gene.info
	} else {
		stop("Must indicate if collection is 'kegg', 'gobp', or 'bader'.")
	}
	if (method=="aSPUpath") {
		for (i in 1:length(gs)) {
			if (length(gs[[i]]) > min & length(gs[[i]]) < max) {
				gene.info <- genes.in.gs(i, gs, all.gene.info)
				snp.info <- snps.in.gs(gene.info)
				snp.pvals <- get.p.values(snp.info)
				ld.matrix <- get.ld.matrix(snp.info)
				results <- aSPUsPath(snp.pvals,
					corrSNP=ld.matrix,
					snp.info=snp.info,
					gene.info=gene.info,
					Ps=T)							# P values instead of Z scores
				results.df[i,1] <- names(gs[i])
				results.df[i,2] <- results[length(results)] #aSPUpath is last
			}
			if (i %% 10 == 0 ) {
				print(paste("Done set", i, "of", length(gs)))	# prints progress
			}
		}
	} else if (method=="HYST") {
		for (i in 1:length(gs)) {
			if (length(gs[[i]]) > min & length(gs[[i]]) < max) {
				gene.info <- genes.in.gs(i, gs, all.gene.info)
				snp.info <- snps.in.gs(gene.info)
				snp.pvals <- get.p.values(snp.info)
				ld.matrix <- get.ld.matrix(snp.info)
				results <- Hyst(snp.pvals,
					ldmatrix=ld.matrix,
					snp.info=snp.info,
					gene.info=gene.info)
				results.df[i,1] <- names(gs[i])
				results.df[i,2] <- results[length(results)]
			}
			if (i %% 10 == 0 ) {
				print(paste("Done set", i, "of", length(gs)))	# prints progress
			}
		}
	} else {
		stop("Must indicate if method of analysis is 'aSPUpath' or 'HYST'.")
	}
	results.df <- na.omit(results.df) %>%
	arrange(Pval)
	return(results.df)
}

# loops through each gene set in the chosen collection and runs either aSPUpath or HYST on the gene set, returning the (adapted) p-values
# MUST ADD CORRECTION FOR MULTIPLE HYPOTHESIS TESTING--not currently included


# Run aSPUpath on gene set collections

kegg.aSPUpath <- run.snp.gsa(collection="kegg", method="aSPUpath")
gobp.aSPUpath <- run.snp.gsa(collection="gobp", method="aSPUpath")
bader.aSPUpath <- run.snp.gsa(collection="bader", method="aSPUpath")


# Run HYST on gene set collections

kegg.hyst <- run.snp.gsa(collection="kegg", method="HYST")
gobp.hyst <- run.snp.gsa(collection="gobp", method="HYST")
bader.hyst <- run.snp.gsa(collection="bader", method="HYST")

