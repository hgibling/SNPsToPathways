##### Getting Pathway p-values using aSPUpath or HYST #####

library(aSPU)
library(gage)
library(biomaRt)
library(topGO)
library(GSA)
library(impute)


### Load sample dataset

assoc.data <- read.table("~/Desktop/gwas.assoc", header=T)


# SNP data needed from Plink:
	# ID (rs number)
	# Chromosome location
	# Base pair location
	# P value
	# Correlation data


### Generate gene set collections

# KEGG

get.kegg.gs <- kegg.gsets(species="hsa", id.type="entrez")
kegg.gs <- get.kegg.gs$kg.sets


# GO Biological Pathways
gobp.gs <- annFUN.org("BP", mapping="org.Hs.eg.db", ID="entrez")


# Import predefined gene sets from Bader Lab

get.bader.gs <- GSA.read.gmt("~/Desktop/Human_AllPathways_January_28_2015_symbol.gmt")


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
		all.gene.info <- getBM(attributes=c("hgnc_symbol", "chromosome_name", 		"start_position", "end_position"), filters=c("chromosome_name", 		"hgnc_symbol"), values=list(chromosome_name=c(1:22, "X", "Y"), 					hgnc_symbol=c(master.list)), mart=mart)
	} else if (type=="ID") {
		all.gene.info <- getBM(attributes=c("entrezgene", "chromosome_name", 		"start_position", "end_position"), filters=c("chromosome_name", 		"entrezgene"), values=list(chromosome_name=c(1:22, "X", "Y"), 					entrezgene=c(master.list)), mart=mart)
	} else {
		stop("Type of gene identification must be either 'ID' (default) or 		'symbol'.")
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


### Extract SNP info from dataset

snp.info <- data.frame(SNP=assoc.data[,2], Chrom=assoc.data[,1], Position=assoc.data[,3])


### Generate linkage disequilibrium matrix (SNP correlation matrix)

control.data <- read.table("~/Desktop/GO_Quad_DATA-clean-CEU.traw", stringsAsFactors=F, header=T)

rownames(control.data) <- control.data$SNP


# Remove uneccesary columns

control.genotypes <- control.data[,-(1:6)]


# Impute missing values

control.imputed <- impute.knn(as.matrix(control.genotypes))
control.rounded <- round(control.imputed$data)


# Transpose so columns are SNPs

control.trans <- t(control.rounded)

ld.matrix <- cor(control.trans)


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
				results <- aSPUsPath(assoc.data$P, 	# P values of SNPs
					corrSNP=ld.matrix,				# correlation of SNPs
					snp.info=snp.info,				# SNP location info
					gene.info=gene.info,			# gene location info
					Ps=T)							# P values instead of Z scores
				results.df[i,1] <- names(gs[i])
				results.df[i,2] <- results[length(results)] #aSPUpath is last
			}
		}
	} else if (method=="HYST") {
		for (i in 1:length(gs)) {
			if (length(gs[[i]]) > min & length(gs[[i]]) < max) {
				gene.info <- genes.in.gs(i, gs, all.gene.info)
				results <- Hyst(assoc.data$P,
					ldmatrix=ld.matrix,
					snp.info=snp.info,
					gene.info=gene.info)
				results.df[i,1] <- names(gs[i])
				results.df[i,2] <- results[21]
			}
		}
	} else {
		stop("Must indicate if method of analysis is 'aSPUpath' or 'HYST'.")
	}
	return(results.df)
}


# Run aSPUpath and order results by significance

kegg.aSPUpath <- run.snp.gsa(collection="kegg", method="aSPUpath")
kegg.aSPUpath.sig <- kegg.aSPUpath[order(kegg.aSPUpath$Pval)]

gobp.aSPUpath <- run.snp.gsa(collection="gobp", method="aSPUpath")
gobp.aSPUpath.sig <- gobp.aSPUpath[order(gobp.aSPUpath$Pval)]

bader.aSPUpath <- run.snp.gsa(collection="bader", method="aSPUpath")
bader.aSPUpath.sig <- bader.aSPUpath[order(bader.aSPUpath$Pval)]


# Run HYST and order results by significance

kegg.hyst <- run.snp.gsa(collection="kegg", method="HYST")
kegg.hyst.sig <- kegg.hyst[order(kegg.hyst$Pval)]

gobp.hyst <- run.snp.gsa(collection="gobp", method="HYST")
gobp.hyst.sig <- gobp.hyst[order(gobp.hyst$Pval)]

bader.hyst <- run.snp.gsa(collection="bader", method="HYST")
bader.hyst.sig <- bader.hyst[order(bader.hyst$Pval)]


### aSPUpath for raw genotype data (not p values)

phenotypes <- c(rep(0, 144), rep(1, 145))	# sample for now

run.aspupath.raw <- function(collection) {
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
	for (i in 1:length(gs)) {
		gene.info <- genes.in.gs(i, gs, all.gene.info)
		results <- aSPUpath(phenotypes, 	# phenotype data
			control.trans,					# genotype data
			snp.info=snp.info,				# SNP location info	
			gene.info=gene.info)			# gene location info			
		results.df[i,1] <- names(gs[i])
		results.df[i,2] <- results[length(results)] #aSPUpath is last
		}
	return(results.df)
}

test <- run.aspupath.raw(collection="kegg")
	
