##### Getting Pathway p-value using aSPUpath and HYST #####

library(aSPU)
library(gage)
library(biomaRt)


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


# get chromosome location

mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

gene.info <- getBM(attributes=c("entrezgene", "chromosome_name", "start_position", "end_position"), filters=c("chromosome_name", "entrezgene"), values=list(chromosome_name=c(1:22, "X", "Y"), entrezgene=c(unique.genes)), mart=mart)

## 108 genes are not in this list--are not mapped to a chromosome in biomart
## some are microRNAs, others unlocalized scaffolds
## some are alternate reference loci?


# extract gene info for one gene set

genes.in.gs <- function(gene.set, master.list) {
	set.genes <- as.numeric(master.list[[gene.set]])
	master.position <- gene.info[,1] %in% set.genes
	return(gene.info[master.position,])
}

