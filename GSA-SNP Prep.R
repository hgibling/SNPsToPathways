##### Preparing Plink Output for GSA-SNP Software #####

### Load sample dataset

assoc.data <- read.table("~/Desktop/gwas.assoc", header=T)
permuted.pvals <- read.csv("~/Desktop/perm_pvalues.txt")


### Combine SNP IDs with P values

gsasnp.input <- data.frame(SNP=assoc.data$SNP, P=assoc.data$P, permuted.pvals[,-1], P.100=permuted.pvals[,1])

# first permuted column renamed 'P.100' so avoid having two columns named 'P'


### Print file

write.table(gsasnp.input, "~/Desktop/GSA-SNP Input.txt", col.names=T, row.names=F, sep="\t", quote=F)