##### Preparing Plink Output for GSA-SNP Software #####

# load sample dataset

assoc.data <- read.table("~/Desktop/gwas.assoc", header=T)
permuted.pvals <- read.csv("~/Desktop/perm_pvalues.txt")


# combine SNP IDs with P values

gsasnp.input <- data.frame(SNP=assoc.data$SNP, P=assoc.data$P, P.0=permuted.pvals$P, permuted.pvals[,-1])


# print file

write.table(gsasnp.input, "~/Desktop/GSA-SNP Input.txt", col.names=T, row.names=F, sep="\t", quote=F)