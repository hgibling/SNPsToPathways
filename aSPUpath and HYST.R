##### Getting Pathway p-value using aSPUpath and HYST #####

library(aSPU)


#generate simulated dataset

sim <- simPathAR1Snp(nGenes=20, 	# 20 genes in gene set
	nGenes1=1, 
	nSNPlim=c(1,20), 				# between 1 and 20 SNPs per gene
	nSNP0=1:3, 
	LOR=0.2, 
	rholim=c(0,0), 
	n=100, 
	MAFlim=c(0.05, 0.4), 
	p0=0.05)


### HYST ###

# get p-values for each SNP from simulated dataset

pvals <- getlogitp(sim$Y, sim$X)



# get correlation matrix for SNPs using controls (to account for linkage disequilibrium)

ld <- cor(sim$X[sim$Y==0,])


# run HYST algorithm

hyst.result <- Hyst(pvec=logitp, ldmatrix=ld, snp.info=sim$snp.info, gene.info=sim$gene.info)

# output is a single p-value for the gene set tested


### aSPUpath ###

aspupath.result <- aSPUpath(sim$Y,
	sim$X,
	snp.info=sim$snp.info,
	gene.info=sim$gene.info,
	model="binomial",			# binomial because using binary traits
	pow=1:8,					# gamma values fpr SNP power
	pow2=c(1, 2, 4, 8),			# gamma values for gene power
	n.perm=1000)				# number of permutations
	
# output has all p-values for SNP and gene power combinations, plus an aSPU p-value (from the adaptive model) for the gene set