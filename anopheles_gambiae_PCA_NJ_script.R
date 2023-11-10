######################## PCA FOR GENES OF INTEREST #########################

## Make PCA for just genes of interest
## Use bcftools to subset the filtered vcf just for the genes of interest and make a PCA of that

# VGSC AgamP4_2L:2358158 to 2431617
# subset VCF

bcftools view -r 2L:2358158-2431617 F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz -Oz -o VGSC_only_F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz


## Shift + Options + R to open a new terminal
## Ctrl + Option + Enter to send commands to terminal 
## ssh sophie@s11
## conda activate fastq2matrix
## R to start working in R
## save your R script using nano as scriptname.R
## if you are sure that your script will run smoothly, you can just type Rscript <filename>.R
## otherwise, you will need to run line by line when you have activated R to make sure it works

## install required packages

install.packages("dplyr")
install.packages("scales")
install.packages("data.table")

## then attach all of these packages by using library(packagename)
library(dplyr)
library(scales)
library(data.table) 

##setwd
setwd('/mnt/storage11/sophie/bijagos_mosq_wgs/2022_gambiae_fq2vcf/gambiae_nov2022_genomicdb/gambiae_nov2022_genotypedvcf/gambiae_nov2022_combinedvcf_filteringsteps')
getwd()

## matrices
#F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.mat.bin
#F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.noniupac.mat

## load binary snp matrix (biallelic nucleotide matrix for making the PCA plots, ending mat.bin)
snp<-fread("VGSC_only_F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.mat.bin",sep="\t",header=T)
snp<-as.data.frame(snp)
count(snp)
nrow(snp)
ncol(snp)

## read in metadata

met <- read.table("metadata_gambiae_2022.csv", sep=",", stringsAsFactors = FALSE, header = TRUE)
count(met)

## Note:For one of my metadata sets I had to use sep=' ' to get R to recognise spaces, and sometimes it has to be sep="\t"

## check that metadata has loaded correctly
head(met,15)
met<-met[match(colnames(snp)[-(1:3)], met$sample),]
count(met)

## allow metadata fields for plotting (aka plot by country, region, site, sample, etc.)

met <- met %>% select(c(sample,year,country,island,phenotype))

## inspect metadata
colnames(met)
table(met$"country")
unique(met[c("island")])
unique(met[c("phenotype")])
count(met)

##Select which countries to include in analysis

#met_subset<-met[which(met$country %in% c("Ivory_Coast","Guinea-Bissau","Gambia","Uganda","Kenya","Gabon","Eritrea","DRC","Madagascar","Malawi")),]
#met_sample<-met_subset$sample_id
#snp_subset <- snp[,which(colnames(snp) %in% c("chr", "pos", "ref", met_sample))]

#count(snp_subset)

#Inspect final metadata

#colnames(met_subset)
#table(met_subset$"country")
#unique(met_subset[c("country")])
#unique(met_subset[c("region")])
#count(met_subset)

##Remove outliers:

#met_eth<-met_subset[which(met_subset$sample_id=="ERR012269"),]
#met_eth<-met_subset$sample_id[which(met_subset$sample_id=="ERR012269")]
#met_eth<-as.vector(met_eth)
#snp_subset<-snp_subset[,-which(colnames(snp_subset)%in%met_eth)]
#met_subset<-met_subset[-which(met_subset$sample_id=="ERR012269"),]
#count(met_subset)

#met_eth<-met_subset[which(met_subset$sample_id=="ERR484640"),]
#met_eth<-met_subset$sample_id[which(met_subset$sample_id=="ERR484640")]
#met_eth<-as.vector(met_eth)
#snp_subset<-snp_subset[,-which(colnames(snp_subset)%in%met_eth)]
#met_subset<-met_subset[-which(met_subset$sample_id=="ERR484640"),]
#count(met_subset)

# can write a new metadata file once you've removed samples
#write.csv(met_subset, file = paste("2021_03_18_pf_metadata_db_all_w_cov_qc_fws_np_plus_bijagos_met_subset_african_regions.csv"))

## CHOOSE PLOT TYPE (you MUST do this before you make the distance matrix!)
#run by site
#met_subset$site<-factor(met_subset$site)
#all(met_subset$sample_id==colnames(snp_subset[-c(1:3)]))
#snp_c<-snp_subset[,-(1:3)]

#run by phenotype
met$phenotype<-factor(met$phenotype)
all(met$sample==colnames(snp[-c(1:3)])) 
snp_c<-snp[,-(1:3)] #this overlays your SNP variable with the binary SNP matrix, and removes anything that is not in both

install.packages("amap")
library(amap)

# calculate the distance matrix used to make the plots

#dist_m<-Dist(t(snp_c), method = "manhattan")
#summary(dist_m)

#print("hello")

## The above line of code did not work for me because the dataset was too large.
## You can create distance matrix using plink on the server instead (not in R).
## Then load it in and combine with slimmed metadata.

## load the distance matrix
library(dplyr)

dist <- read.table("VGSC_only_dist_m.dist", header = F)
id <- read.table("VGSC_only_dist_m.dist.id", header = F)

ncol(dist)
ncol(id)

names(id)
names(met)

#samplestoinclude <- which(id[,1] %in% met$sample) #takes sample IDs from met_subset and matches with id variable
#length(samplestoinclude) #prints the variable so you can check it's the right number of samples to be in the distance matrix
#dist_filt <- dist[samplestoinclude,samplestoinclude] #filters the distance matrix to only include the coordinates which contain distances for these samples
#dim(dist_filt) #check that it contains the same numbers of rows and columns, and that the number of distances = the number of samples

#id_filt <-id[samplestoinclude,] #filters the id variable to only include the sample id's we care about
#dim(id_filt) #check that it's the right number of samples
#desc <- id_filt %>% left_join(met, by = c("V1" = "sample_id")) #desc is a function that joins up the id_filt variable with the rest of the metadata information for the samples
#head(desc) #check and see that the metadata has joined

##read in new metadata file here, the one that has been created after slimming down the samples
##need to subset id variable so that it only has the sample_ids that are in met_subset

#dist_m <- as.matrix(dist_filt)
#colnames(dist_m) <- desc$V1
#rownames(dist_m) <- desc$V1
#dim(dist_m)

## save dist_m to start from next time, as this takes a while to produce

#saveRDS(dist_m, file = "/mnt/storage9/sophie/basebijagospopgen/wgs_global_database/filtered/pca/distancematrix_20222311.rds")

## to read this RDS file in future to start from this point, make sure you're in the right directory and then use this code:
# readRDS(distancematrix.rds, refhook = NULL)


## Now you have made the distance matrix, calculate components with cmdscale

cmd_m<-cmdscale(dist, k = 10, eig = TRUE, x.ret = TRUE)

# Calculate the proportion of variance
# include the first and second values in the axes titles, first value = PC1, second value = PC2, etc...

#vars <- round(cmd_pv$eig / sum(cmd_pv$eig)*100, 1)
variance <- (cmd_m$eig / sum(cmd_m$eig)*100)
var1<-variance[[1]]
var2<-variance[[2]]

island<-met$island
island<-as.character(island)

phenotype<-met$phenotype
phenotype<-as.character(phenotype)

year<-met$year
year<-as.character(year)

sample<-met$sample
sample<-as.character(sample)

# colour palette

install.packages("RColorBrewer")
install.packages("colorRamps")
library(RColorBrewer)
library(colorRamps)
library(scales)
library(unikn)

pal<-c("purple", "red", "green")

# generate plots 

# PCA by phenotype
# PCA by island with sample names
pdf("VGSC_gambiae_2022_plot_by_phenotype.pdf", width = 8)
par(mar = c(5.1, 4.1, 4.1, 11.1), xpd = TRUE)
plot(
  cmd_m$points[, 1], cmd_m$points[, 2],
  col = pal[as.numeric(as.factor(phenotype))], pch = 20, cex = 0.9,
  xlab = paste("PC1 ", round(var1, digits = 2), "%"),
  ylab = paste("PC2 ", round(var2, digits = 2), "%")
)

# Add sample names as labels
text(
 cmd_m$points[, 1], cmd_m$points[, 2], labels = sample,
pos = 3, col = "black", cex = 0.7
)

legend("bottomright", inset = c(-0.43, 0), col = unique(pal[as.factor(phenotype)]), pch = 20, leg = unique(as.factor(phenotype)))
dev.off()



# generate a tree file that corresponds to the PCA plots, can be visualised in iTOL or in R
install.packages("ape")
library(ape)
tree_m<-nj(dist_m)

# Writing to a newick tree file
write.tree(phy=tree_m, file="/mnt/storage9/sophie/basebijagospopgen/database/filtered/pca/2022_02_21_2_tree.newick")


#### A.Miles thesis methodology for PCA:
# SNPs for inclusion in principal components analysis were chosen by selecting biallelic
# variants from within the regions 3R:1-37 Mbp and 3L:15-41 Mbp. Only variants with
# minor allele frequency >= 1% were retained and each chromosome arm was randomly
# down-sampled to 100,000 variants. I then pruned to remove SNPs in linkage disequilibrium,
# excluding SNPs above an r2 threshold of 0.01 in moving windows of 500 SNPs with a step
# size of 250 SNPs. SNPs from both chromosome arms were then concatenated, and PCA
# was run following methods described in Patterson et al. (2006).