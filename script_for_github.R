library(phyloseq)
library(DESeq2)
packageVersion(DESeq2)
library(ggplot2)
theme_set(theme_bw())
library(dplyr)
library(tidyr)

otutable <- import_biom(BIOMfilename = 'biom_for_phyloseq.biom', 
                        treefilename = 'rep_set_no_chimeras.tre', 
                        parseFunction = parse_taxonomy_greengenes)
#warnings()

mapping <- import_qiime_sample_data(mapfilename = 'map.txt')

phylo <- merge_phyloseq(otutable, mapping)

###################################
####Testing for effect of phenotype
###################################
phylo_subset <- subset_samples(phylo, Phenotype != "Dead")
sample_data(phylo_subset)$individual <- factor(sample_data(phylo_subset)$individual)

pheno <- phyloseq_to_deseq2(phylo_subset, ~ individual + Phenotype)

pheno_results <- DESeq(pheno, test="Wald")

pheno_res <- results(pheno_results)
summary(pheno_res)

alpha <- 0.05
pheno_sigtab <- pheno_res[which(pheno_res$padj < alpha), ]
pheno_sigtab <- cbind(as(pheno_sigtab, "data.frame"), as(tax_table(phylo)[rownames(pheno_sigtab), ], "matrix"))
write.table(pheno_sigtab, "pheno_padj05.txt", sep="\t")

##########################################
####Compare samples at day 0
##########################################
day0_samples<-prune_samples(sample_data(phylo)$Day==0, phylo)

day0_test <- phyloseq_to_deseq2(day0_samples, ~Final_phenotype)
day0_test_deseq <- DESeq(day0_test, test="Wald")
day0_results <- results(day0_test_deseq)
summary(day0_results)
head(day0_results)

alpha <- 0.05
day0_sigtab <- day0_results[which(day0_results$padj < alpha), ]
day0_sigtab <- cbind(as(day0_sigtab, "data.frame"), as(tax_table(phylo)[rownames(day0_sigtab), ], "matrix"))
write.table(day0_sigtab, "day0_padj05.txt", sep="\t")

##########################################
####Before and after symptom onset
##########################################
##anova_samples.txt is a list of samples included in this test: all samples before and after symptom
##onset in individuals that got sick and corresponding samples pairs from healthy individuals
prune <- read.table ("anova_samples.txt", header=FALSE) 
prune<-as.character(prune[,1])

pruned<-  prune_samples(prune, phylo)

sample_data(pruned)$individual<-factor(sample_data(pruned)$individual)

##time is a factor, either 'before' or 'after' to describe if that samle was taken before or after
##symptom onset
full_model <- phyloseq_to_deseq2(pruned, ~ Final_phenotype + time + Final_phenotype:time)
full_model_deseq<-DESeq(full_model, test="Wald")
full_model_less_reduced_model<-DESeq(full_model_deseq, test="LRT", reduced= ~ Final_phenotype + time)
resultsINT<-results(full_model_less_reduced_model)
summary(resultsINT)
head(resultsINT)

alpha <- 0.05
int_sigtab_05 <- resultsINT[which(resultsINT$padj < alpha), ]
int_sigtab_05 <- cbind(as(int_sigtab_05, "data.frame"), as(tax_table(phylo)[rownames(int_sigtab_05), ], "matrix"))
nrow(int_sigtab_05)
write.table(int_sigtab_05, "ANOVA_padj05.txt", sep="\t")

##########################################
####Testing for effect of phenotype number
##########################################
sample_data(phylo)$individual<-factor(sample_data(phylo)$individual)
sample_data(phylo)$pn<-factor(sample_data(phylo)$pn)

pheno_num <- phyloseq_to_deseq2(phylo, ~ individual + pn)
pheno_num_results <- DESeq(pheno_num, test="Wald")
pheno_num_res <- results(pheno_num_results)

pheno_num_res_0_1<- results(pheno_num_results, contrast=c("pn","0","1"))
pheno_num_res_0_2<- results(pheno_num_results, contrast=c("pn","0","2"))
pheno_num_res_1_2<- results(pheno_num_results, contrast=c("pn","1","2"))
head(pheno_num_res_0_1)
summary(pheno_num_res_0_1)

alpha <- 0.05
results_table_0_1 <- pheno_num_res_0_1[which(pheno_num_res_0_1$padj < alpha), ]
results_table_0_1 <- cbind(as(results_table_0_1, "data.frame"), as(tax_table(phylo)[rownames(results_table_0_1), ], "matrix"))
write.table(results_table_0_1, "H0vsS1.txt", sep="\t")

results_table_0_2 <- pheno_num_res_0_2[which(pheno_num_res_0_2$padj < alpha), ]
results_table_0_2 <- cbind(as(results_table_0_2, "data.frame"), as(tax_table(phylo)[rownames(results_table_0_2), ], "matrix"))
write.table(results_table_0_2, "H0vsS2.txt", sep="\t")

results_table_1_2 <- pheno_num_res_1_2[which(pheno_num_res_1_2$padj < alpha), ]
results_table_1_2 <- cbind(as(results_table_1_2, "data.frame"), as(tax_table(phylo)[rownames(results_table_1_2), ], "matrix"))
write.table(results_table_1_2, "S1vsS2_regrouped.txt", sep="\t")



