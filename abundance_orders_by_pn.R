library(phyloseq)
library(DESeq2)
library(ggplot2)
theme_set(theme_bw())
library(dplyr)
library(tidyr)
library(RColorBrewer)

otutable <- import_biom(BIOMfilename = 'biom_for_phyloseq.biom', 
                        treefilename = 'rep_set_no_chimeras.tre', 
                        parseFunction = parse_taxonomy_greengenes)
#warnings()

mapping <- import_qiime_sample_data(mapfilename = 'map.txt')

phylo <- merge_phyloseq(otutable, mapping)
phylo


###############################################
#First we must rarefy
set.seed(28132)
phyloR = rarefy_even_depth(phylo, sample.size = 20000)

#Merge on basis of pheno number so we can make relative abundance graph
phyloRm = merge_samples(phyloR, "Pheno_num")

#get rid of very low abundance taxa
tax_table(phyloRm)[is.na(tax_table(phyloRm))] <- 0

glom <- tax_glom(phyloRm, "Order")
glommed <- filter_taxa(glom, function(x) sum(x > 800) > (0.33*length(x)), TRUE)

#turn into proportional counts
phyloRt <- transform_sample_counts(glommed, function(x) 1 * x/sum(x))

sample_data(phyloRt)$Pheno_num <- as.factor(sample_data(phyloRt)$Pheno_num)

#pick colors
orderPalette2<-c("#081D58","#57BEC0","#2167AC","#1D8CBD","#DBF1B2","#23479D","#1D2D83",
"#B9E3B5","#FFFFD9","#85CFBA","#F1F9BB", "#1D2D83","#33A8C2")

labels=c("Unknown","Alteromonadales","Campylobacterales", "Cytophagales",
"Flavobacteriales","Oceanospirillales", "Pseudomonadales","Rhizobiales","Rhodobacterales","Spirochaetales","UA01","Vibrionales")

pdf("order_abundance_sym_num.pdf")
p = plot_bar(phyloRt, "Pheno_num", fill="Order")
p + geom_bar(aes(color=Order, fill=Order), stat="identity", position="stack")+xlab("Symptom progression number")+
ylab("Proportion")  +
theme(axis.title = element_text(color="#666666", face="bold", size=16)) +
theme(axis.text.x = element_text(size=12, angle=360)) +
theme(legend.text=element_text(size=12)) + 
scale_color_manual(values=orderPalette2, labels=labels) + 
scale_fill_manual(values= orderPalette2,labels=labels)  +
theme(legend.title=element_text(colour="#666666", size=12, face="bold"))
dev.off()
