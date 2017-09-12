library(phyloseq)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(reshape2)

otutable <- import_biom(BIOMfilename = 'biom_for_phyloseq.biom', 
                        treefilename = 'rep_set_no_chimeras.tre', 
                        parseFunction = parse_taxonomy_greengenes)
#warnings()


mapping <- import_qiime_sample_data(mapfilename = 'map.txt')

phylo <- merge_phyloseq(otutable, mapping)
phylo

phylo_subset = subset_samples(phylo, Phenotype != "Dead")

#First we must rarefy
set.seed(28132)
phyloR = rarefy_even_depth(phylo_subset, sample.size = 20000)

pdf("pheno_ord.pdf")
phylo_ord = ordinate(phyloR, "PCoA", "bray")
p = plot_ordination(phyloR, phylo_ord, color = "Phenotype")
p1= p + geom_point(size = 2, alpha = 0.7) + scale_colour_manual(values = c("orange","purple"))+
theme(axis.title = element_text(color="#666666", size=14)) +
	theme(legend.text=element_text(size=12)) + 
	theme(legend.title = element_text(colour="#666666", size=12))
p1+
  stat_ellipse(type = "t")
dev.off()

