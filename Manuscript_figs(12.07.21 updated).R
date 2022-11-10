### These are the graphs for the scRNA seq figure for the manuscript
library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(cowplot)
setwd("~/Desktop/CTCF mice/scRNA")
#load in the data
load("from Will/toxo.merge.fixed")


###overlay umap 

Idents(toxo.merge) <- "Genotype"
DimPlot(toxo.merge, reduction = "umap", pt.size = 0.00000001, cols = c("#f47b00","#4876FF"),
        shuffle=T, order=c("WT","KO")) + 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
  theme(axis.ticks = element_blank(), axis.line = element_blank()) +
  theme(axis.title = element_blank()) + NoLegend()

#ggsave("overlay.pdf", width = 10, height = 10, units = "in", dpi = "screen")
#system('open \"overlay.pdf"')

#ggsave("Graphs/overlay.png", width = 10, height = 10, units = "in", dpi = "print")
system('open \"Graphs/overlay.png"')

###regular umap
Idents(toxo.merge) <- "seurat_clusters"
cluster_umap <- DimPlot(toxo.merge, reduction = "umap", label = F, label.size = 12) + 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
  theme(axis.ticks = element_blank(), axis.line = element_blank()) +
  theme(axis.title = element_blank()) + NoLegend()
cluster_umap

#ggsave("Graphs/cluster_umap.png", width = 10, height = 10, units = "in", dpi = "print")
#system('open \"Graphs/cluster_umap.png"')

###sample uma
table(toxo.merge$CT.Day.Tissue)
     
Idents(toxo.merge) <- "CT.Day.Tissue"
cluster_umap <- DimPlot(toxo.merge, reduction = "umap", label = F, pt.size=0.001) + 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
  theme(axis.ticks = element_blank(), axis.line = element_blank()) +
  theme(axis.title = element_blank()) + #+ NoLegend()
  theme(legend.text=element_text(size=12))
cluster_umap 


ggsave("Graphs/CT.Day.Tissue_umap.png", width = 8, height = 6, units = "in", dpi = "print")
system('open \"Graphs/CT.Day.Tissue_umap.png"')


##cell cycle umap

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

#Need to change the feature labels to uppercase to match the list of cell cycle genes
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

g2m.genes <- tolower(g2m.genes)
s.genes <- tolower(s.genes)

g2m.genes <- capwords(g2m.genes)
s.genes <- capwords(s.genes)

#Now run the cell cycle function
toxo.merge <- CellCycleScoring(toxo.merge, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

#Now do a feature plot showing the cell cycle idents
cell.cycle <- DimPlot(toxo.merge, group.by = "Phase") +
  theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20)) + 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
  theme(axis.ticks = element_blank(), axis.line = element_blank()) +
  theme(axis.title = element_blank()) + #+ NoLegend()
  theme(legend.text=element_text(size=12)) +
  labs(title = element_blank())
plot_grid(cell.cycle)

ggsave("Graphs/cell_cycle_umap.png", width = 4.5, height = 3.6, units = "in", dpi = "print")
system('open \"Graphs/cell_cycle_umap.png"')


### Proportion Graphs

Idents(toxo.merge) <- "seurat_clusters"
Cell_Proportion.clusters <- prop.table(table(Idents(toxo.merge), toxo.merge$Genotype), margin = 2)
print(Cell_Proportion.clusters)
DF.t <- as.data.frame(Cell_Proportion.clusters)
print(DF.t)
names(DF.t)[1] <- "Cluster"
names(DF.t)[2] <- "Genotype"
names(DF.t)[3] <- "Percentage"
DF.t$Percentage <- DF.t$Percentage *100
DF.t$Percentage <- round(DF.t$Percentage, 1)

NK.clusters <- c(1, 2, 3, 5, 9, 12, 16, 27, 28, 29, 31, 35, 38, 42)
t.D3.clustes <- c(6, 14, 15, 17, 18, 20, 21, 24, 32, 40, 41, 43, 44, 46, 47, 50)
t.d7.clusters <- c(4, 7, 8, 10, 11, 13, 22, 23, 25, 30, 33, 34, 36, 37, 39, 48, 49, 51)

DF.nk <- DF.t[NK.clusters, ]
DF.t.d3 <- DF.t[t.D3.clustes, ]  
DF.t.d7 <- DF.t[t.d7.clusters, ]

proportion_graph.nk <- 
  ggplot(DF.nk, aes(x = Genotype, y = Percentage, fill = Cluster)) +
  geom_bar(color = "black", position = "fill", stat = "identity", width = 0.3, size = 0.25) + 
  theme_classic() +
  geom_text(aes(label = Cluster), position = position_fill(vjust = 0.5), size = 2.5) +
  coord_flip() + 
  NoLegend() +
  theme(axis.title = element_text(size = 9), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7)) + 
  theme(axis.text.x = element_text(vjust = 0.5)) + 
  theme(axis.line = element_line(size = 0.25)) +
  ylab("Proportion") 


proportion_graph.nk
ggsave("Graphs/NK_by_genotype.pdf",width = 2.5, height = 1.2, units = "in", dpi = "print")

#Now let's make a graph for the D3 Th cells
proportion_graph.t.d3 <- 
  ggplot(DF.t.d3, aes(x = Genotype, y = Percentage, fill = Cluster)) +
  geom_bar(color = "black", position = "fill", stat = "identity", width = 0.3, size = 0.25) + 
  theme_classic() +
  geom_text(aes(label = Cluster), position = position_fill(vjust = 0.5), size = 2.5) +
  coord_flip() + 
  NoLegend() +
  theme(axis.title = element_text(size = 9), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7)) + 
  theme(axis.text.x = element_text(vjust = 0.5)) + 
  theme(axis.line = element_line(size = 0.25)) +
  ylab("Proportion") 



proportion_graph.t.d3


ggsave("Graphs/T.D3_by_genotype.pdf",width = 2.5, height = 1.2, units = "in", dpi = "print")

#Now let's make a graph for the D7 Th cells
proportion_graph.t.d7 <- 
  ggplot(DF.t.d7, aes(x = Genotype, y = Percentage, fill = Cluster)) +
  geom_bar(color = "black", position = "fill", stat = "identity", width = 0.3, size = 0.15) + 
  theme_classic() +
  geom_text(aes(label = Cluster), position = position_fill(vjust = 0.5), size = 2.5) +
  coord_flip() + 
  NoLegend() +
  theme(axis.title = element_text(size = 9), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7)) + 
  theme(axis.text.x = element_text(vjust = 0.5)) + 
  theme(axis.line = element_line(size = 0.25)) +
  ylab("Proportion") 


proportion_graph.t.d7

ggsave("Graphs/T.D7_by_genotype.pdf",width = 2.5, height = 1.2, units = "in", dpi = "print")



#######################################################################
#going to redo the violin plots and add in a boxplot to each
#######################################################################
#setwd("/Users/montgomerywf/Desktop/2020_05 CTCF Manuscript Figs/VLN_w_boxplots")

Idents(toxo.merge) <- "Cell Type"
#split the original object into three
obj_list <- SplitObject(toxo.merge, split.by = "Cell Type")

#extract them from the resulting list
th_cells <- obj_list$Th1
nk_cells <- obj_list$NK
as15_cells <- obj_list$as15


Idents(th_cells) <- "Tissue"
th_tis <- Idents(th_cells)
Idents(th_cells) <- "Day"
th_day <- Idents(th_cells)
Idents(th_cells) <- "Genotype"
th_gen <- Idents(th_cells)
th_combo <- paste(th_tis, th_day, th_gen, sep = " ")
th_cells <- AddMetaData(th_cells, metadata = th_combo, col.name = "Manuscript")
th_cells$Manuscript <- factor(x = th_cells$Manuscript, 
                              levels = c("PEC Day7 WT", "PEC Day7 KO", 
                                         "Spleen Day7 WT", "Spleen Day7 KO", 
                                         "PEC Day3 WT", "PEC Day3 KO", 
                                         "Spleen Day3 WT", "Spleen Day3 KO"))


th_plot <- function(gene, save_name) {
  plot <- VlnPlot(th_cells, features = gene, pt.size = 0, group.by = "Manuscript", 
                  cols = c("#4299FB", "#FFA54F", "#4299FB", "#FFA54F", "#4299FB", "#FFA54F", "#4299FB", "#FFA54F")) + 
    theme(axis.title = element_text(size = 9), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7)) +
    ylab("") +
    ggtitle(label = element_blank()) + 
    theme(axis.title.x = element_blank()) + 
    theme(axis.line = element_line(size = 0.2)) + 
    theme(axis.ticks = element_line(size = 0.2)) + 
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "lightgrey",) + 
    NoLegend()
  
  plot$layers[[1]]$aes_params$size = 0.2
  plot$layers[[2]]$aes_params$size = 0.2
  plot
  
  ggsave(file = save_name, width = 2.8, height = 1.8, units = "in", dpi = "retina")
}




Idents(nk_cells) <- "Tissue"
nk_tis <- Idents(nk_cells)
Idents(nk_cells) <- "Day"
nk_day <- Idents(nk_cells)
Idents(nk_cells) <- "Genotype"
nk_gen <- Idents(nk_cells)
nk_combo <- paste(nk_tis, nk_day, nk_gen, sep = " ")
nk_cells <- AddMetaData(nk_cells, metadata = nk_combo, col.name = "Manuscript")

nk_cells$Manuscript <- factor(x = nk_cells$Manuscript, 
                              levels = c("PEC Day3 WT", "PEC Day3 KO", "PEC Day7 WT", "PEC Day7 KO",
                                         "Spleen Day3 WT", "Spleen Day3 KO"))


nk_plot <- function(gene, save_name) {
  plot <- VlnPlot(nk_cells, features = gene, pt.size = 0, group.by = "Manuscript", 
                  cols = c("#4299FB", "#FFA54F","#4299FB", "#FFA54F","#4299FB", "#FFA54F")) + 
    theme(axis.title = element_text(size = 9), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7)) +
    ylab("") +
    theme(legend.text = element_text(size = 7)) + 
    ggtitle(label = element_blank()) + 
    theme(axis.title.x = element_blank()) +
    theme(axis.line = element_line(size = 0.2)) + 
    theme(axis.ticks = element_line(size = 0.2)) + 
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "lightgrey") + 
    NoLegend()
  
  plot$layers[[1]]$aes_params$size = 0.2
  plot$layers[[2]]$aes_params$size = 0.2
  plot
  
  ggsave(file = save_name, width = 2.4, height = 1.8, units = "in", dpi = "retina")
}




Idents(as15_cells) <- "Tissue"
as_tis <- Idents(as15_cells)
Idents(as15_cells) <- "Day"
as_day <- Idents(as15_cells)
Idents(as15_cells) <- "Genotype"
as_gen <- Idents(as15_cells)
as_combo <- paste(as_tis, as_day, as_gen, sep = " ")
as15_cells <- AddMetaData(as15_cells, metadata = as_combo, col.name = "Manuscript")

as15_cells$Manuscript <- factor(x = as15_cells$Manuscript, levels = c("Spleen Day7 WT", "Spleen Day7 KO"))

as_plot <- function(gene, save_name) {
  plot <- VlnPlot(as15_cells, features = gene, pt.size = 0, group.by = "Manuscript", cols = c("#4299FB", "#FFA54F")) + 
    theme(axis.title = element_text(size = 9), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7)) +
    ylab("") +
    ggtitle(label = element_blank()) + 
    theme(axis.title.x = element_blank()) +
    #  scale_x_discrete(labels=c("Spleen Day7" = "Day 7")) +
    theme(axis.line = element_line(size = 0.2)) + 
    theme(axis.ticks = element_line(size = 0.2)) + 
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "lightgrey") + 
    NoLegend()
  
  plot$layers[[1]]$aes_params$size = 0.2
  plot$layers[[2]]$aes_params$size = 0.2
  
  plot
  
  
  ggsave(file = save_name, width = 1, height = 1.8, units = "in", dpi = "retina")
}


#Now that we've defined the functions, let's do them for Ifng and the Ifng response genes
#setwd("~/Desktop/2020_05 CTCF Manuscript Figs/VLN_w_boxplots/")

##This is for the Ifng plots
th_plot("Ifng", "Graphs/th_IFNG_VLN_w_boxplot.pdf")
nk_plot("Ifng", "Graphs/nk_IFNG_VLN_w_boxplot.pdf")
as_plot("Ifng", "Graphs/as15_IFNG_VLN_w_boxplot.pdf")

system('open \"Graphs/th_IFNG_VLN_w_boxplot.pdf"')
system('open \"Graphs/nk_IFNG_VLN_w_boxplot.pdf"')
system('open \"Graphs/as15_IFNG_VLN_w_boxplot.pdf"')

#Now for Igtp
th_plot("Igtp", "Graphs/th_IGTP_w_boxplot.pdf")
nk_plot("Igtp", "Graphs/nk_IGTP_w_boxplot.pdf")
as_plot("Igtp", "Graphs/as15_IGTP_w_boxplot.pdf")

system('open \"Graphs/th_IGTP_w_boxplot.pdf"')
system('open \"Graphs/nk_IGTP_w_boxplot.pdf"')
system('open \"Graphs/as15_IGTP_w_boxplot.pdf"')

#Now for Ifitm3
th_plot("Ifitm3", "Graphs/th_IFITM3_w_boxplot.pdf")
nk_plot("Ifitm3", "Graphs/nk_IFITM3_w_boxplot.pdf")
as_plot("Ifitm3", "Graphs/as15_IFITM3_w_boxplot.pdf")




####Now to create the feature plots for Th17 genes
toxo.merge$Genotype <- factor(x = toxo.merge$Genotype, levels = c("WT", "KO"))



FeaturePlot(toxo.merge, features = "Rorc", split.by = "Genotype", sort.cell = TRUE) 
ggsave(filename = "Rorc.pdf", width = 9, height = 4.5, dpi = "print")

FeaturePlot(toxo.merge, features = "Ccr6", split.by = "Genotype", sort.cell = TRUE) 
ggsave(filename = "Ccr6.pdf", width = 9, height = 4.5, dpi = "print")

FeaturePlot(toxo.merge, features = "Il17f", split.by = "Genotype", sort.cell = TRUE)
ggsave(filename = "Il17f.pdf", width = 9, height = 4.5, dpi = "print")

FeaturePlot(toxo.merge, features = "Tmem176a", split.by = "Genotype", sort.cell = TRUE)
ggsave(filename = "Tmem176a.pdf", width = 9, height = 4.5, dpi = "print")

FeaturePlot(toxo.merge, features = "Tmem176b", split.by = "Genotype", sort.cell = TRUE)
ggsave(filename = "Tmem176b.pdf", width = 9, height = 4.5, dpi = "print")

###Now to create the feature plots for cluster 12

genes <- list("Mt1", "Mt3", "Ccne1", "Mt2", "Ybx3")

FeaturePlot(toxo.merge, features = "Ccne1", split.by = "Genotype")
ggsave(filename = "Ccne1.pdf", width = 9, height = 4.5, dpi = "print")

FeaturePlot(toxo.merge, features = "Ccne2", split.by = "Genotype")
ggsave(filename = "Ccne2.pdf", width = 9, height = 4.5, dpi = "print")

FeaturePlot(toxo.merge, features = c("Mt1", "Mt3"), split.by = "Genotype")
ggsave(filename = "MT1_3.pdf", width = 9, height = 9, dpi = "print")

FeaturePlot(toxo.merge, features = "Mt2", split.by = "Genotype", sort.cell = TRUE)
ggsave(filename = "MT2.pdf", width = 9, height = 4.5, dpi = "print")


##########################################################################
# Build a cluster tree for toxo.merge
##########################################################################

cluster_tree_pca <- BuildClusterTree(toxo.merge, dims = 1:75)
PlotClusterTree(cluster_tree_pca)


###################################################################################################
#Going to make a new version of this where the clusters are numbered by the hierarchical tree 
###################################################################################################
Idents(toxo.merge) <- "seurat_clusters"
cluster_tree_pca_nn <- BuildClusterTree(toxo.merge, dims = 1:75)
cluster_tree_pca_nn <- RenameIdents(cluster_tree_pca_nn, '10' = '1', '24' = '2', '6' = '3', '9' = '4',
                                 '3' = '5', '12' = '6', '22' = '7', '7' = '8', '21' = '9', '5' = '10', 
                                 '14' = '11', '17' = '12', '18' = '13', '20' = '14', '16' = '15', 
                                 '19' = '16', '0' = '17', '2' = '18', '1' = '19', '4' = '20', '13' = '21',
                                 '8' = '22', '11' = '23', '15' = '24', '23' = '25', '25' = '26')
cluster_tree_pca_nn <- BuildClusterTree(cluster_tree_pca_nn, dims = 1:75)
PlotClusterTree(cluster_tree_pca_nn)


#####Now to make a new umap to match the above
Idents(toxo.merge) <- "seurat_clusters"
cluster_umap <- DimPlot(cluster_tree_pca_nn, reduction = "umap", label = TRUE, label.size = 12) + 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
  theme(axis.ticks = element_blank(), axis.line = element_blank()) +
  theme(axis.title = element_blank()) + NoLegend()

cluster_umap
#ggsave("Graphs/cluster_umap_renamed.png", width = 10, height = 10, units = "in", dpi = "print")
#system('open \"Graphs/cluster_umap_renamed.png"')
#ggsave("cluster_umap_renamed.pdf", width = 15, height = 15, units = "in", dpi = "print")


#################################################################################################################
# Now to make a bar graph to show comparison of the clusters like they did in the paper Nilisha presented
#################################################################################################################
Cell_Proportion.Genotype <- prop.table(table(Idents(cluster_tree_pca_nn), cluster_tree_pca_nn$Genotype), margin = 2)
print(Cell_Proportion.Genotype)
DF.g <- as.data.frame(Cell_Proportion.Genotype)
print(DF.g)
names(DF.g)[1] <- "Cluster"
names(DF.g)[2] <- "Genotype"
names(DF.g)[3] <- "Percentage"
DF.g$Percentage <- DF.g$Percentage * 100

DF.g$Genotype <- factor(x = DF.g$Genotype, levels = c("KO", "WT"))
DF.g$Cluster <- factor(x = DF.g$Cluster, levels = rev(as.character(1:27)))

proportion_graph.genotype <-
  ggplot() +
  geom_bar(data = DF.g, aes(y = Percentage, x = Cluster, fill = Genotype), stat = "identity", position = 
             position_dodge(), width = 0.8) + theme_classic() + scale_y_continuous(expand = c(0,0), limits = c(0, 12)) + 
  scale_fill_manual("Legend", values = c("WT" = "#4662E5", "KO" = "#FFA54F")) + 
  theme(axis.title.x = element_text(size = 10), axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 8)) + 
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5))  +
  theme(axis.line.x = element_line(size = 0.2), axis.ticks.x = element_line(size = 0.2)) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank()) +
  theme(legend.key.height = unit(0.25, "line")) + 
  theme(legend.key.width = unit(0.25, "line")) +
  theme(legend.text = element_text(size = 8)) + 
  theme(legend.title = element_blank()) + 
  coord_flip() + NoLegend()

  

proportion_graph.genotype

ggsave("Graphs/proportion_graph.pdf", width = 1.6, height = 3, units = "in", dpi = "print")
system('open \"Graphs/proportion_graph.pdf"')

#################################################################################################################
# Now to make a new log ratio graph for the reorganized data
#################################################################################################################


prop_table <- prop.table(table(Idents(cluster_tree_pca_nn), cluster_tree_pca_nn$Genotype), margin = 2)
prop_KO <- as.numeric(prop_table[, 1])
prop_WT <- as.numeric(prop_table[, 2])
prop_KO <- (prop_KO * 100)+1
prop_WT <- (prop_WT * 100)+1
prop_log2 <- log2(prop_KO/prop_WT)
prop_DF <- data.frame("Cluster" = as.character(1:26), "Log2 Proportion" = prop_log2)
positions <- rev(as.character(1:26))
range(prop_DF[,2])

prop_graph <- ggplot() +
  geom_bar(aes(y = Log2.Proportion, x = Cluster), data = prop_DF, stat = "identity", position = position_dodge(width = 0.7)) +
  coord_flip() +
  theme_classic() + 
  theme(axis.line.x = element_line(size = 0.2), axis.ticks.x = element_line(size = 0.2)) +
  theme(axis.line.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_text(size = 8)) + 
  scale_x_discrete(limits = positions) +
  scale_y_continuous(limits = c(-1, 2)) + 
  geom_abline(slope=0, intercept=0.58,  col = "grey",lty=2) +
  geom_abline(slope=0, intercept=-0.58,  col = "grey",lty=2) +
  theme(axis.title = element_text(size = 10), axis.text.x = element_text(size = 8)) + 
  labs(y = "Log2 of KO/WT")

prop_graph

# identify interest clusters (at least 1% and 1.5 FC for the proportion)
abs(prop_log2)> 0.58

ggsave(filename = "Graphs//log2_proportion_graph_renamed_new.pdf", width = 1.6, height = 3, units = "in", dpi = "print")
system('open \"Graphs/log2_proportion_graph_renamed_new.pdf"')


#################################################################################################################
#Plot gene expression
#################################################################################################################

Th1_genes = c("Ifng","Ifitm2","Ifitm3", "Tbx21", "Ecm1", "Il12rb2","Slamf1","Cxcr6", "Gstt1")
featureplot_mygenes_Th1=FeaturePlot(cluster_tree_pca_nn, features = Th1_genes, by.col=	3,pt.size=0.05,
                                          # split.by = "orig.ident"
                                           )
ggsave("Graphs/featureplot_Th1_genes.png", plot = featureplot_mygenes_Th1, width = 25, height = 20, units = "cm")
system('open \"Graphs/featureplot_Th1_genes.png"')


Th17_genes = c("Il17a","Il17f","Rorc", "Tmem176a", "Tmem176b", "Ccr6", "Acsbg1", "Actn2","Rora")
featureplot_mygenes_Th17=FeaturePlot(cluster_tree_pca_nn, features = Th17_genes, by.col=	3,pt.size=0.05,
                                           # split.by = "orig.ident"
)
ggsave("Graphs/featureplot_Th17_genes.png", plot = featureplot_mygenes_Th17, width = 25, height = 20, units = "cm")
system('open \"Graphs/featureplot_Th17_genes.png"')


Th17_genes_2 = c("Ifit1", "Ifit3", "Actn2", "Slfn5", "Oasl2", "S100a4", "Igtp", "Isg15", "Stat1") 

featureplot_mygenes_Th17=FeaturePlot(cluster_tree_pca_nn, features = Th17_genes_2, by.col=	3,pt.size=0.05,
                                     # split.by = "orig.ident"
)

ggsave("Graphs/featureplot_Th17_genes_2.png", plot = featureplot_mygenes_Th17, width = 25, height = 20, units = "cm")
system('open \"Graphs/featureplot_Th17_genes_2.png"')



featureplot_mygenes_KOvsWT=FeaturePlot(cluster_tree_pca_nn, features = c("Ifitm3"), #by.col=	3,
                                           split.by = "Genotype")
featureplot_mygenes_KOvsWT

ridgeplot_lineage_genes=RidgePlot(cluster_tree_pca_nn, features = Th17_genes)
ridgeplot_lineage_genes
#ggsave("ridgeplot_lineage_genes.png", plot = ridgeplot_lineage_genes, width = 20, height = 40, units = "cm")

mygenes = c("Il17a","Il17f","Rorc", "Ifng", "Ifitm3", "Cxcr6")
dotplot_mygenes=DotPlot(cluster_tree_pca_nn, features = mygenes) + #+ coord_flip()
  theme(axis.ticks.y = element_blank(), axis.title.y = element_blank()) +
  scale_y_discrete(limits=rev) +
  theme(legend.key.height = unit(0.25, "line")) + 
  theme(legend.key.width = unit(0.25, "line")) +
  theme(legend.text = element_text(size = 8)) + 
  theme(legend.title = element_blank())
dotplot_mygenes


Th_d7_subsets <- subset(cluster_tree_pca_nn, idents = c("1", "2", "3", "4", "5" ,"6", "7", "8", "9"))
mygenes = c("Ifng", "Tbx21", "Ifitm2", "Ifitm3", "Cxcr6", "Il12rb2")
dotplot_mygenes=DotPlot(Th_d7_subsets, features = mygenes, #split.by= "Genotype", 
                       ) +  
  coord_flip() +
  theme(axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(),
        axis.line = element_line(size = 0.2),axis.ticks = element_line(size = 0.2), 
        axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 9)) +
  scale_x_discrete(limits=rev) +
  theme(legend.key.height = unit(0.25, "line"), 
        legend.key.width = unit(0.25, "line"),
        legend.text = element_text(size = 8), 
        legend.title = element_blank())
dotplot_mygenes

ggsave("Graphs/dotplot_Th1_genes.pdf", plot = dotplot_mygenes, width = 12, height = 5, units = "cm")
system('open \"Graphs/dotplot_Th1_genes.pdf"')

Th_d3_subsets <- subset(cluster_tree_pca_nn, idents = c("10", "11", "12", "13", "14", "15" ,"16", "21"))
mygenes_1 = c("Il17a","Il17f","Rorc")

dotplot_mygenes=DotPlot(Th_d3_subsets, features = mygenes_1, #split.by= "Genotype", 
) +  
  coord_flip() +
  theme(axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(),
        axis.line = element_line(size = 0.2),axis.ticks = element_line(size = 0.2), 
        axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 9)) +
  scale_x_discrete(limits=rev) +
  theme(legend.key.height = unit(0.25, "line"), 
        legend.key.width = unit(0.25, "line"),
        legend.text = element_text(size = 7), 
        legend.title = element_blank())
dotplot_mygenes

ggsave("Graphs/dotplot_Th17_genes_1.pdf", plot = dotplot_mygenes, width = 9, height = 3, units = "cm")
system('open \"Graphs/dotplot_Th17_genes_1.pdf"')

mygenes_2 = c("Rora", "Tmem176a","Tmem176b")

dotplot_mygenes=DotPlot(Th_d3_subsets, features = mygenes_2, #split.by= "Genotype", 
) +  
  coord_flip() +
  theme(axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(),
        axis.line = element_line(size = 0.2),axis.ticks = element_line(size = 0.2), 
        axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 9)) +
  scale_x_discrete(limits=rev) +
  theme(legend.key.height = unit(0.25, "line"), 
        legend.key.width = unit(0.25, "line"),
        legend.text = element_text(size = 7), 
        legend.title = element_blank())
dotplot_mygenes

ggsave("Graphs/dotplot_Th17_genes_2.pdf", plot = dotplot_mygenes, width = 10, height = 3, units = "cm")
system('open \"Graphs/dotplot_Th17_genes_2.pdf"')

#################################################################################################################
#Create stacked violin plot (https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/)
#################################################################################################################

#VlnPlot(all_gex_s_Filt_Scaled, features = top10)
library(patchwork)

## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "lightgrey")
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

Th_d7_subsets <- subset(cluster_tree_pca_nn, idents = c("1", "2", "3", "4", "5" ,"6", "7", "8", "9"))
mygenes = c("Ifng", "Ifitm3", "Cxcr6")
StackedVlnPlot(obj = Th_d7_subsets, features = mygenes)

Th_d3_subsets <- subset(cluster_tree_pca_nn, idents = c("10", "11", "12", "13", "14", "15" ,"16", "21"))
mygenes = c("Il17a","Il17f","Rorc")
StackedVlnPlot(obj = Th_d3_subsets, features = mygenes)



#ggsave("vlnplot_mygenes_all.png", plot = vlnplot_mygenes_all, width = 50, height = 20, units = "cm")

dotplot_mygenes_KOvsWT_all=DotPlot(cluster_tree_pca_nn, features = mygenes, split.by = "orig.ident")

#ggsave("dotplot_mygenes_KOvsWT_all.png", plot = dotplot_mygenes_KOvsWT_all, width = 30, height = 40, units = "cm")

doheatmap_mygenes_all=DoHeatmap(subset(Th_d3_subsets), features = mygenes, size = 3)

#ggsave("doheatmap_mygenes_all.png", plot = doheatmap_mygenes_all, width = 50, height = 30, units = "cm")


#################################################################################################################
#Create markers for each cluster
#################################################################################################################

# all.markers <- FindAllMarkers(cluster_tree_pca_nn, only.pos = FALSE)
# write.csv(all.markers,"all.markers.csv", row.names = TRUE)
all.markers<-read.csv(file = "all.markers.csv")

markers_table_top5 <- all.markers %>% group_by(cluster) %>% top_n(n=5, wt = avg_logFC)
write.csv(markers_table_top10,"all.markers_top5.csv", row.names = TRUE)
markers_table_top10 <- all.markers %>% group_by(cluster) %>% top_n(n=10, wt = avg_logFC)
write.csv(markers_table_top10,"all.markers_top10.csv", row.names = TRUE)
markers_table_top20 <- all.markers %>% group_by(cluster) %>% top_n(n=20, wt = avg_logFC)
write.csv(markers_table_top20,"all.markers_top20.csv", row.names = TRUE)

#################################################################################################################
# Make Heatmap table
#################################################################################################################

cluster_tree_pca_nn
png("Graphs/Cluster_heatmap_top5.png", height=20, width=16, units = "in", res = 300)
DoHeatmap(cluster_tree_pca_nn, features = paste(markers_table_top5$gene), label = FALSE) + NoLegend() + 
  theme(axis.text.y = element_text(size = 14))
dev.off()
system('open \"Graphs/Cluster_heatmap_top5.png"')

png("Graphs/Cluster_heatmap_top10.png", height=24, width=16, units = "in", res = 300)
DoHeatmap(cluster_tree_pca_nn, features = paste(markers_table_top10$gene), label = FALSE) + NoLegend() + 
  theme(axis.text.y = element_text(size = 12))
dev.off()
system('open \"Graphs/Cluster_heatmap_top10.png"')

png("Graphs/Cluster_heatmap_top20.png", height=40, width=16, units = "in", res = 300)
DoHeatmap(cluster_tree_pca_nn, features = paste(markers_table_top20$gene), label = FALSE) + NoLegend() + 
  theme(axis.text.y = element_text(size = 12))
dev.off()
system('open \"Graphs/Cluster_heatmap_top20.png"')

#################################################################################################################
#Create markers for T cell d7 clusters
#################################################################################################################

Th_d7_subsets <- subset(cluster_tree_pca_nn, idents = c("1", "2", "3", "4", "5" ,"6", "7", "8", "9"))
# Th.d7.markers <- FindAllMarkers(Th_d7_subsets, only.pos = FALSE)
# write.csv(Th.d7.markers,"Th.d7.markers.csv", row.names = TRUE)
Th.d7.markers<-read.csv(file = "Th.d7.markers.csv")

Th.d7.markers_top10 <- Th.d7.markers %>% group_by(cluster) %>% top_n(n=10, wt = avg_logFC)
write.csv(Th.d7.markers_top10,"Th.d7.markers_top10.csv", row.names = TRUE)
Th.d7.markers_top20 <- Th.d7.markers %>% group_by(cluster) %>% top_n(n=20, wt = avg_logFC)
write.csv(Th.d7.markers_top20,"Th.d7.markers_top20.csv", row.names = TRUE)



#################################################################################################################
# Make Heatmap table for Th d7
#################################################################################################################
Th_d7_subsets
png("Graphs/Th.d7.markers_top10_heatmap.png", height=12, width=16, units = "in", res = 300)
DoHeatmap(Th_d7_subsets, features = paste(Th.d7.markers_top10$gene), label = FALSE) + NoLegend() + 
  theme(axis.text.y = element_text(size = 12))
dev.off()
system('open \"Graphs/Th.d7.markers_top10_heatmap.png"')

png("Graphs/Th.d7.markers_top20_heatmap.png", height=24, width=16, units = "in", res = 300)
DoHeatmap(Th_d7_subsets, features = paste(Th.d7.markers_top20$gene), label = FALSE) + NoLegend() + 
  theme(axis.text.y = element_text(size = 12))
dev.off()
system('open \"Graphs/Th.d7.markers_top20_heatmap.png"')

#################################################################################################################
#Create markers for T cell d3 clusters
#################################################################################################################

Th_d3_subsets <- subset(cluster_tree_pca_nn, idents = c("10", "11", "12", "13", "14", "15" ,"16", "21"))
Th.d3.markers <- FindAllMarkers(Th_d3_subsets, only.pos = FALSE)
# write.csv(Th.d3.markers,"Th.d3.markers.csv", row.names = TRUE)
Th.d3.markers<-read.csv(file = "Th.d3.markers.csv")

Th.d3.markers_top10 <- Th.d3.markers %>% group_by(cluster) %>% top_n(n=10, wt = avg_logFC)
write.csv(Th.d3.markers_top10,"Th.d3.markers_top10.csv", row.names = TRUE)
Th.d3.markers_top20 <- Th.d3.markers %>% group_by(cluster) %>% top_n(n=20, wt = avg_logFC)
write.csv(Th.d3.markers_top20,"Th.d3.markers_top20.csv", row.names = TRUE)



#################################################################################################################
# Make Heatmap table for Th d3
#################################################################################################################
Th_d3_subsets
png("Graphs/Th.d3.markers_top10_heatmap.png", height=12, width=16, units = "in", res = 300)
DoHeatmap(Th_d3_subsets, features = paste(Th.d3.markers_top10$gene), label = FALSE) + NoLegend() + 
  theme(axis.text.y = element_text(size = 12))
dev.off()
system('open \"Graphs/Th.d3.markers_top10_heatmap.png"')

png("Graphs/Th.d3.markers_top20_heatmap.png", height=24, width=16, units = "in", res = 300)
DoHeatmap(Th_d3_subsets, features = paste(Th.d3.markers_top20$gene), label = FALSE) + NoLegend() + 
  theme(axis.text.y = element_text(size = 12))
dev.off()
system('open \"Graphs/Th.d3.markers_top20_heatmap.png"')


#################################################################################################################
# Compare WT and KO gene expression in single cluster
#################################################################################################################

markersingroup14 <- FindMarkers(cluster_tree_pca_nn, ident.1 = "WT", ident.2 = "KO", group.by = 'Genotype', subset.ident = 14)
markersingroup14[order(markersingroup14$avg_logFC),]
markersingroup6 <- FindMarkers(cluster_tree_pca_nn, ident.1 = "WT", ident.2 = "KO", group.by = 'Genotype', subset.ident = 6)
markersingroup6[order(markersingroup6$avg_logFC),]
