# Loading Packages
library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggpubr)
library(broom)
library(RColorBrewer)
library(RgoogleMaps)
library(metagMisc)
library(shiny)
library(viridis)
library(svglite)
library(gdtools)
library(VennDiagram)
library(indicspecies)
library(geosphere)
library(ecodist)

# Color palette
pal = c("#6b5456","#6abf2a","#70acbe","#01c95b","#c00014","#31332f","#f7d000","#abba00")
col_island <- viridis(7)
names(col_island) <- c("Hantu", "Jong", "Kusu", "Raffles Lighthouse", "Semakau", "Sisters", "Sultan Shoal")

theme_set(theme_bw())

#### Analysis ####
ps1 <- readRDS("./Output/Pl/final_phyloseq_object.RDS")
taxtable.update <- as.data.frame(ps1@tax_table)
for(j in 1:ncol(taxtable.update)){
  taxtable.update[,j] <- as.character(taxtable.update[,j])
  for(i in 1:nrow(taxtable.update)){
    if(!is.na(taxtable.update[i, j])){
      textLabel <- strsplit(as.character(taxtable.update[i, j]), split = '_')[[1]][3:length(strsplit(as.character(taxtable.update[i, j]), split = '_')[[1]])]
      taxtable.update[i, j] <- paste(textLabel, collapse='_')
    }
  }
}


# Data Prep for FUNGuild
taxafinal <- taxtable.update
taxafinal$taxonomy <- taxafinal$Species
for(i in 1:nrow(taxafinal)){
  taxafinal$taxonomy[i] <- paste(taxafinal[i,-ncol(taxafinal)], collapse = ' ')
}
OTUTable <- as.data.frame(ps1@otu_table)
OTUTable <- as.data.frame(t(OTUTable))
OTUTable$OTU_ID <- row.names(OTUTable)
OTUTable$taxonomy <- taxafinal$taxonomy
OTUTable <- OTUTable[,98:99]
write.csv(OTUTable, "./Output/Pl/FUNGuild_OTU_Table.csv")


# Normalize (relative abundance) ####
ps1ra <- transform_sample_counts(ps1, function(otu){otu/sum(otu)})



### Rarefaction Curves ###
grp <- factor(ps1@sam_data$Island)
cols <- col_island[grp]
rarefaction_curve <- rarecurve(ps1@otu_table, step = 20, col = cols, label = FALSE)
Nmax <- sapply(rarefaction_curve, function(x) max(attr(x, "Subsample")))
Smax <- sapply(rarefaction_curve, max)
svg(filename = "./Output/Pl/Analysis/Rarefaction Curves.svg", )
plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = "Number of Sequences",
     ylab = "Number of ASVs", type = "n",
     main = "Rarefaction Curves")
for (i in seq_along(rarefaction_curve)) {
  N <- attr(rarefaction_curve[[i]], "Subsample")
  lines(N, rarefaction_curve[[i]], col = cols[i])
}
legend(60000, 20, legend = names(col_island), col = col_island, lty = 1, cex = 0.8, box.lty = 1)
dev.off()





### Relative Abundance Bar Plots by Location and Structure ###

# Creating a dataframe of the otu table which also includes the Island variables
combineddf <- as.data.frame(ps1@otu_table)
combineddf$Island <- ps1@sam_data$Island
##############
# Creating a phyloseq object where samples are grouped and merged according to Island
island_otu <- combineddf %>%
  group_by(Island) %>%
  summarise_all(.funs=sum)
island_otu <- as.data.frame(island_otu)
row.names(island_otu) <- c("Samples_Hantu", "Samples_Jong", "Samples_Kusu", "Samples_Raffles", "Samples_Semakau", "Samples_Sisters", "Samples_Sultan")
island_otu <- island_otu[,-1]
island_meta <- data.frame(Island = c("Hantu", "Jong", "Kusu", "Raffles Lighthouse", "Semakau", "Sisters", "Sultan Shoal"))
row.names(island_meta) <- NULL
row.names(island_meta) <- c("Samples_Hantu", "Samples_Jong", "Samples_Kusu", "Samples_Raffles", "Samples_Semakau", "Samples_Sisters", "Samples_Sultan")
speciesList <- taxtable.update$Species
genusList <- taxtable.update$Genus
genusSpeciesList <- speciesList
for(i in 1:length(speciesList)){
  if(!is.na(speciesList[i])){
    genusSpeciesList[i] <- paste(genusList[i], speciesList[i], sep = ' ')
  }
}
taxtable.update$Genus_species <- genusSpeciesList

island_ps <- phyloseq(otu_table(island_otu, taxa_are_rows=FALSE), 
                      sample_data(island_meta), 
                      tax_table(ps1@tax_table))

# Calculating Relative abundances and plotting bar plots according to Location and Structure
ps_Island_ra <- transform_sample_counts(island_ps, function(otu){otu/sum(otu)})


# Defining function to make bar charts without black lines separating samples. Based on phyloseq function "plot_bar".
simple_plot_bar = function (physeq, x = "Sample", y = "Abundance", fill = NULL, title = NULL, 
                            facet_wrap = NULL) {
  mdf = psmelt(physeq)
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  p = p + geom_bar(stat = "identity", position = "stack")
  p = p + theme_bw() + theme(axis.text=element_text(size=15), axis.title=element_text(size=17,face="bold"), 
                             axis.text.x = element_text(angle = 70, hjust = 1))
  p = p + labs(y = "Relative Abundance")
  p = p + guides(guide_legend(ncol = 1), fill = guide_legend(ncol = 3))
  if (!is.null(facet_wrap)) {
    p <- p + facet_wrap(facet_wrap, nrow = 1) + theme(strip.text = element_text(size=15))
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

# Making stacked bar charts for relative abundance of taxa

# According to Phylum
ra_Phylum_barplot_island <- simple_plot_bar(ps_Island_ra, x="Island", fill="Phylum") + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = sort(unique(taxtable.update$Phylum)))
ggsave(ra_Phylum_barplot_island, filename = "./Output/Pl/Analysis/Relative Abundance of Phylum by Island - Bar Plot.svg", dpi=300, width = 12, height = 10)

ra_Phylum_barplot_all <- simple_plot_bar(ps1ra, x="SampleID", fill="Phylum") + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = sort(unique(taxtable.update$Phylum)))
ggsave(ra_Phylum_barplot_all, filename = "./Output/Pl/Analysis/Relative Abundance of Phylum by Sample - Bar Plot.svg", dpi=300, width = 25, height = 10)

ra_Phylum_barplot_combined <- simple_plot_bar(ps1ra, x="sample_Species", fill="Phylum") + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = sort(unique(taxtable.update$Phylum)))
ggsave(ra_Phylum_barplot_combined, filename = "./Output/Pl/Analysis/Relative Abundance of Phylum Combined - Bar Plot.svg", dpi=300, width = 12, height = 10)

# According to Class
ra_Class_barplot_island <- simple_plot_bar(ps_Island_ra, x="Island", fill="Class") + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = sort(unique(taxtable.update$Class)))
ggsave(ra_Class_barplot_island, filename = "./Output/Pl/Analysis/Relative Abundance of Class by Island - Bar Plot.svg", dpi=300, width = 12, height = 10)

ra_Class_barplot_all <- simple_plot_bar(ps1ra, x="SampleID", fill="Class") + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = sort(unique(taxtable.update$Class)))
ggsave(ra_Class_barplot_all, filename = "./Output/Pl/Analysis/Relative Abundance of Class by Sample - Bar Plot.svg", dpi=300, width = 25, height = 10)

ra_Class_barplot_combined <- simple_plot_bar(ps1ra, x="sample_Species", fill="Class") + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = sort(unique(taxtable.update$Class)))
ggsave(ra_Class_barplot_combined, filename = "./Output/Pl/Analysis/Relative Abundance of Class Combined - Bar Plot.svg", dpi=300, width = 12, height = 10)

# According to Order
ra_Order_barplot_island <- simple_plot_bar(ps_Island_ra, x="Island", fill="Order") + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = sort(unique(taxtable.update$Order)))
ggsave(ra_Order_barplot_island, filename = "./Output/Pl/Analysis/Relative Abundance of Order by Island - Bar Plot.svg", dpi=300, width = 12, height = 10)

ra_Order_barplot_all <- simple_plot_bar(ps1ra, x="SampleID", fill="Order") + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = sort(unique(taxtable.update$Order)))
ggsave(ra_Order_barplot_all, filename = "./Output/Pl/Analysis/Relative Abundance of Order by Sample - Bar Plot.svg", dpi=300, width = 25, height = 10)

ra_Order_barplot_combined <- simple_plot_bar(ps1ra, x="sample_Species", fill="Order") + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = sort(unique(taxtable.update$Order)))
ggsave(ra_Order_barplot_combined, filename = "./Output/Pl/Analysis/Relative Abundance of Order Combined - Bar Plot.svg", dpi=300, width = 12, height = 10)

# According to Family
ra_Family_barplot_island <- simple_plot_bar(ps_Island_ra, x="Island", fill="Family") + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = sort(unique(taxtable.update$Family)))
ggsave(ra_Family_barplot_island, filename = "./Output/Pl/Analysis/Relative Abundance of Family by Island - Bar Plot.svg", dpi=300, width = 12, height = 10)

ra_Family_barplot_all <- simple_plot_bar(ps1ra, x="SampleID", fill="Family") + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = sort(unique(taxtable.update$Family)))
ggsave(ra_Family_barplot_all, filename = "./Output/Pl/Analysis/Relative Abundance of Family by Sample - Bar Plot.svg", dpi=300, width = 25, height = 10)

ra_Family_barplot_combined <- simple_plot_bar(ps1ra, x="sample_Species", fill="Family") + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = sort(unique(taxtable.update$Family)))
ggsave(ra_Family_barplot_combined, filename = "./Output/Pl/Analysis/Relative Abundance of Family Combined - Bar Plot.svg", dpi=300, width = 12, height = 10)

# According to Genus
ra_Genus_barplot_island <- simple_plot_bar(ps_Island_ra, x="Island", fill="Genus") + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = sort(unique(taxtable.update$Genus)))
ggsave(ra_Genus_barplot_island, filename = "./Output/Pl/Analysis/Relative Abundance of Genus by Island - Bar Plot.svg", dpi=300, width = 12, height = 10)

ra_Genus_barplot_all <- simple_plot_bar(ps1ra, x="SampleID", fill="Genus") + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = sort(unique(taxtable.update$Genus)))
ggsave(ra_Genus_barplot_all, filename = "./Output/Pl/Analysis/Relative Abundance of Genus by Sample - Bar Plot.svg", dpi=300, width = 25, height = 10)

ra_Genus_barplot_combined <- simple_plot_bar(ps1ra, x="sample_Species", fill="Genus") + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = sort(unique(taxtable.update$Genus)))
ggsave(ra_Genus_barplot_combined, filename = "./Output/Pl/Analysis/Relative Abundance of Genus Combined - Bar Plot.svg", dpi=300, width = 12, height = 10)

# According to Species
sorter <- c(sort(unique(taxtable.update$Species), index.return = TRUE)$ix + 1, 1)

ra_Species_barplot_island <- simple_plot_bar(ps_Island_ra, x="Island", fill="Species") + theme(legend.text=element_text(size=rel(1.3)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = unique(taxtable.update$Genus_species)[sorter])
ggsave(ra_Species_barplot_island, filename = "./Output/Pl/Analysis/Relative Abundance of Species by Island - Bar Plot.svg", dpi=300, width = 12, height = 10)

ra_Species_barplot_all <- simple_plot_bar(ps1ra, x="SampleID", fill="Species") + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = unique(taxtable.update$Genus_species)[sorter])
ggsave(ra_Species_barplot_all, filename = "./Output/Pl/Analysis/Relative Abundance of Species by Sample - Bar Plot.svg", dpi=300, width = 25, height = 10)

ra_Species_barplot_combined <- simple_plot_bar(ps1ra, x="sample_Species", fill="Species") + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = unique(taxtable.update$Genus_species)[sorter])
ggsave(ra_Species_barplot_combined, filename = "./Output/Pl/Analysis/Relative Abundance of Species Combined - Bar Plot.svg", dpi=300, width = 12, height = 10)




### Heatmap ###
melted_ps <- psmelt(ps1ra)


## By Phylum
factor_melted_ps_by_sample <- unique(melted_ps$SampleID) #factoring by sample
dataframe <- data.frame()  #creating dataframe to fill in information

for (i in factor_melted_ps_by_sample) {                                    #For each sample,
  sub_ps <- melted_ps[melted_ps$SampleID == i,]  #subset dataframe.
  
  per_phylum_abundance <- aggregate(Abundance ~ Phylum, data = sub_ps, sum)  #Aggregate abundance based on each phylum
  
  per_phylum_abundance$SampleID <- sub_ps$SampleID[1:nrow(per_phylum_abundance)] #Adding sample name to each phylum row
  per_phylum_abundance$Island <- sub_ps$Island[1:nrow(per_phylum_abundance)]
  
  dataframe <- rbind(dataframe, per_phylum_abundance)                       #store this in a dataframe for each row
}

#Sorting dataframe in order of island
ordered_df <- dataframe[order(dataframe$Island),]

#Plotting
table(ordered_df$Island) / sum(unique(ordered_df$Phylum) == unique(ordered_df$Phylum))
heatplot <- ggplot(ordered_df, aes(reorder(SampleID, -desc(Island)), reorder(Phylum, desc(Phylum)))) +
  geom_tile(aes(fill = Abundance)) + 
  labs(y = "Phylum", x = "Samples") + 
  theme(axis.text.x = element_blank(), axis.title = element_text(size = 20), axis.text.y = element_text(size = 13), legend.text = element_text(size = 13), legend.title = element_text(size = 15)) + 
  scale_fill_gradient(low = "#FFFFFF", high = "#680000") +
  geom_vline(xintercept = c(58, 117, 177), alpha = 0.2)

#Adding labels
y.min <- 2.5
y.max <- 3.5
y.mid <- (y.min + y.max) / 2
heatplot <- heatplot +
  annotate("rect", xmin = 0.5, xmax = 5.5, ymin = y.min, ymax = y.max,
           alpha = 0.9, fill = col_island[1]) +
  annotate("rect", xmin = 5.5, xmax = 14.5, ymin = y.min, ymax = y.max,
           alpha = 0.95, fill = col_island[2]) +
  annotate("rect", xmin = 14.5, xmax = 34.5, ymin = y.min, ymax = y.max,
           alpha = 0.95, fill = col_island[3]) +
  annotate("rect", xmin = 34.5, xmax = 54.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[4]) + 
  annotate("rect", xmin = 54.5, xmax = 65.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[5]) + 
  annotate("rect", xmin = 65.5, xmax = 85.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[6]) + 
  annotate("rect", xmin = 85.5, xmax = 86.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[7]) + 
  annotate("rect", xmin = 86.5, xmax = 87.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[8]) + 
  annotate("rect", xmin = 87.5, xmax = 97.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[9]) + 
  annotate("text", x = 3, y = y.mid, label = names(col_island[1]), size = 7, srt = 90) +
  annotate("text", x = 10, y = y.mid, label = names(col_island[2]), size = 7, srt = 90) +
  annotate("text", x = 24.5, y = y.mid, label = names(col_island[3]), size = 7, srt = 90) +
  annotate("text", x = 44.5, y = y.mid, label = names(col_island[4]), size = 7, srt = 90) +
  annotate("text", x = 60, y = y.mid, label = names(col_island[5]), size = 7, srt = 90) +
  annotate("text", x = 75.5, y = y.mid, label = names(col_island[6]), size = 7, srt = 90) +
  annotate("text", x = 86, y = y.mid, label = names(col_island[7]), size = 4, srt = 90) +
  annotate("text", x = 87, y = y.mid, label = names(col_island[8]), size = 4, srt = 90) +
  annotate("text", x = 92.5, y = y.mid, label = names(col_island[9]), size = 7, srt = 90)
heatplot
ggsave(heatplot, filename = "./Output/Pl/Analysis/Heatmap of Phylum Grouped by Island.svg", dpi = 300, width = 18, height = 10)


## By Class
factor_melted_ps_by_sample <- unique(melted_ps$Sample) #factoring by sample
dataframe <- data.frame()  #creating dataframe to fill in information

for (i in factor_melted_ps_by_sample) {                                    #For each sample,
  sub_ps <- melted_ps[melted_ps$SampleID == i,]  #subset dataframe.
  
  per_class_abundance <- aggregate(Abundance ~ Class, data = sub_ps, sum)  #Aggregate abundance based on each class
  
  per_class_abundance$SampleID <- sub_ps$SampleID[1:nrow(per_class_abundance)] #Adding sample name to each class row
  per_class_abundance$Island <- sub_ps$Island[1:nrow(per_class_abundance)]
  
  dataframe <- rbind(dataframe, per_class_abundance)                       #store this in a dataframe for each row
}

#Sorting dataframe in order of structure, sub-ordered by location
ordered_df <- dataframe[order(dataframe$Island),]

#Plotting
table(ordered_df$Island) / sum(unique(ordered_df$Class) == unique(ordered_df$Class))
heatplot <- ggplot(ordered_df, aes(reorder(SampleID, -desc(Island)), reorder(Class, desc(Class)))) +
  geom_tile(aes(fill = Abundance)) + 
  labs(y = "Class", x = "Samples") + 
  theme(axis.text.x = element_blank(), axis.title = element_text(size = 20), axis.text.y = element_text(size = 10), legend.text = element_text(size = 9), legend.title = element_text(size = 10)) + 
  scale_fill_gradient(low = "#FFFFFF", high = "#680000") +
  geom_vline(xintercept = c(58, 117, 177), alpha = 0.2)

#Adding labels
y.min <- 4.5
y.max <- 6
y.mid <- (y.min + y.max) / 2
heatplot <- heatplot +
  annotate("rect", xmin = 0.5, xmax = 5.5, ymin = y.min, ymax = y.max,
           alpha = 0.9, fill = col_island[1]) +
  annotate("rect", xmin = 5.5, xmax = 14.5, ymin = y.min, ymax = y.max,
           alpha = 0.95, fill = col_island[2]) +
  annotate("rect", xmin = 14.5, xmax = 34.5, ymin = y.min, ymax = y.max,
           alpha = 0.95, fill = col_island[3]) +
  annotate("rect", xmin = 34.5, xmax = 54.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[4]) + 
  annotate("rect", xmin = 54.5, xmax = 65.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[5]) + 
  annotate("rect", xmin = 65.5, xmax = 85.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[6]) + 
  annotate("rect", xmin = 85.5, xmax = 86.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[7]) + 
  annotate("rect", xmin = 86.5, xmax = 87.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[8]) + 
  annotate("rect", xmin = 87.5, xmax = 97.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[9]) + 
  annotate("text", x = 3, y = y.mid, label = names(col_island[1]), size = 7, srt = 90) +
  annotate("text", x = 10, y = y.mid, label = names(col_island[2]), size = 7, srt = 90) +
  annotate("text", x = 24.5, y = y.mid, label = names(col_island[3]), size = 7, srt = 90) +
  annotate("text", x = 44.5, y = y.mid, label = names(col_island[4]), size = 7, srt = 90) +
  annotate("text", x = 60, y = y.mid, label = names(col_island[5]), size = 7, srt = 90) +
  annotate("text", x = 75.5, y = y.mid, label = names(col_island[6]), size = 7, srt = 90) +
  annotate("text", x = 86, y = y.mid, label = names(col_island[7]), size = 4, srt = 90) +
  annotate("text", x = 87, y = y.mid, label = names(col_island[8]), size = 4, srt = 90) +
  annotate("text", x = 92.5, y = y.mid, label = names(col_island[9]), size = 7, srt = 90)
heatplot
ggsave(heatplot, filename = "./Output/Pl/Analysis/Heatmap of Class Grouped by Island.svg", dpi = 300, width = 18, height = 10)


## By Species
factor_melted_ps_by_sample <- unique(melted_ps$Sample) #factoring by sample
dataframe <- data.frame()  #creating dataframe to fill in information

for (i in factor_melted_ps_by_sample) {                                    #For each sample,
  sub_ps <- melted_ps[melted_ps$SampleID == i,]  #subset dataframe.
  
  per_species_abundance <- aggregate(Abundance ~ Species, data = sub_ps, sum)  #Aggregate abundance based on each species
  
  per_species_abundance$SampleID <- sub_ps$SampleID[1:nrow(per_species_abundance)] #Adding sample name to each species row
  per_species_abundance$Island <- sub_ps$Island[1:nrow(per_species_abundance)]
  
  dataframe <- rbind(dataframe, per_species_abundance)                       #store this in a dataframe for each row
}

#Sorting dataframe in order of structure, sub-ordered by location
ordered_df <- dataframe[order(dataframe$Island),]

#Plotting
table(ordered_df$Island) / sum(unique(ordered_df$Species) == unique(ordered_df$Species))
genus.species.labels <- unique(taxtable.update$Genus_species)[rev(sorter)]
genus.species.labels <- genus.species.labels[-1]
heatplot <- ggplot(ordered_df, aes(reorder(SampleID, -desc(Island)), reorder(Species, desc(Species)))) +
  geom_tile(aes(fill = Abundance)) + 
  labs(y = "Species", x = "Samples") + 
  theme(axis.text.x = element_blank(), axis.title = element_text(size = 20), axis.text.y = element_text(size = 10), legend.text = element_text(size = 9), legend.title = element_text(size = 10)) + 
  scale_fill_gradient(low = "#FFFFFF", high = "#680000") +
  geom_vline(xintercept = c(58, 117, 177), alpha = 0.2) +
  scale_y_discrete(labels = genus.species.labels)

#Adding labels
y.min <- 21.5
y.max <- 29
y.mid <- (y.min + y.max) / 2
heatplot <- heatplot +
  annotate("rect", xmin = 0.5, xmax = 5.5, ymin = y.min, ymax = y.max,
           alpha = 0.9, fill = col_island[1]) +
  annotate("rect", xmin = 5.5, xmax = 14.5, ymin = y.min, ymax = y.max,
           alpha = 0.95, fill = col_island[2]) +
  annotate("rect", xmin = 14.5, xmax = 34.5, ymin = y.min, ymax = y.max,
           alpha = 0.95, fill = col_island[3]) +
  annotate("rect", xmin = 34.5, xmax = 54.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[4]) + 
  annotate("rect", xmin = 54.5, xmax = 65.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[5]) + 
  annotate("rect", xmin = 65.5, xmax = 85.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[6]) + 
  annotate("rect", xmin = 85.5, xmax = 86.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[7]) + 
  annotate("rect", xmin = 86.5, xmax = 87.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[8]) + 
  annotate("rect", xmin = 87.5, xmax = 97.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[9]) + 
  annotate("text", x = 3, y = y.mid, label = names(col_island[1]), size = 7, srt = 90) +
  annotate("text", x = 10, y = y.mid, label = names(col_island[2]), size = 7, srt = 90) +
  annotate("text", x = 24.5, y = y.mid, label = names(col_island[3]), size = 7, srt = 90) +
  annotate("text", x = 44.5, y = y.mid, label = names(col_island[4]), size = 7, srt = 90) +
  annotate("text", x = 60, y = y.mid, label = names(col_island[5]), size = 7, srt = 90) +
  annotate("text", x = 75.5, y = y.mid, label = names(col_island[6]), size = 7, srt = 90) +
  annotate("text", x = 86, y = y.mid, label = names(col_island[7]), size = 4, srt = 90) +
  annotate("text", x = 87, y = y.mid, label = names(col_island[8]), size = 4, srt = 90) +
  annotate("text", x = 92.5, y = y.mid, label = names(col_island[9]), size = 7, srt = 90)
heatplot
ggsave(heatplot, filename = "./Output/Pl/Analysis/Heatmap of Species Grouped by Island.svg", dpi = 300, width = 18, height = 10)




### Shannon Diversity Plots ###
div <- data.frame(Island = ps1ra@sam_data$Island,
                  Shannon = diversity(otu_table(ps1ra)))
Richness = colSums(decostand(otu_table(ps1ra), method = "pa"))
write.csv(div, file = "./Output/Pl/Analysis/Diversity_table.csv", quote = FALSE)
# By Island
div %>% group_by(Island) %>%
  summarise(N = n(), Mean = mean(Shannon))
ggplot(div, aes(x=Island, y=Shannon, fill = Island)) + 
  geom_boxplot() + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 25), axis.title.y = element_text(vjust = 2), axis.title = element_text(size = 30), title = element_text(size = 17), legend.text = element_text(size = 25), legend.title = element_text(size = 30)) + 
  labs(y="Shannon Diversity") +
  scale_fill_manual(values = col_island)
ggsave(filename = "./Output/Pl/Analysis/Shannon_Diversity_by_Location.svg", dpi = 300, width = 12, height = 10)




### Ordination(s) ###
ps.nmds <- subset_taxa(ps1ra, !is.na(Genus)) # Getting rid off non identified taxa at genus level
ps.nmds <- prune_samples(sample_sums(ps.nmds) > 0, ps.nmds)

NMDS = ordinate(ps.nmds, method = "NMDS", distance = "bray", trymax = 100)
PCoA = ordinate(ps.nmds, method = "PCoA", distance = "bray", trymax = 100)

# Stress Plot
svg("./Output/Pl/Analysis/Full_NMDS_Stress_Plot.svg")
p_stress_full <- stressplot(NMDS, title("Stress Plot for NMDS"))
dev.off()

# NMDS with/without ellipses
NMDS2 = data.frame(NMDS1 = NMDS$points[,1], NMDS2 = NMDS$points[,2],group=ps.nmds@sam_data$Island)
NMDS2.mean=aggregate(NMDS2[,1:2],list(group=ps.nmds@sam_data$Island),mean)

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell <- data.frame()
for(g in levels(NMDS2$group)){
  if(g != "St. John" & g != "Sultan Shoal")
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS2[NMDS2$group==g,],
                                                     veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2)))))
                                  ,group=g))
}
p_NMDS1 <- ggplot(data = NMDS2, aes(NMDS1, NMDS2)) + geom_point(aes(color = group), size = 4) +
  ggtitle(paste("NMDS (Stress Value = ", toString(round(NMDS$stress, digits = 3)), ")", sep = "")) + theme_bw() + scale_color_manual(values = col_island) + 
  theme(axis.title = element_text(size = 20), title = element_text(size = 20), 
        legend.text = element_text(size = 20)) + 
  guides(shape = guide_legend(override.aes = list(size = 5)), color = guide_legend(override.aes = list(size = 5))) + 
  labs(color = "Island")
p_NMDS1_ell <- ggplot(data = NMDS2, aes(NMDS1, NMDS2)) + geom_point(aes(color = group), size = 4) +
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,colour=group), size=2, linetype=2) +
  ggtitle(paste("NMDS (Stress Value = ", toString(round(NMDS$stress, digits = 3)), ")", sep = "")) + theme_bw() + scale_color_manual(values = col_island) + 
  theme(axis.title = element_text(size = 20), title = element_text(size = 20), 
        legend.text = element_text(size = 20)) + 
  guides(shape = guide_legend(override.aes = list(size = 5)), color = guide_legend(override.aes = list(size = 5))) + 
  labs(color = "Structure")
p_NMDS1_ell

# PCoA
p_PCoA <- plot_ordination(ps.nmds, PCoA, color = "Island") +  
  theme_bw() + scale_color_manual(values = col_island) + geom_point(size = 3) +
  theme(axis.title = element_text(size = 10), title = element_text(size = 12), 
        legend.text = element_text(size = 10)) + 
  guides(color = guide_legend(override.aes = list(size = 4)))

ggsave(p_NMDS1, filename = "./Output/Pl/Analysis/Full_NMDS_w_Island_colored.svg", dpi = 300, width = 12, height = 10)
ggsave(p_NMDS1_ell, filename = "./Output/Pl/Analysis/Full_NMDS_w_Island_colored_and_ellipses.svg", dpi = 300, width = 12, height = 10)
ggsave(p_PCoA, filename = "./Output/Pl/Analysis/Full_PCoA_w_Island_colored.svg", dpi = 300, width = 12, height = 10)




### PERMANOVA Test ###
ps.permanova <- subset_taxa(ps.nmds, !is.na(Genus))

otu = as.data.frame(otu_table(ps.permanova))
meta = as.data.frame(sample_data(ps.permanova))
df = data.frame(SampleID = meta$SampleID, Island = meta$Island)
# Island
permanova_Island <- adonis(otu ~ Island, data = df)

sink("./Output/Pl/Analysis/adonis_Island_table.txt")
noquote(print("PermANOVA Table:"))
permanova_Island
sink(NULL)

pairwise.adonis <- function(x,factors, sim.method = 'bray', p.adjust.m ='bonferroni')
{
  library(vegan)
  co = combn(unique(factors),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(co[1,elem],co[2,elem]),] ~ factors[factors %in% c(co[1,elem],co[2,elem])] , method =sim.method);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}
padonis_Island <- pairwise.adonis(otu,as.character(meta$Island))
sink("./Output/Pl/Analysis/adonis_Island_table.txt", append = TRUE)
noquote(print("Pairwise adonis between islands (Bonferroni corrected Pvalues):"))
padonis_Island
sink(NULL)




### Mantel Test ###
ps.mantel <- ps1ra

## Extracting Longitude and Latitude data
meta <- as.data.frame(ps.mantel@sam_data)
meta$lon <- sapply(strsplit(meta$GPS, " "), `[`, 2)
meta$lon <- as.double(sapply(strsplit(meta$lon, "E"), `[`, 1), length=10)
meta$lat <- as.double(sapply(strsplit(meta$GPS, "N"), `[`, 1), length=10)

geo_full <- data.frame(meta$lon, meta$lat)

## Preparing asv tables
otu_full <- as.data.frame(ps.mantel@otu_table)

## Making distance matrices
# Adundance data frames - bray curtis dissimilarity
dist.otu_full <- vegdist(otu_full, method = "bray")

# Geographic data frame - haversie distance
d.geo_full <- distm(geo_full, fun = distHaversine)

dist.geo_full <- as.dist(d.geo_full)

## Running Mantel Test
# Abundance vs Geographic
abund_geo_full <- vegan::mantel(dist.otu_full, dist.geo_full, method = "spearman", permutations = 9999, na.rm = TRUE)

## Saving Output
sink("./Output/Pl/Analysis/Mantel_Test.txt")
noquote("Mantel Test on All Samples")
abund_geo_full
sink(NULL)




### Multiple Regression on distance matrices ###
dist_MRM <- MRM(dist.otu_full ~ dist.geo_full,  nperm = 9999)
dist_MRM

sink("./Output/Pl/Analysis/MRM_Table.txt")
print("Bray-Curtis distance regressed against spatial distance (Multiple regression on matrices) (All Samples):")
print(dist_MRM)
sink(NULL)

