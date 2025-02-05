library(tidyverse)
library(devtools)
library(qiime2R)

# NOTE: all .qza files used were accessed in R project. Path may vary if 
# code is replicated elsewhere 

# BRAY-CURTIS PCoA plots 
## Reading in data 
braycurtis <- read_qza("Full Sample_Alt_Filtering_take2/core-metrics-results/bray_curtis_pcoa_results.qza") 
metadata<-read_q2metadata("Full Sample_Alt_Filtering_take2/metadata-44.tsv")

## Grabbing PCs and merging with metadata to create one complete dataset
braycurtisPCOA <- braycurtis$data$Vectors %>%
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata)

## Mapping PCoA for sex and plate 
### Creating centroids for sexes
centroid_sex <- braycurtisPCOA %>% 
  group_by(Sex) %>% 
  summarize(PC1=mean(PC1), PC2=mean(PC2)) 

### Creating PCoA plot 
braycurtisPCOA %>% 
  ggplot(aes(x=PC1, y=PC2, color=Sex, shape=Plate)) +
  geom_point(size = 2.5) + 
  theme_q2r() + 
  scale_fill_manual(values=c("#0000fe", "#7da9d7")) + 
  scale_color_manual(values=c("#0000fe", "#7da9d7")) +
  theme(
  text = element_text(size = 14), # Adjusts all text size
  axis.title = element_text(size = 16), # Axis title font size
  axis.text = element_text(size = 14), # Axis text font size
  legend.title = element_text(size = 14), # Legend title font size
  legend.text = element_text(size = 12) # Legend text font size
)
## Mapping PCoA for county 
### Creating centroids for counties
centroid_county <- braycurtisPCOA %>% 
  group_by(County) %>% 
  summarize(PC1=mean(PC1), PC2=mean(PC2)) 

### Creating PCoA plot 
braycurtisPCOA %>% 
  ggplot(aes(x=PC1, y=PC2, color=County)) +
  geom_point(size=2.5) + 
  theme_q2r() + 
  scale_fill_manual(values=c("#fc8e59", "#b4b4b4","#99d594","#cc6766","#9999cc","#65cc99")) + 
  scale_color_manual(values=c("#fc8e59", "#b4b4b4","#99d594","#cc6766","#9999cc","#65cc99")) +
  theme(
    text = element_text(size = 14), # Adjusts all text size
    axis.title = element_text(size = 16), # Axis title font size
    axis.text = element_text(size = 14), # Axis text font size
    legend.title = element_text(size = 14), # Legend title font size
    legend.text = element_text(size = 12) # Legend text font size
  )

# Unweighted Unifrac 
## Reading in data 
unweighted_unifrac <- read_qza("Full Sample_Alt_Filtering_take2/core-metrics-results/unweighted_unifrac_pcoa_results.qza") 

## Grabbing PCs and merging with metadata to create one complete dataset
unweighted_unifrac_PCOA <- unweighted_unifrac$data$Vectors %>%
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata)

## Creating centroids for sexes
centroid_sex_uu <- unweighted_unifrac_PCOA %>% 
  group_by(Sex) %>% 
  summarize(PC1=mean(PC1), PC2=mean(PC2)) 

## Creating PCoA plots 
unweighted_unifrac_PCOA %>% 
  ggplot(aes(x=PC1, y=PC2, color=Sex, shape=Plate)) +
  geom_point() + 
  xlab("PC1 (28.48%)") + ylab("PC2 (9.87%)") +
  theme_q2r() + 
  scale_fill_manual(values=c("#0000fe", "#7da9d7")) + 
  scale_color_manual(values=c("#0000fe", "#7da9d7")) 

unweighted_unifrac_PCOA %>% 
  ggplot(aes(x=PC1, y=PC2, color=County)) +
  geom_point() + 
  theme_q2r() + 
  xlab("PC1 (28.48%)") + ylab("PC2 (9.87%)") +
  scale_fill_manual(values=c("#fc8e59", "#b4b4b4","#99d594","#cc6766","#9999cc","#65cc99")) + 
  scale_color_manual(values=c("#fc8e59", "#b4b4b4","#99d594","#cc6766","#9999cc","#65cc99"))

# Weighted Unifrac 
## Reading in data 
weighted_unifrac <- read_qza("Full Sample_Alt_Filtering_take2/core-metrics-results/weighted_unifrac_pcoa_results.qza") 

## Grabbing PCs and merging with metadata to create one complete dataset
weighted_unifrac_PCOA <- weighted_unifrac$data$Vectors %>%
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata)

## Creating centroids for sexes
centroid_sex_wu <- weighted_unifrac_PCOA %>% 
  group_by(Sex) %>% 
  summarize(PC1=mean(PC1), PC2=mean(PC2)) 

## Creating PCoA plots 
weighted_unifrac_PCOA %>% 
  ggplot(aes(x=PC1, y=PC2, color=Sex, shape=Plate)) +
  geom_point() + 
  xlab("PC1 (53.58%)") + ylab("PC2 (16.15%)") +
  theme_q2r() + 
  theme_q2r() + 
  scale_fill_manual(values=c("#0000fe", "#7da9d7")) + 
  scale_color_manual(values=c("#0000fe", "#7da9d7")) 

weighted_unifrac_PCOA %>% 
  ggplot(aes(x=PC1, y=PC2, color=County)) +
  geom_point() + 
  theme_q2r() + 
  xlab("PC1 (53.58%)") + ylab("PC2 (16.15%)") +
  scale_fill_manual(values=c("#fc8e59", "#b4b4b4","#99d594","#cc6766","#9999cc","#65cc99")) + 
  scale_color_manual(values=c("#fc8e59", "#b4b4b4","#99d594","#cc6766","#9999cc","#65cc99")) 
