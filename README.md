# ASEplot: Allele Specific Expression data plot

## Load data and filter
```
library(ggplot2)
library(Gviz)
library(GenomicRanges)
library(biomaRt)
library(ggridges)
library(dplyr)
library(tidyverse)
library(pheatmap)
library(ggrepel)
library(ASEplot)

data('ase_data.test')
ase_df = ase_data$ase_df
exons = ase_data$union_exons_per_gene

# Select only the lines with unique genes
ase_df_uniqGene = ase_df %>% filter(! grepl(';', genes_exonic))

# Filter data
ase_selc = ase_df_uniqGene %>% filter( 
    (totalCount >= 10) & 
    (!is.na(exons_merged)) & 
    (!grepl(';', exons_merged)) &
    (nonAltFreq_perRNAid < 0.05)) %>% dplyr::select(
    RNAid,variantID,contig,position,strand,refCount,altCount,
    totalCount,rawASE,exons_merged,gene_type_exonic,
    PatAllele,MatAllele,PatDepth,MatDepth,PatFreq,genes_exonic)
```

## Plots

### Check contamination

- Check the distribution of sample contamination measured by the non-Alt-Freq at 1/1 sites, and non-Ref-Freq at 0/0 sites

```
contam = unique(ase_df_uniqGene %>% select(RNAid, nonAltFreq_perRNAid, nonRefFreq_perRNAid))
ggplot(contam, aes(x=nonAltFreq_perRNAid, y=nonRefFreq_perRNAid)) + 
   geom_point(alpha=0.6) + 
   theme_bw() + 
   geom_vline(xintercept = 0.05, color='red', linetype='dashed',alpha=0.6) + 
   geom_text_repel(data = contam %>% filter(nonAltFreq_perRNAid > 0.05), aes(label = RNAid))
```
![](figures/contam.png)


- Check contamination per sample and per gene

```
gene_contam = unique(ase_df_uniqGene %>% 
   filter(! is.na(homRef_nonRefFreq_atMatAlt_mean_perGene_perRNAid)) %>% 
   select(RNAid, genes_exonic,homRef_nonRefFreq_atMatAlt_mean_perGene_perRNAid))

gene_contam = gene_contam %>% 
   pivot_wider(id_cols = genes_exonic, 
               names_from = RNAid, 
               values_from = homRef_nonRefFreq_atMatAlt_mean_perGene_perRNAid) %>% 
   column_to_rownames('genes_exonic')

pheatmap(gene_contam[1:40,], 
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         na_col ='white',color = colorRampPalette(c("skyblue", "red"))(500),
         breaks = seq(0, 0.05, 0.05/500))
```
![](figures/contam_per_gene.png)


### SNP location relative to a given gene

- With transcripts split
```
snp_location(ase_selc, exons, 'RHOBTB3', 'split', 'hg38', '123884')
```
![](figures/snp_location.png)

- With transcripts collapsed
```
snp_location(ase_selc, exons, 'RHOBTB3', 'collapse', 'hg38', '123884')
```
![](figures/snp_location_collapsed.png)


### Gene-level average POE (Parent-Of-Origin) ASE for a given gene across samples
```
gene_poe_histogram(ase_selc, 'RHOBTB3', 'pat-freq', sample_name = '123884')
```
![](figures/histogram.png)


### Gene-level average POE ridge plot comparing multiple genes
```
gene_poe_ridge(ase_selc, c('MEG8', 'CYB5R2', 'IGF2', 'RHOBTB3', 'THEGL', 'GNAS', 'PEG3'), 'pat-freq')
```
![](figures/ridges.png)


### SNP-level ASE for a given gene in a heatmap
```
snp_gene_ase_heatmap(ase_selc, exons, 'RHOBTB3', 'pat-freq')
```
![](figures/heatmap.png)


### SNP-level ASE for a given gene in a box plot
```
snp_gene_ase_boxplot(ase_selc, 'RHOBTB3', 'pat-freq')
```
![](figures/boxplot.png)

### SNP-level ASE for a given gene in a scatter plot
```
snp_gene_ase_scatter(ase_selc, exons, 'RHOBTB3','pat-freq','123884')
```
![](figures/scatter.png)


