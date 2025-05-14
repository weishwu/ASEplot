#' Plot SNP-level ASE or parental allele frequency distribution as a scatter plot for a given gene.
#' @param ase_data ASE data frame, generated from ASEprep.
#' @param exons_all Exon coordinates, generated from ASEprep
#' @param gene_symbol Gene symbol.
#' @param data_type Either 'ase' (|0.5-RefFreq|), 'pat-freq' or 'mat-freq'.
#' @param sample_name If a sample name is given, the data from this sample will be indicated.
#' @return A ggplot that shows each SNP-level ASE data point as a dot. Different exons are indicated by vertical lines. Introns are not shown.
#' @examples
#' data('ase_data.test')
#' ase_df = ase_data$ase_df
#' exons = ase_data$union_exons_per_gene
#' snp_gene_ase_scatter(ase_df, exons, "RHOBTB3", 'pat-freq', '123882')
#' @export
snp_gene_ase_scatter = function(ase_data, exons_all, gene_symbol, data_type, sample_name=NULL) {

if (! data_type %in% c('ase', 'pat-freq', 'mat-freq')) {
stop(paste0(data_type, " is not 'ase', 'pat-freq', or 'mat-freq"))}

snp_data_exons = ase_data_in_gene(ase_data, exons_all, gene_symbol)
snp_data = snp_data_exons[[1]]
exons_nospace = snp_data_exons[[2]]

gene_type = unique(matrix(unlist(strsplit(as.character(snp_data$exons_merged),':')),nrow(snp_data),7,byrow=T)[,7])
gene_strand = unique(matrix(unlist(strsplit(as.character(snp_data$exons_merged),':')),nrow(snp_data),7,byrow=T)[,4])

exon_blocks = data.frame(x1=c(exons_nospace[,3]), x2=exons_nospace[,4], y1=-0.01, y2=1.01, t=rep(c('a','b'),1e6)[1:nrow(exons_nospace)])

if (data_type == 'ase') {
snp_data$data_to_plot = snp_data$rawASE
} else if (data_type == 'pat-freq') {
snp_data = subset(snp_data, !is.na(PatFreq))
snp_data$data_to_plot = snp_data$PatFreq
} else if (data_type == 'mat-freq') {
snp_data = subset(snp_data, !is.na(PatFreq))
snp_data$data_to_plot = 1 - snp_data$PatFreq}

if (data_type == 'ase') {
hmarks = c(0, 0.5)
y_lab = 'ASE (|0.5 - RefFreq|)'
} else if (data_type == 'pat-freq') {
hmarks = c(0, 0.5, 1)
y_lab = 'Paternal Allele Freq'
} else if (data_type == 'mat-freq') {
hmarks = c(0, 0.5, 1)
y_lab = 'Maternal Allele Freq'}


if (!is.null(sample_name)) {

sample_name = as.character(sample_name)
if (! sample_name %in% snp_data$RNAid) {
stop(paste0(sample_name, " is not found in data"))}

snp_data = snp_data %>% mutate(
    samples = case_when(RNAid == sample_name ~ sample_name, TRUE ~ "others"))
snp_data$samples = factor(snp_data$samples, levels = c('others',sample_name))
snpCount = sum(snp_data$samples == sample_name)
median_val = format(round(median(snp_data$data_to_plot[snp_data$samples == sample_name], na.rm = T),2),nsmall=2)

scatterplot = ggplot(data = snp_data %>% filter(samples == 'others'), aes(x=interpolated, y=data_to_plot, col=samples, shape=samples, alpha=samples)) + 
   geom_point(size=2) + 
   geom_point(data=snp_data %>% filter(samples != 'others'), aes(x=interpolated, y=data_to_plot, col=samples, shape=samples, alpha=samples), inherit.aes = FALSE) +
   theme_classic() + 
   scale_color_manual(values=c("black", "red")) +
   scale_shape_manual(values=c(16, 17)) +
   scale_alpha_manual(values=c(0.2, 0.8)) +
   geom_vline(xintercept = exon_blocks$x1,col='gray60',alpha=0.2) +
   geom_hline(yintercept = hmarks, color = "gray60", alpha=0.2) +
   theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right") + 
   labs(title=paste("Gene symbol: ",gene_symbol, "; Gene type: ", gene_type, "; strand: ", gene_strand, "\nRNA sample: ", sample_name, "; SNP count: ", snpCount, "; median: ", median_val, sep=""), y = y_lab, x = 'exon boundaries (introns not shown)') +
   theme(axis.title.x = element_text(margin = margin(t = 0)))
return(scatterplot)

} else {

scatterplot = ggplot(data=snp_data, aes(x=interpolated, y=data_to_plot)) + 
   geom_point(col="black", shape=16, alpha=0.2, size=2) + 
   theme_classic() + 
   geom_vline(xintercept = exon_blocks$x1,col='gray60',alpha=0.2) +
   geom_hline(yintercept = hmarks, color = "gray60", alpha=0.2) +
   theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right") + 
   labs(title=paste0("Gene symbol: ",gene_symbol, "; Gene type: ", gene_type, "; strand: ", gene_strand), y = y_lab, x = 'exon boundaries (introns not shown)') +
   theme(axis.title.x = element_text(margin = margin(t = 0)))
return(scatterplot)
}
}

