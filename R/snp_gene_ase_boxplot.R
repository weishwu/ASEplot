#' Plot SNP-level ASE, or parental allele frequency distribution as a boxplot for a given gene for all samples.
#' @param ase_data ASE data frame, generated from ASEprep.
#' @param gene_symbol Gene symbol.
#' @param data_type Either 'ase' (|0.5-RefFreq|), 'pat-freq' or 'mat-freq'.
#' @return A ggplot.
#' @examples
#' data('ase_data.test')
#' ase_df = ase_data$ase_df
#' snp_gene_ase_boxplot(ase_df, 'RHOBTB3', 'pat-freq')
#' @export
snp_gene_ase_boxplot = function(ase_data, gene_symbol, data_type) {

gene_symbol = as.character(gene_symbol)

ase_data$genes_exonic_symbol = sapply(ase_data$genes_exonic, get_gene_symbol)

snp_data = subset(ase_data, (genes_exonic_symbol == gene_symbol))

if (length(unique(as.character(snp_data$genes_exonic)))!=1) {
stop(paste0(gene_symbol," is not found or matches multiple gene IDs"))}

if (data_type == 'ase') {
snp_data$data_to_plot = snp_data$rawASE
} else if (data_type == 'pat-freq') {
snp_data = subset(snp_data, !is.na(PatFreq))
snp_data$data_to_plot = snp_data$PatFreq
} else if (data_type == 'mat-freq') {
snp_data = subset(snp_data, !is.na(PatFreq))
snp_data$data_to_plot = 1 - snp_data$PatFreq}

agg_ = aggregate(snp_data$data_to_plot,by=list(snp_data$RNAid),median)
agg_ = cbind(agg_,order(order(agg_[,2])))
colnames(agg_)=c('RNAid','median_data','median_order')
snp_data = snp_data %>% inner_join(agg_,by='RNAid') %>% mutate(SNP_position = paste(contig, position, strand,sep='_')) %>% mutate(RNAid_2 = as.character(RNAid))
snp_data[,'SNP_position']=gsub('_minus', '', gsub('_plus', '', snp_data[,'SNP_position']))

x_axis_text_size = 10
if (nrow(agg_) > 50) {x_axis_text_size=7}

if (data_type == 'ase') {
y_lim = c(0, 0.5)
y_lab = 'ASE (|0.5 - RefFreq|)'
} else if (data_type == 'pat-freq') {
y_lim = c(0, 1)
y_lab = 'Paternal Allele Freq'
} else if (data_type == 'mat-freq') {
y_lim = c(0, 1)
y_lab = 'Maternal Allele Freq'}

ase_boxplot = ggplot(snp_data, aes(x=reorder(RNAid_2, data_to_plot, FUN=median), y=data_to_plot)) + 
geom_boxplot() + 
theme_classic() +
theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, size=x_axis_text_size),
    axis.text.y=element_text(size=10),
      axis.title=element_text(size=20),
      plot.title=element_text(size=20)) +
xlab('') + ylab(y_lab) + ggtitle(as.character(gene_symbol)) +
coord_cartesian(ylim = y_lim) +
geom_point(data=snp_data, aes(x=median_order, y=data_to_plot, color=SNP_position)) +
geom_hline(yintercept=0.5, linetype='dashed')

return(ase_boxplot)
}

