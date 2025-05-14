#' Plot gene-level parental allele frequency distribution as histogram.
#' @param ase_data ASE data frame, generated from ASEprep.
#' @param gene_symbol Gene symbol.
#' @param data_type Either 'pat-freq' or 'mat-freq'.
#' @param binwidth Histogram binwidth. Defaults to 0.01.
#' @param sample_name If a sample name is provided, a vertical line will be drawn to indicate its position relative to the histogram.
#' @return A ggplot.
#' @examples
#' data('ase_data.test')
#' ase_df = ase_data$ase_df
#' gene_poe_histogram(ase_df, 'RHOBTB3', 'pat-freq', '123882')
#' @export
gene_poe_histogram = function(ase_data, gene_symbol, data_type, binwidth = NULL, sample_name = NULL) {

if (! data_type %in% c('pat-freq', 'mat-freq')) {
stop(paste0(data_type, " is not 'pat-freq', or 'mat-freq"))}

if (is.null(binwidth)) {binwidth = 0.01}

gene_symbol = as.character(gene_symbol)
ase_data$genes_exonic_symbol = sapply(ase_data$genes_exonic, get_gene_symbol)
snp_data = subset(ase_data, ((genes_exonic_symbol == gene_symbol) & (!is.na(PatFreq))))

if (length(unique(as.character(snp_data$genes_exonic)))!=1) {
stop(paste0(gene_symbol," is not found or matches multiple gene IDs"))}

counts_bygene = aggregate(subset(snp_data, select=c(PatDepth,MatDepth,totalCount)),by=list(snp_data$RNAid,snp_data$genes_exonic),sum)
colnames(counts_bygene)[1:2] = c('RNAid','gene_exonic')
counts_bygene$PatFreq_PerGenePerRNAid = counts_bygene$PatDepth/counts_bygene$totalCount

if (data_type == 'pat-freq') {
counts_bygene$data_to_plot = counts_bygene$PatFreq
x_lab = 'Gene level PatFreq'
} else if (data_type == 'mat-freq') {
counts_bygene$data_to_plot = 1 - counts_bygene$PatFreq
x_lab = 'Gene level MatFreq'}

# histogram
if (!is.null(sample_name)) {
poe_histogram = ggplot(counts_bygene,aes(x= data_to_plot)) +
geom_histogram(binwidth = binwidth) +
coord_cartesian(xlim = c(0,1)) +
theme_bw() +
geom_vline(data=counts_bygene[counts_bygene$RNAid == sample_name,], aes(xintercept = data_to_plot,col=RNAid), linewidth=1, linetype='dashed') + 
xlab(x_lab) + ylab('Counts') + ggtitle(gene_symbol)
} else {

poe_histogram = ggplot(counts_bygene,aes(x= data_to_plot)) +
geom_histogram(binwidth=0.01) +
coord_cartesian(xlim = c(0,1)) +
theme_bw() + 
xlab(x_lab) + ylab('Counts') + ggtitle(gene_symbol)
}
return(poe_histogram)
}

