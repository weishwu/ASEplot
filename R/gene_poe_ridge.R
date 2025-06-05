#' Plot gene-level parental allele frequency distribution as a ridge plot for given genes.
#' @param ase_data ASE data frame, generated from ASEprep.
#' @param gene_symbols Gene symbols.
#' @param data_type Either 'pat-freq' or 'mat-freq'.
#' @return A ggplot.
#' @examples
#' data('ase_data.test')
#' ase_df = ase_data$ase_df
#' gene_poe_ridge(ase_df, c("MEG3", "PEG3", "RHOBTB3"), 'pat-freq')
#' @export
gene_poe_ridge = function(ase_data, gene_symbols, data_type)
{

if (! data_type %in% c('pat-freq', 'mat-freq')) {
stop(paste0(data_type, " is not 'pat-freq', or 'mat-freq"))}

gene_symbols = as.character(gene_symbols)
ase_data$genes_exonic_symbol = sapply(ase_data$genes_exonic, get_gene_symbol)
ase_data = subset(ase_data, ((genes_exonic_symbol %in% gene_symbols) & (!is.na(PatFreq))))

if (length(unique(as.character(ase_data$genes_exonic))) == 0) {
stop(paste0("None of the genes is not found"))}

counts_bygene = aggregate(subset(ase_data, select=c(PatDepth,MatDepth,totalCount)),by=list(ase_data$RNAid,ase_data$genes_exonic_symbol),sum)
colnames(counts_bygene)[1:2] = c('RNAid','gene_exonic')
counts_bygene$PatFreq_PerGenePerRNAid = counts_bygene$PatDepth/counts_bygene$totalCount

if (data_type == 'pat-freq') {
counts_bygene$data_to_plot = counts_bygene$PatFreq
x_lab = 'Gene level PatFreq'
} else if (data_type == 'mat-freq') {
counts_bygene$data_to_plot = 1 - counts_bygene$PatFreq
x_lab = 'Gene level MatFreq'}

poe_ridge = ggplot(counts_bygene,aes(x= data_to_plot,y=gene_exonic,fill=after_stat(x))) +
geom_density_ridges_gradient(scale = 1, rel_min_height = 0.01) +
xlim(0,1) +
theme_classic() +
xlab(x_lab) + ylab('Gene') +
scale_fill_gradientn(name = x_lab, limits = c(0,1),colours=c("deeppink", "white", "deepskyblue")) +
theme(plot.margin=unit(c(0.5,0,0.5,0), "cm"), axis.text=element_text(size=10),
      axis.title=element_text(size=15),
      plot.title=element_text(size=20))  
return(poe_ridge)                 
}
 
