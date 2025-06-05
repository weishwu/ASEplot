#' Plot SNP-level ASE or parental allele frequency distribution as a heatmap for a given gene. 
#' @param ase_data ASE data frame, generated from ASEprep.
#' @param exons_all Exon coordinates, generated from ASEprep
#' @param gene_symbol Gene symbol.
#' @param data_type Either 'ase' (|0.5-RefFreq|), 'pat-freq' or 'mat-freq'.
#' @return A ggplot in which each column is a SNP site and eeach row is a sample.
#' @examples
#' data('ase_data.test')
#' ase_df = ase_data$ase_df
#' exons = ase_data$union_exons_per_gene
#' gene_ase_heatmap(ase_df, exons, "RHOBTB3", 'pat-freq')
#' @export
gene_ase_heatmap = function(ase_data, exons_all, gene_symbol, data_type) {

if (! data_type %in% c('ase', 'pat-freq', 'mat-freq')) {
stop(paste0(data_type, " is not 'ase', 'pat-freq', or 'mat-freq"))}

gene_symbol = as.character(gene_symbol)
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

snp_matrix = with(snp_data, tapply(snp_data$data_to_plot, list(as.character(snp_data$RNAid),snp_data$interpolated),sum))

htdata = expand.grid(X=as.character(colnames(snp_matrix)), Y=as.character(rownames(snp_matrix)))
htdata$data_to_plot = as.vector(t(snp_matrix))
htdata = merge(htdata, aggregate(rep(1,nrow(snp_data)),by=list(snp_data$RNAid),sum),by.x=2,by.y=1)
htdata$Y = factor(htdata$Y, levels=unique((htdata$Y)[order(htdata[,4])]))

if (data_type == 'ase') {
color_range = c(0, 0.5)
color_key = c("white", "deepskyblue")
leg_lab = 'ASE (|0.5 - RefFreq|)'
} else if (data_type == 'pat-freq') {
color_range = c(0, 1)
color_key = c("deeppink", "white", "deepskyblue")
leg_lab = 'Paternal Allele Freq'
} else if (data_type == 'mat-freq') {
color_range = c(0, 1)
color_key = c("deepskyblue", "white", "deeppink")
leg_lab = 'Maternal Allele Freq'}

ase_heatmap = ggplot(htdata, aes(X, Y)) + 
    geom_tile(aes(fill = data_to_plot)) + 
    scale_fill_gradientn(limits = color_range,colours = color_key,na.value = 'black') + 
    theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) + 
    labs(title=paste("sample count: ", length(unique(as.character(snp_data$RNAid))), "; SNP count: ", length(unique(as.numeric(snp_data$position))), sep=""), y = "RNA samples", x = "SNPs", fill = leg_lab)
return(ase_heatmap)
}

