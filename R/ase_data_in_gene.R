#' Subset ASE data by a gene name and create an interpolated exon coordinate file.
#' @param ase_data ASE data frame, generated from ASEprep
#' @param exons_all Exon coordinates, generated from ASEprep
#' @param gene_symbol Gene symbol
#' @return A list that contains the subsetted ASE data frame for the given gene, and the interpolated exon coordinates for plotting
#' @examples
#' data('ase_data.test')
#' ase_df = ase_data$ase_df
#' exons = ase_data$union_exons_per_gene
#' ase_data_in_gene(ase_df, exons, 'RHOBTB3')
#' @export
ase_data_in_gene = function(ase_data, exons_all, gene_symbol)
{
gene_symbol = as.character(gene_symbol)
#gene_data_allsamples = ase_data %>% filter(gene_name_from_exons == gene_symbol)
gene_data_allsamples = subset(ase_data, gene_name_from_exons == gene_symbol)
if (length(unique(as.character(gene_data_allsamples$gene_id_from_exons)))!=1) {
stop(paste0(gene_symbol," is not found or matches multiple gene IDs"))}

exons = as.data.frame(exons_all[exons_all[,6] == gene_symbol, , drop=F],stringsAsFactors = F)
exons[,2:3]=lapply(exons[,2:3], function(x) if(is.character(x)) as.numeric(x) else x)
exons = exons[order(exons[,2]),]
strand = unique(as.character(exons[,4]))
      
approxfun_list = as.list(1:nrow(exons))
exons_nospace = cbind(exons[,2:3],exons[,2:3])
for (i in 1:nrow(exons_nospace))
{
if (i>1)
{
exon_len=exons_nospace[i,2]-exons_nospace[i,1]
exons_nospace[i,3]=exons_nospace[i-1,4]+1
exons_nospace[i,4]=exons_nospace[i,3]+exon_len
}
approxfun_list[[i]]=approxfun(range(c(exons_nospace[i,1],exons_nospace[i,2])), c(exons_nospace[i,3],exons_nospace[i,4]))
}

positions = unique(c(as.numeric(gene_data_allsamples$position),as.numeric(exons[,2]),as.numeric(exons[,3])))

position_map = cbind(positions,positions)
colnames(position_map) = c('original','interpolated')

for (k in 1:nrow(position_map)) {
new_pos = c()
for (j in 1:length(approxfun_list)) {
new_pos = c(new_pos,approxfun_list[[j]](position_map[k,1]))}
position_map[k,2] = as.numeric(na.omit(new_pos))}

interpolate_fun = function(p){as.numeric(position_map[position_map[,1]==p,2])}
gene_data_allsamples$interpolated = unlist(lapply(as.numeric(gene_data_allsamples$position),interpolate_fun))

return(list(gene_data_allsamples,exons_nospace))
}

