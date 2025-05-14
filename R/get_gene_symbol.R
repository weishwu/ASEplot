#' Get gene symbos from the genes_exonic column.
#' @param ase_data ASE data frame, generated from ASEprep
#' @return A vector with the gene symbols split from the genes_exonic column
#' @examples
#' data('ase_data.test')
#' ase_df = ase_data$ase_df
#' ase_df$genes_exonic_symbol = sapply(ase_df$genes_exonic, get_gene_symbol)
#' @export
get_gene_symbol = function(x) {
    if (is.na(x)) {
        return(NA)
    } else {
        x_split = strsplit(x, split = ";")[[1]]
        return(paste(sapply(x_split, function(i){strsplit(i,':')[[1]][2]}), collapse=';'))
    }}
