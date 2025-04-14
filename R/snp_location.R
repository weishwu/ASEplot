#' Plot SNP positions relative to the given gene.
#' @param ase_data ASE data frame, generated from ASEprep.
#' @param exons_all Exon coordinates, generated from ASEprep
#' @param gene_symbol Gene symbol.
#' @param transcripts_display Either 'split' (keeps all the transcripts displayed separately) or 'collapse' (collapses the transcripts into a single model).
#' @param sample_name If a sample name is given, the SNP positions from this sample will be indicated.
#' @return A plot that shows the positions of SNPs in the top track (or two tracks, when a sample_name is given).
#' @examples
#' data('ase_data.test')
#' ase_df = ase_data$ase_df
#' exons = ase_data$union_exons_per_gene
#' snp_location(ase_df, exons, "RHOBTB3", 'collapse', '123882')
#' @importFrom GenomicRanges GRanges
#' @importFrom Gviz AnnotationTrack plotTracks GenomeAxisTrack BiomartGeneRegionTrack
#' @importFrom biomaRt useEnsembl
#' @export
snp_location = function(ase_data, exons_all, gene_symbol, transcripts_display, sample_name = NULL) {

snp_data_exons = ase_data_in_gene(ase_data, exons_all, gene_symbol)
snp_data = snp_data_exons[[1]]
exons_nospace = snp_data_exons[[2]]
exons=as.data.frame(exons_all[exons_all[,6]==gene_symbol, , drop=F],stringsAsFactors = F)
exons[,2:3]=lapply(exons[,2:3], function(x) if(is.character(x)) as.numeric(x) else x)
exons = exons[order(exons[,2]),]

snps_loc_all=unique(subset(snp_data,select=c(contig,position)))
snps_loc_all=snps_loc_all[order(snps_loc_all$position),]
snps_grange_all=GRanges(seqnames = snps_loc_all$contig,ranges = IRanges(start = snps_loc_all$position, end = snps_loc_all$position))
snps_all_track=AnnotationTrack(snps_grange_all, name = "SNPs_all", stacking="dense")

if (!is.null(sample_name)) {
    
if (! sample_name %in% snp_data$RNAid) {
stop(paste0(sample_name, "is not found in data"))}

snps_loc_sample=unique(snp_data %>% filter(RNAid == sample_name) %>% select(c(contig,position)))
snps_loc_sample=snps_loc_sample[order(snps_loc_sample$position),]
snps_grange_sample=GRanges(seqnames = snps_loc_sample$contig,ranges = IRanges(start = snps_loc_sample$position, end = snps_loc_sample$position))
snps_sample_track=AnnotationTrack(snps_grange_sample, name = paste("SNPs_", sample_name, sep=""), stacking="dense")

gtrack=GenomeAxisTrack()
mart=useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
biomTrack=BiomartGeneRegionTrack(start=min(exons[,2])-1000, end=max(exons[,3])+1000, chromosome=unique(exons[,1]), genome = "hg38", name = "ENSEMBL", biomart = mart, transcriptAnnotation = "symbol")

if (transcripts_display == 'split') {
plotTracks(list(gtrack,snps_all_track,snps_sample_track,biomTrack),from=min(exons[,2])-1000,to=max(exons[,3])+1000)
} else if (transcripts_display == 'collapse') {
plotTracks(list(gtrack,snps_all_track,snps_sample_track,biomTrack), collapseTranscripts="meta",from=min(exons[,2])-1000,to=max(exons[,3])+1000)}

} else {

gtrack=GenomeAxisTrack()
mart=useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
biomTrack=BiomartGeneRegionTrack(start=min(exons[,2])-1000, end=max(exons[,3])+1000, chromosome=unique(exons[,1]), genome = "hg38", name = "ENSEMBL", biomart = mart, transcriptAnnotation = "symbol")

if (transcripts_display == 'split') {
plotTracks(list(gtrack,snps_all_track,biomTrack),from=min(exons[,2])-1000,to=max(exons[,3])+1000)
} else if (transcripts_display == 'collapse') {
plotTracks(list(gtrack,snps_all_track,biomTrack), collapseTranscripts="meta",from=min(exons[,2])-1000,to=max(exons[,3])+1000)}
}
}

