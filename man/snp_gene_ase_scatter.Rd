% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/snp_gene_ase_scatter.R
\name{snp_gene_ase_scatter}
\alias{snp_gene_ase_scatter}
\title{Plot SNP-level ASE or parental allele frequency distribution as a scatter plot for a given gene.}
\usage{
snp_gene_ase_scatter(
  ase_data,
  exons_all,
  gene_symbol,
  data_type,
  sample_name = NULL
)
}
\arguments{
\item{ase_data}{ASE data frame, generated from ASEprep.}

\item{exons_all}{Exon coordinates, generated from ASEprep}

\item{gene_symbol}{Gene symbol.}

\item{data_type}{Either 'ase' (|0.5-RefFreq|), 'pat-freq' or 'mat-freq'.}

\item{sample_name}{If a sample name is given, the data from this sample will be indicated.}
}
\value{
A ggplot that shows each SNP-level ASE data point as a dot. Different exons are indicated by vertical lines. Introns are not shown.
}
\description{
Plot SNP-level ASE or parental allele frequency distribution as a scatter plot for a given gene.
}
\examples{
data('ase_data.test')
ase_df = ase_data$ase_df
exons = ase_data$union_exons_per_gene
snp_gene_ase_scatter(ase_df, exons, "RHOBTB3", 'pat-freq', '123882')
}
