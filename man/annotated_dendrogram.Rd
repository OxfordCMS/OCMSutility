% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotated_dendrogram.R
\name{annotated_dendrogram}
\alias{annotated_dendrogram}
\title{annotated_dendrogram}
\usage{
annotated_dendrogram(
  dist,
  met,
  id,
  method = "complete",
  coord = NULL,
  pal = NULL
)
}
\arguments{
\item{dist}{class \code(dist). distance matrix}

\item{met}{dataframe. samples in rows, metadata variables in columns.}

\item{id}{string. sample identifier. must be a column in met}

\item{method}{string. agglomeration method used by \code(hclust).}

\item{coord}{numeric vector. optional. allows you to specify coordinates
of annotations. Vector must be same length as number of variables supplied.
Leave as \code(NULL) for default coordinates.}

\item{pal}{list. Colours to be used for each variable. see details section.
Set to \code(NULL) to use default (RColorBrewer) colours.}
}
\value{
ggplot
}
\description{
Creates an annotated dendrogram from a distance matrix.
Dendrogram annotations are based on metadata supplied
}
\details{
Colour values should be supplied \code(pal)
}
\examples{
set.seed(1)
# get relative abundance data
data(dss_example)
ddata <- dss_example$merged_abundance_id[,2:26]
rownames(ddata) <- dss_example$merged_abundance_id[,1]
ddata <- t(OCMSutility::relab(ddata))
# distance matrix
mydist <- vegan::vegdist(ddata, method='bray')

# metdata variable
mdata <- dss_example$metadata
mdata <- mdata[,c('sampleID','Genotype','Phenotype')]
annotated_dendrogram(mydist, mdata, 'sampleID')
# custom colours
col_geno <- RColorBrewer::brewer.pal(9, "Paired")[1:2]
col_phen <- RColorBrewer::brewer.pal(9, "Paired")[3:4]
annotated_dendrogram(mydist, mdata, 'sampleID', pal=list(col_geno, col_phen))
}