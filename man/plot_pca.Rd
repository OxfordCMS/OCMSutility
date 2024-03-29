% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_pca.R
\name{plot_pca}
\alias{plot_pca}
\title{plot_pca}
\usage{
plot_pca(
  pc,
  metadata,
  colourby = "none",
  shapeby = "none",
  group = "none",
  continuous = FALSE,
  pcs = c("PC1", "PC2")
)
}
\arguments{
\item{pc}{prcomp object}

\item{metadata}{metadata dataframe}

\item{colourby}{colour by. default "none"}

\item{shapeby}{shape by. default "none"}

\item{group}{group aes in ggplot. default "none"}

\item{continuous}{is colourby a continuous variable? default FALSE}

\item{pcs}{principal components to plot. default \code{c("PC1", "PC2")}}
}
\value{
list of ggplot pca plot and corresponding dataframe used to make the plot
}
\description{
plot PCA
}
\examples{
pca_result <- prcomp(USArrests, scale = TRUE)
state_data <- data.frame(abb = state.abb, area = state.area,
                         center = state.center, region = state.region,
                         division = state.division)
rownames(state_data) <- state.name
plot_pca(pca_result, state_data, colourby='division', shapeby='region')
}
