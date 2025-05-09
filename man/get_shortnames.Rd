% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_shortnames.R
\name{get_shortnames}
\alias{get_shortnames}
\title{get_shortnames
return the shortnames of taxonomy labels}
\usage{
get_shortnames(longnames)
}
\arguments{
\item{longnames}{vector; vector of long taxonomic labels}
}
\value{
vector
}
\description{
get_shortnames
return the shortnames of taxonomy labels
}
\examples{
longnames <- c("k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Ruminococcus",
               "k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides")
get_shortnames(longnames)
}
