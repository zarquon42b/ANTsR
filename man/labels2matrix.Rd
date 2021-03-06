% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/labels2matrix.R
\name{labels2matrix}
\alias{labels2matrix}
\title{Convert label image to a matrix}
\usage{
labels2matrix(img, mask)
}
\arguments{
\item{img}{input label image}

\item{mask}{defines domain of interest}
}
\value{
matrix is output
}
\description{
Convert a labeled image to an n x m binary matrix where n = number of voxels
and m = number of labels. Only includes values inside the provided mask while
including background ( img == 0 ) for consistency with timeseries2matrix and
other image to matrix operations.
}
\examples{

fi = antsImageRead(getANTsRData("r16") ,2) \%>\% resampleImage(c(60,60),1,0)
mask = getMask( fi )
labs = kmeansSegmentation( fi, 3 )$segmentation
labmat = labels2matrix( labs, mask )

}
\author{
Avants BB
}

