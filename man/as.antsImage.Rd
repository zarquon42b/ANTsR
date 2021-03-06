% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/antsImage_class.R
\name{as.antsImage}
\alias{as.antsImage}
\title{as.antsImage}
\usage{
as.antsImage(object, ...)
}
\arguments{
\item{object}{An object}

\item{...}{Extra named arguments passed to FUN}

\item{pixeltype}{a character string e.g. "float"}

\item{spacing}{numeric vector matching image dimensionality e.g. c(1.2,1.2)}

\item{origin}{numeric vector matching image dimensionality e.g. c(0,0)}

\item{direction}{numeric matrix matching image dimensionality eg diag(2)}

\item{components}{number of components per voxel, e.g. 1}

\item{reference}{optional reference antsImage providing all header info}
}
\description{
convert types to antsImage
}

