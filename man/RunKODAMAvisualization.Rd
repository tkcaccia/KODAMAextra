\name{RunKODAMAvisualization}
\alias{RunKODAMAvisualization}
\alias{RunKODAMAvisualization.default}
\alias{RunKODAMAvisualization.SingleCellExperiment}
\alias{RunKODAMAvisualization.Seurat}
\alias{RunKODAMAvisualization.SpatialExperiment}
\alias{RunKODAMAvisualization.giotto}
\title{Generic KODAMA Visualization Method}
\description{
  Embeds KODAMA dissimilarity into a low-dimensional space using class-specific visualization methods.
}
\usage{
  RunKODAMAvisualization(object, ...)
}
\arguments{
  \item{object}{
    An object containing KODAMA dissimilarity results.
  }
  \item{...}{
    Additional visualization arguments.
  }
}
\value{
  Modified object with 2D embedding stored.
}
\author{
  Stefano Cacciatore
}
\seealso{
  \code{\link{KODAMA.matrix.parallel}}, \code{\link{RunKODAMAmatrix}}
}
\examples{
  sce <- RunKODAMAvisualization(sce, dims = 2)
}
\keyword{methods}
