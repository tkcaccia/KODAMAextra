\name{RunKODAMAmatrix}
\alias{RunKODAMAmatrix}
\alias{RunKODAMAmatrix.default}
\alias{RunKODAMAmatrix.SingleCellExperiment}
\alias{RunKODAMAmatrix.Seurat}
\alias{RunKODAMAmatrix.SpatialExperiment}
\alias{RunKODAMAmatrix.giotto}
\title{Generic KODAMA Feature Extraction for Spatial Omics Data}
\description{
  Dispatches KODAMA dimensionality reduction on spatial omics data across supported object types.
}
\usage{
  RunKODAMAmatrix(object, reduction = NULL, dims = NULL, assay = NULL, ...)
}
\arguments{
  \item{object}{
    A spatial omics object: SingleCellExperiment, Seurat, SpatialExperiment, or Giotto.
  }
  \item{reduction}{
    Character. Name of the dimensionality reduction (e.g., 'PCA').
  }
  \item{dims}{
    Integer vector. Dimensions of the reduction to use.
  }
  \item{assay}{
    Character. Assay to extract if no reduction is specified.
  }
  \item{...}{
    Additional arguments passed to KODAMA.matrix.parallel or class-specific methods.
  }
}
\value{
  Modified object with KODAMA coordinates added to the appropriate slot.
}
\author{
  Stefano Cacciatore
}
\seealso{
  \code{\link{KODAMA.matrix.parallel}}
}
\examples{
  sce <- RunKODAMAmatrix(sce, reduction = 'PCA', dims = 1:30)
}
\keyword{methods}
