\name{refinecluster}
\alias{refinecluster}
\alias{refinecluster.default}
\alias{refinecluster.Seurat}
\alias{refinecluster.SpatialExperiment}
\alias{refinecluster.giotto}
\title{Spatially Informed Cluster Refinement}
\description{
  Refines clusters by smoothing labels using neighborhood information.
}
\usage{
  refinecluster(object, shape = 'square', ...)
}
\arguments{
  \item{object}{
    Spatial object or cluster labels and coordinates.
  }
  \item{shape}{
    Neighborhood shape: 'square', 'circle', etc.
  }
  \item{...}{
    Further parameters to refinement logic.
  }
}
\value{
  Refined cluster labels or object with updated labels.
}
\author{
  Stefano Cacciatore
}
\seealso{
  \code{\link{refine_cluster}}
}
\examples{
  refined <- refinecluster(object, shape = 'circle')
}
\keyword{methods}
