\name{refine_SVM}
\alias{refine_SVM}
\title{Support Vector Machine–Based Spatial Label Refinement}
\description{
  Performs SVM-based label refinement across spatial neighborhoods or tiles.
}
\usage{
  refine_SVM(xy, labels, samples = NULL, newdata = NULL, newsamples = NULL, tiles = NULL, ...)
}
\arguments{
  \item{xy}{
    Spatial coordinates matrix.
  }
  \item{labels}{
    Vector of original class labels.
  }
  \item{samples}{
    Stratum vector.
  }
  \item{newdata}{
    New data to predict refined labels.
  }
  \item{newsamples}{
    Stratum vector for new data.
  }
  \item{tiles}{
    Grid configuration for tiling.
  }
  \item{...}{
    Arguments passed to svm.
  }
}
\value{
  Factor of refined labels.
}
\author{
  Stefano Cacciatore
}
\seealso{
  \code{\link[e1071]{svm}}
}
\examples{
  refined <- refine_SVM(coords, labels, tiles = c(10, 10))
}
\keyword{methods}
