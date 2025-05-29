\name{multi_SPARKX}
\alias{multi_SPARKX}
\alias{multi_SPARKX.default}
\alias{multi_SPARKX.SpatialExperiment}
\title{Run SPARK-X Differential Expression Across Strata}
\description{
  Executes SPARK-X marker detection for each sample or stratum and aggregates the results.
}
\usage{
  multi_SPARKX(data, spatial, samples = NULL, n.cores = 1)
}
\arguments{
  \item{data}{
    Matrix of gene expression values (spots × features).
  }
  \item{spatial}{
    Spatial coordinates aligned to `data`.
  }
  \item{samples}{
    Optional factor for strata (e.g. sample IDs).
  }
  \item{n.cores}{
    Number of CPU cores to use.
  }
}
\value{
  Character vector of ranked feature names by significance.
}
\author{
  Stefano Cacciatore
}
\seealso{
  \code{\link{sparkx}}
}
\examples{
  markers <- multi_SPARKX(data_mat, coords, samples, n.cores = 2)
}
\keyword{methods}
