\name{KODAMA.matrix.parallel}
\alias{KODAMA.matrix.parallel}
\title{Compute KODAMA Dissimilarity Matrix in Parallel}
\description{
  Executes the KODAMA algorithm with multiple permutations in parallel, optionally incorporating spatial constraints.
}
\usage{
  KODAMA.matrix.parallel(data,
    spatial = NULL,
    samples = NULL,
    M = 100,
    Tcycle = 20,
    FUN = c("fastpls", "simpls"),
    ncomp = 50,
    W = NULL,
    metrics = "euclidean",
    constrain = NULL,
    fix = NULL,
    landmarks = 10000,
    splitting = ifelse(nrow(data) < 40000, 100, 300),
    spatial.resolution = 0.3,
    n.cores = 1,
    seed = 1234)
}
\arguments{
  \item{data}{
    A numeric matrix of features (rows = samples, columns = features).
  }
  \item{spatial}{
    Optional spatial coordinates matrix aligned to `data`.
  }
  \item{samples}{
    Optional grouping vector defining strata for batch-specific processing.
  }
  \item{M}{
    Number of permutations (default 100).
  }
  \item{Tcycle}{
    Number of cycles per permutation to tune PLS components.
  }
  \item{FUN}{
    PLS implementation: either "fastpls" or "simpls".
  }
  \item{ncomp}{
    Number of PLS components (default 50).
  }
  \item{W}{
    Optional weights or labels for semi-supervised learning.
  }
  \item{metrics}{
    Distance metric for nearest neighbor search, e.g., "euclidean".
  }
  \item{constrain}{
    Optional grouping constraints for post-processing.
  }
  \item{fix}{
    Logical vector marking fixed samples.
  }
  \item{landmarks}{
    Number of landmark points to sample per permutation.
  }
  \item{splitting}{
    Number of clusters used when landmarking.
  }
  \item{spatial.resolution}{
    Proportion of spatially informed landmarks.
  }
  \item{n.cores}{
    Number of CPU cores to use (parallel execution).
  }
  \item{seed}{
    Random seed for reproducibility.
  }
}
\value{
  A list with dissimilarity matrix, permutation labels, accuracy trajectories, and knn indices.
}
\author{
  Stefano Cacciatore
}
\seealso{
  \code{\link{RunKODAMAmatrix}}
}
\examples{
  data_mat <- matrix(rnorm(2000), nrow = 100, ncol = 20)
  res <- KODAMA.matrix.parallel(data_mat, M = 50, n.cores = 2)
}
\keyword{methods}
