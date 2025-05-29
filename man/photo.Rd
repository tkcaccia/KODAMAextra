\name{photo}
\alias{photo}
\title{Extract Local KODAMA Feature Image Tiles}
\description{
  Creates fixed-size image tiles of KODAMA embeddings based on local spatial neighborhoods.
}
\usage{
  photo(vis, xy, range = 0.05, n_pixels = 25)
}
\arguments{
  \item{vis}{
    KODAMA embeddings (rows = spots, cols ≥ 2).
  }
  \item{xy}{
    Spatial coordinates aligned to `vis`.
  }
  \item{range}{
    Neighborhood window (default 0.05).
  }
  \item{n_pixels}{
    Number of pixels per image side (default 25).
  }
}
\value{
  Matrix where each row is a flattened image of local embedding.
}
\author{
  Stefano Cacciatore
}
\seealso{
  \code{\link{KODAMA.matrix.parallel}}
}
\examples{
  imgs <- photo(vis, xy, range = 0.1, n_pixels = 10)
}
\keyword{methods}
