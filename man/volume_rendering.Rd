\name{volume_rendering}
\alias{volume_rendering}
\title{3D Volume Rendering of Tissue Segments}
\description{
  Generates isosurface contours for tissue segments using kernel density estimation and renders them with rgl.
}
\usage{
  volume_rendering(xyz, tissue_segments, selection = NULL, alpha = NULL, colors = NULL, cells = c(20,20,20), level = exp(-3))
}
\arguments{
  \item{xyz}{
    3D coordinates of observations.
  }
  \item{tissue_segments}{
    Factor vector of segment labels.
  }
  \item{selection}{
    Optional list of segments to render.
  }
  \item{alpha}{
    Transparency values for each segment.
  }
  \item{colors}{
    Colors for each segment.
  }
  \item{cells}{
    Density grid resolution (default c(20,20,20)).
  }
  \item{level}{
    Isosurface threshold (default exp(-3)).
  }
}
\value{
  An rgl WebGL widget showing rendered surfaces.
}
\author{
  Stefano Cacciatore
}
\seealso{
  \code{\link[MASS]{kde3d}}, \code{\link[misc3d]{contour3d}}, \code{\link[rgl]{rglwidget}}
}
\examples{
  volume_rendering(xyz, seg, colors = c('red','blue'), alpha = c(0.5, 0.3))
}
\keyword{methods}
