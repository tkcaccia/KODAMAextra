\name{plot_slide}
\alias{plot_slide}
\title{Faceted Spatial Domain Plot by Slide}
\description{
  Generates faceted scatter plots of spatial coordinates labeled by domains per slide/sample.
}
\usage{
  plot_slide(xy, slide, labels, col = NULL, nrow = 1, scales = 'free', size.dot = 3, size.legend.text = 10, size.legend.title = 10, size.legend.dot = 5, size.strip.text = 10)
}
\arguments{
  \item{xy}{
    2D coordinates matrix.
  }
  \item{slide}{
    Vector of slide/sample identities.
  }
  \item{labels}{
    Domain or cluster labels per point.
  }
  \item{col}{
    Optional vector of colors.
  }
  \item{nrow}{
    Number of rows in facet grid.
  }
  \item{scales}{
    Scaling option for facets.
  }
  \item{size.dot}{
    Dot size.
  }
  \item{size.legend.text}{
    Legend text size.
  }
  \item{size.legend.title}{
    Legend title size.
  }
  \item{size.legend.dot}{
    Dot size in legend.
  }
  \item{size.strip.text}{
    Facet label text size.
  }
}
\value{
  A ggplot object.
}
\author{
  Stefano Cacciatore
}
\seealso{
  \code{\link[ggplot2]{facet_wrap}}
}
\examples{
  p <- plot_slide(coords, slides, domains); print(p)
}
\keyword{methods}
