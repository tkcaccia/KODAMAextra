\name{new_trajectory}
\alias{new_trajectory}
\title{Interactively Define a Two-Dimensional Trajectory}
\description{
Defines a smooth trajectory through a two-dimensional embedding and optionally
summarizes data over nearest observations along that trajectory.
}
\usage{
new_trajectory(
  input, plotting = FALSE, n = 20, data = NULL, draw = TRUE,
  knn = 10, trace = NULL, ...
)
}
\arguments{
\item{input}{Numeric matrix with two embedding coordinates per row.}
\item{plotting}{Whether to plot \code{input} before interactive selection.}
\item{n}{Number of points retained along the trajectory.}
\item{data}{Optional observation-by-variable matrix to summarize.}
\item{draw}{Whether to draw the fitted spline.}
\item{knn}{Number of neighboring observations summarized per trajectory point.}
\item{trace}{Optional precomputed list with numeric \code{x} and \code{y}.}
\item{...}{Additional arguments passed to \code{plot}.}
}
\value{
A list containing the fitted coordinates, optional summarized trajectory,
neighbor indices, and settings.
}
\examples{
\dontrun{
set.seed(1)
embedding <- cbind(seq(0, 1, length.out = 100), stats::rnorm(100, sd = 0.05))
new_trajectory(embedding, plotting = TRUE)
}
}
