\name{read_annotations}
\alias{read_annotations}
\title{Read Compact Annotation Labels}
\description{
Reads the two-column annotation CSV format used by KODAMAextra and returns one
label for each unique identifier.
}
\usage{
read_annotations(address)
}
\arguments{
\item{address}{Path to a two-column CSV annotation file.}
}
\value{
A named character vector containing one parsed annotation per identifier.
}
