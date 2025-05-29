\name{RunKODAMAmatrix}
\alias{RunKODAMAmatrix}
\alias{RunKODAMAmatrix,default-method}
\title{Generic KODAMA Feature Extraction for Spatial Omics Data}
\description{
  KODAMA (KnOwledge Discovery by Accuracy MAximization) is an unsupervised and  
  semi-supervised learning algorithm for feature extraction from noisy Spatial  
  Omics datasets. \code{RunKODAMAmatrix} is the generic S3 method that dispatches  
  to class-specific implementations based on the input spatial object.
}
\usage{
  RunKODAMAmatrix(object, reduction = NULL, dims = NULL, assay = NULL, \dots)
}
\arguments{
  \item{object}{A spatial omics data object, such as a \linkS4class{SingleCellExperiment},  
    Giotto object, Seurat object, or \linkS4class{SpatialExperiment}, containing  
    precomputed latent-space embeddings or assays.}
  \item{reduction}{Character. Name of the dimensionality reduction embedding to use  
    (e.g., \code{"PCA"}, \code{"UMAP"}, \code{"ICA"}). Defaults to \code{NULL},  
    in which case the primary assay is used.}
  \item{dims}{Integer vector. Which dimensions of the chosen reduction to retrieve  
    (e.g., \code{1:10}). Defaults to \code{NULL}, meaning all available dimensions.}
  \item{assay}{Character. Name of the assay to extract from \code{object} if  
    \code{reduction = NULL}. Defaults to \code{NULL}, which uses the default assay.}
  \item{\dots}{Additional arguments passed to \code{\link{KODAMA.matrix.parallel}}  
    and to class-specific \code{RunKODAMAmatrix} methods.}
}
\value{
  An object of the same class as \code{object}, with KODAMA-derived feature  
  coordinates or reconstructed embeddings added to its \code{reducedDims} slot  
  (or equivalent metadata slot).}
\author{
  Stefano Cacciatore
}
\seealso{
  \code{\link{KODAMA.matrix.parallel}},  
  \code{\link{RunKODAMAmatrix.SingleCellExperiment}},  
  \code{\link{RunKODAMAmatrix.Seurat}},  
  \code{\link{RunKODAMAmatrix.SpatialExperiment}}
}
\examples{
  ## Default dispatch on a Seurat object
  # sce <- RunKODAMAmatrix(seurat_obj, reduction = "PCA", dims = 1:20)
  
  ## On a SingleCellExperiment
  # sce <- RunKODAMAmatrix(sce, reduction = "UMAP", dims = 1:10)
  
  ## Passing additional KODAMA parameters
  # sce <- RunKODAMAmatrix(sce, reduction = "PCA", dims = 1:15,
  #                        nPerm = 100, cores = 4)
}
\keyword{methods}



