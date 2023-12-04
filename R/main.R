RunKODAMAmatrix <- function(...) {
    UseMethod("RunKODAMAmatrix")
}

RunKODAMAvisualization <- function(...) {
    UseMethod("RunKODAMAvisualization")
}


refinecluster <- function(...) {
    UseMethod("refinecluster")
}


RunKODAMAmatrix.default = function(data, ...) {
  kk=KODAMA.matrix.parallel(data = data, ...)
  return(kk)
}

RunKODAMAmatrix.SingleCellExperiment = function(object, reduction= "PCA", dims=50, ...) {
  if (!is(object, "SingleCellExperiment")) {
    stop("object is not a SingleCellExperiment object")
  }
  if(reduction =="PCA"){
    data <- SingleCellExperiment::reducedDim(object, "PCA")
    nc <- ncol(data)
    if(nc < dims){
      dims = nc
      message("dims is set higher than the number of principal components")
    }
    data= data[ , 1:nc]
  }
  kk <- KODAMA.matrix.parallel(data = data, ...)
  object@int_colData@listData[["reducedDims"]]@listData[["KODAMA"]] <- kk
  return(object)
}

RunKODAMAmatrix.SpatialExperiment = function(object, reduction= "PCA", dims=50, ...) {
  if (!is(object, "SpatialExperiment")) {
    stop("object is not a spatialExperiment object")
  }
  if(reduction =="PCA"){
    data <- SingleCellExperiment::reducedDim(object, "PCA")
    nc <- ncol(data)
    if(nc < dims){
      dims = nc
      message("dims is set higher than the number of principal components")
    }
    
    data= data[ , 1:nc]
  }
  #spat_coord = NULL
  # In this dataset, the names of the assays are "counts" and "logcounts"
  #if("Spatial"  %in%  names(object@assays)){
    spat_coord<- as.matrix(SpatialExperiment::spatialCoords(object))
   # }

  kk=KODAMA.matrix.parallel(data = data, spatial = spat_coord, ...)
  #I have this error when I run the following line
  
  #SingleCellExperiment::reducedDim(object, "KODAMA") <- kk
  #Error in .set_internal_character(x, type, value, getfun = int_colData,  : 
  #invalid 'value' in 'reducedDim(<SpatialExperiment>, type="character") <- value':
  # 'value' should have number of rows equal to 'ncol(x)'
  
  # X FARAG: a class is organize with specific variables. You cannot simply assign a list to the class reducedDim.
    
 
    
  
  #So I assigned KODAMA object manually
  object@int_colData@listData[["reducedDims"]]@listData[["KODAMA"]] <- kk
  return(object)
}



RunKODAMAmatrix.giotto = function(object,reduction="pca",dims=50, ...) {
  if (!is(object, "giotto")) {
    stop("object is not a giotto object")
  }
  if(reduction=="pca"){
    data=Giotto::getDimReduction(object,reduction = c("cells", "feats"),reduction_method = "pca",name = "pca",output =  "matrix",set_defaults = TRUE)
    nc=ncol(data)
    if(nc<dims){
      dims=nc
      message("dims is set higher than the number of principal components")
    }
    data=data[,1:nc]
  }
    #expression_data= Seurat::GetAssayData(brain, assay = assay)
  spat_coord = NULL

   # Giotto pipeline deals only with spatial omics data
   spat_coord=getSpatialLocations(object,spat_unit = NULL,name = NULL,output = "data.table",copy_obj = TRUE,verbose = TRUE,set_defaults = TRUE)
   xy_names=spat_coord$cell_ID 
   spat_coord=as.matrix(spat_coord[,-3])
   rownames(spat_coord)=xy_names
   
      
  kk=KODAMA.matrix.parallel(data = data, spatial = spat_coord, ...)

dimObject=createDimObj(
  coordinates=matrix(0),
  name = "KODAMA",
  spat_unit = "cell",
  feat_type = "rna",
  method = "KODAMA",
  reduction = "cells",
  provenance = NULL,
  misc = kk,
  my_rownames = NULL)

  object = set_dimReduction(gobject = object, dimObject = dimObject)

  return(object)
}

RunKODAMAmatrix.Seurat = function(object,reduction="pca",dims=50,assay="counts",  ...) {
  if (!is(object, "Seurat")) {
    stop("object is not a Seurat object")
  }
  if(reduction=="pca"){
    data=Seurat::Embeddings(object, reduction = "pca")
    nc=ncol(data)
    if(nc<dims){
      dims=nc
      message("dims is set higher than the number of principal components")
    }
    data=data[,1:nc]
  }
  if(is.null(reduction)){
      data= Seurat::GetAssayData(object, assay = assay)
  }  
    
  spat_coord = NULL
  if(names(object@assays)=="Spatial")
      spat_coord<- as.matrix(Seurat::GetTissueCoordinates(object))
      
  kk=KODAMA.matrix.parallel(data = data, spatial = spat_coord, ...)
  KODAMA=CreateDimReducObject(
    embeddings = data[,1:2],
    key = "Dimensions_",
    assay = "RNA",
    misc=kk
  )
  object@reductions$KODAMA=KODAMA
  return(object)
}

RunKODAMAvisualization.default = function(kk, ...) {

  vis <- KODAMA.visualization(kk, ...)
  return(vis)
}

RunKODAMAvisualization.SingleCellExperiment <- function(object, ...) {
  if (!is(object, "SingleCellExperiment")) {
    stop("object is not a SingleCellExperiment object")
  }
  reducedDims_KODAMA <- object@int_colData@listData[["reducedDims"]]@listData[["KODAMA"]]
  vis <- KODAMA.visualization(reducedDims_KODAMA, ...)
  object@int_colData@listData[["reducedDims"]]@listData[["KODAMA"]] <- vis
  return(object)
}

RunKODAMAvisualization.SpatialExperiment = function(object, ...) {
  if (!is(spe, "SpatialExperiment")) {
    stop("object is not a SpatialExperiment object")
  }
  reducedDims_KODAMA <- object@int_colData@listData[["reducedDims"]]@listData[["KODAMA"]]
  vis <- KODAMA.visualization(reducedDims_KODAMA, ...)
  SingleCellExperiment::reducedDim(object, "KODAMA") <- vis
  return(object)
}

RunKODAMAvisualization.Seurat = function(object, ...) {
  if (!is(object, "Seurat")) {
    stop("object is not a Seurat object")
  }
  vis=KODAMA.visualization(object@reductions$KODAMA@misc, ...)
  KODAMA=CreateDimReducObject(
    embeddings = vis,
    key = "Dimensions_",
    assay = "RNA",
    misc=object@reductions$KODAMA@misc
  )
  object@reductions$KODAMA=KODAMA
  return(object)
}


RunKODAMAvisualization.giotto = function(object, ...) {
  if (!is(object, "giotto")) {
    stop("object is not a giotto object")
  }
    
  vis=KODAMA.visualization(object@dimension_reduction$cells$cell$rna$KODAMA$KODAMA@misc, ...)

dimObject=createDimObj(
  coordinates=vis,
  name = "KODAMA",
  spat_unit = "cell",
  feat_type = "rna",
  method = "KODAMA",
  reduction = "cells",
  provenance = NULL,
  misc = object@dimension_reduction$cells$cell$rna$KODAMA$KODAMA@misc,
  my_rownames = NULL)

  object = set_dimReduction(gobject = object, dimObject = dimObject,verbose = FALSE)



    
  return(object)
}


refinecluster.default = function (clusterlabels, location, shape = shape){
  t <- refine_cluster(clusterlabels, location, shape = shape) 
  return(t)
}

refinecluster.SpatialExperiment = function (object, shape = "square", assay = "Spatial"){
  if (!is(object, "SpatialExperiment")) {
    stop("object is not a SpatialExperiment object")
  }
  #   clusterlabels <- SingleCellExperiment::clusterCells(object, use.dimred="KODAMA", assay.type="logcounts")
  # should we add assay.type to be spatial?
  clusterlabels <- scran::clusterCells(object, use.dimred="KODAMA")
  location <- as.matrix(SpatialExperiment::spatialCoords(object))
  t <- refine_cluster(clusterlabels, location, shape = shape) 
  colLabels(object) <- t
  return(object)
}


refinecluster.Seurat = function (object, shape = "square"){
  if (!is(object, "Seurat")) {
    stop("object is not a Seurat object")
  }
  t=refine_cluster(object@active.ident, as.matrix(Seurat::GetTissueCoordinates(object)), shape = shape) 
  object@active.ident=as.factor(t)
  return(object)
}


refinecluster.giotto = function (object,name, shape = "square"){
  if (!is(object, "giotto")) {
    stop("object is not a Seurat object")
  }
  spat_coord=getSpatialLocations(object,spat_unit = NULL,name = NULL,output = "data.table",copy_obj = TRUE,verbose = TRUE,set_defaults = TRUE)
   xy_names=spat_coord$cell_ID 
   spat_coord=as.matrix(spat_coord[,-3])
   rownames(spat_coord)=xy_names
   cluster=object@cell_metadata$cell$rna@metaDT[,..name][[name]]
  t=refine_cluster(cluster, spat_coord, shape = shape) 
  t=data.frame(t)
  colnames(t)="refined"
  object = addCellMetadata(gobject = object,new_metadata = t)

  return(object)
}
    

 refine_SVM = function(xy,labels,samples,cost =  c(0.001,0.01,0.1,1,10,100), gamma=10^(-2:0), ...){
   samples=as.factor(samples)
   labels=as.factor(labels)
   sa=levels(samples)
   refine=rep(NA,length(labels))
   for(s in sa){
     print(s)
     sel=samples==s
     xr = xy[sel,]
     yr <- as.factor(as.vector(labels[sel]))
     tuned <- tune.svm(x=xr, y=yr, 
                       cost =  cost,
                       gamma = gamma)
     
     model <-svm(x=xr, y=yr,
                 cost = tuned$best.parameters$cost, 
                 gamma = tuned$best.parameters$gamma,
     )
     refine[sel] <- as.vector(fitted(model))
   }
   refine=factor(refine,levels=levels(labels))
   refine
 }
 

new_trajectory = function(x,y,n=20,data=NULL,knn=10,FUN=mean){
  ii=identify(x,y,order = TRUE)
  ii=ii$ind[order(ii$order)]
  dd= xspline(x[ii], y[ii], shape = c(0,rep(-1, 10-2),0), border="red",draw = FALSE)
  ll=length(dd$x)
  sel=seq(1,ll,length.out =n)
  
 # points(dd,col=2,bg="#eeeeee",lwd=2,pch=21)
  dd$x=dd$x[sel]
  dd$y=dd$y[sel]
  xy=cbind(dd$x,dd$y)
  xy_total=cbind(x,y)
  selection=knn_Armadillo(xy_total,xy,k = knn)$nn_index
  if(!is.null(data)){
     trajectory=apply(selection,1,function(z) apply(data[z,],2,FUN))
  }
  points(dd,col=2,bg="#eeeeee",lwd=2,pch=21)
  list(xy=dd,selection=selection,trajectory=trajectory,
       settings=list(x=x,y=y,n=n,data=data,knn=knn,FUN=FUN))
}

add_branch = function(dd){
  n_start=identify(dd$xy,n=1)
  start_x=dd$xy$x[n_start]
  start_y=dd$xy$y[n_start]
  ii=identify(x,y,order = TRUE)
  ii=ii$ind[order(ii$order)]
  
  branch= xspline(c(start_x,x[ii]), c(start_y,y[ii]), shape = c(0,rep(-1, 10-2),0), border="red",draw = FALSE)
  ll=length(branch$x)
  sel=seq(1,ll,length.out =dd$settings$n-n_start+1)
  
  # points(dd,col=2,bg="#eeeeee",lwd=2,pch=21)
  branch$x=branch$x[sel]
  branch$y=branch$y[sel]
  
  
  xy=cbind(branch$x,branch$y)
  xy_total=cbind(dd$settings$x,dd$settings$y)
  if(!is.null(data)){
    selection=knn_Armadillo(xy_total,xy,k = dd$settings$knn)$nn_index
    trajectory=apply(selection,1,function(z) apply(data[z,],2,FUN))
  }
  points(branch,col=3,bg="#eeeeee",lwd=2,pch=21)
  
}













                              
KODAMA.matrix.parallel =
function (data,                       # Dataset
          spatial = NULL,             # In spatial are conteined the spatial coordinates of each entries
          M = 100, Tcycle = 20, 
          FUN = c("PLS","PK","KNN"), 
          f.par.knn = 5, f.par.pls = 5,
          W = NULL, 
          constrain = NULL, fix = NULL, epsilon = 0.05, landmarks = 10000,  
          splitting = 50, spatial.resolution = 0.3, n.cores = 1) 
{
  neighbors = min(c(landmarks, nrow(data)/3)) + 1
  if (sum(is.na(data)) > 0) {
    stop("Missing values are present")
  } 
  data = as.matrix(data)
  nsample = nrow(data)
  nvariable = ncol(data)
  nsample_spatial= nrow(spatial)
  
  if (is.null(spatial)) {
    spatial_flag = FALSE
  }  else {
    spatial_flag = TRUE
  }
  if (is.null(fix)) 
    fix = rep(FALSE, nsample)
  if (is.null(constrain)) 
    constrain = 1:nsample
  shake = FALSE
  
  if(nsample<=landmarks){
    landmarks=ceiling(nsample*0.75)
    simm_dissimilarity_matrix=TRUE
  } else{
    simm_dissimilarity_matrix=FALSE
  }

  nspatialclusters=round(landmarks*spatial.resolution)
  
  QC=quality_control(data_row = nsample,
                     data_col = nvariable,
                     spatial_row = nsample_spatial,
                     FUN = FUN,
                     data = data,
                     f.par.knn = f.par.knn,
                     f.par.pls = f.par.pls)
  matchFUN=QC$matchFUN

  f.par.pls=QC$f.par.pls



  res = matrix(nrow = M, ncol = nsample)


  vect_acc = matrix(NA, nrow = M, ncol = Tcycle)
  accu = NULL

    

  my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "PSOCK"
  )
  print(my.cluster)
  doParallel::registerDoParallel(cl = my.cluster)
  foreach::getDoParRegistered()
  foreach::getDoParWorkers()



doSNOW::registerDoSNOW(my.cluster)
pb <- txtProgressBar(min = 1, max = M, style = 1)
  

res_parallel <- foreach(k = 1:M, 
                  .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dopar%
{

    library("KODAMA")
  

    

  
    # The landmarks samples are chosen in a way to cover all different profile
    # The data are divided in a number of clusters equal to the number of landmarks
    # A landpoint is chosen randomly from each cluster

    landpoints=NULL
    clust = as.numeric(kmeans(data, landmarks)$cluster)
    for (ii in 1:landmarks) {
      www = which(clust == ii)
      lwww=length(www)
      landpoints = c(landpoints,www[sample.int(lwww, 1, FALSE, NULL)])
    }

    # Variables are splitted in two where 
    # X variables are the variables for cross-validation accuracy maximizzation
    # T variables are the variable for the projections
    
    Tdata = data[-landpoints, , drop = FALSE]
    Xdata = data[landpoints, , drop = FALSE]
    Tfix = fix[-landpoints]
    Xfix = fix[landpoints]
    whF = which(!Xfix)
    whT = which(Xfix)
    Tconstrain = as.numeric(as.factor(constrain[-landpoints]))
    Xconstrain = as.numeric(as.factor(constrain[landpoints]))
    Xspatial = spatial[landpoints, , drop = FALSE]
    Tspatial = spatial[-landpoints, , drop = FALSE]

    
 

    if (spatial_flag) {
      spatialclusters=as.numeric(kmeans(Xspatial, nspatialclusters)$cluster)
      tab = apply(table(spatialclusters, Xconstrain), 2,which.max)
      Xconstrain = as.numeric(as.factor(tab[as.character(Xconstrain)]))  
    }
  
    if (landmarks<200) {
      XW = Xconstrain
    } else {
      clust = as.numeric(kmeans(Xdata, splitting)$cluster)
      tab = apply(table(clust, Xconstrain), 2, which.max)
      XW = as.numeric(as.factor(tab[as.character(Xconstrain)]))
    }

    if(!is.null(W)){
      SV_startingvector = W[landpoints][ssa]
      unw = unique(SV_startingvector)
      unw = unw[-which(is.na(unw))]
      ghg = is.na(SV_startingvector)
      SV_startingvector[ghg] = as.numeric(as.factor(SV_startingvector[ghg])) + length(unw)
      tab = apply(table(SV_startingvector,Xconstrain), 2,  which.max)
      XW = as.numeric(as.factor(tab[as.character(Xconstrain)]))
    }
   
    
    clbest = XW
    options(warn = -1)
    yatta = 0
    attr(yatta, "class") = "try-error"
    while (!is.null(attr(yatta, "class"))) {
      yatta = try(core_cpp(Xdata, Tdata, clbest, Tcycle, FUN, 
                           f.par.knn,f.par.pls,
                           Xconstrain, Xfix, shake, Xspatial, 
                           Tspatial), silent = FALSE)

    }
    options(warn = 0)
    res_k=rep(NA,nsample)
    if (is.list(yatta)) {
      clbest = as.vector(yatta$clbest)
      accu = yatta$accbest
      yatta$vect_acc = as.vector(yatta$vect_acc)
      yatta$vect_acc[yatta$vect_acc == -1] = NA
      vect_acc[k, ] = yatta$vect_acc
      
      yatta$vect_proj = as.vector(yatta$vect_proj)
      yatta$vect_proj[Tfix] = W[-landpoints][Tfix]

      res_k[landpoints] = clbest
      res_k[-landpoints] = yatta$vect_proj



      
    }


    list(res_k=res_k)

  }
  
#  parallel::stopCluster(cl = my.cluster)
  
  for(k in 1:M){
    res[k,] = res_parallel[[k]]$res_k
  }
  close(pb)

  

    
#  dissimilarity=NULL
#  ma=NULL
#  if(simm_dissimilarity_matrix){
#    ma = matrix(0, ncol = nsample, nrow = nsample)
#    for(k in 1:M){
#      uni = unique(res[k,])
#      nun = length(uni)
#      res_k=res[k,]
#      for (ii in 1:nun) 
#        ma[res[k,] == uni[ii], res_k ==  uni[ii]] = ma[res_k == uni[ii], res_k == uni[ii]] + 1
#    }
#    ma = ma/M
#    Edist = as.matrix(dist(data))
#    ma[ma < epsilon] = 0
#
#  #  Entropy calculation
#  #    y = ma
#  #    diag(y) = NA
#  #    yy = as.numeric(y)
#  #    yy = yy[!is.na(yy)]
#  #    yy = yy/sum(yy)
#  #    H = -sum(ifelse(yy > 0, yy * log(yy), 0))
#    
#    mam = (1/ma) * Edist
#    mam[is.na(mam)] <- .Machine$double.xmax
#    mam[is.infinite(mam) & mam > 0] <- .Machine$double.xmax
#    mam = floyd(mam)
#    mam[mam == .Machine$double.xmax] <- NA#
#    prox = Edist/mam
#    diag(prox) = 1
#    prox[is.na(prox)] = 0
#    maxvalue = max(mam, na.rm = TRUE)
#    mam[is.na(mam)] = maxvalue
#
#    dissimilarity = mam
#  }

  knn_Armadillo = knn_Armadillo(data, data, neighbors + 1)
  knn_Armadillo$distances = knn_Armadillo$distances[, -1]
  knn_Armadillo$nn_index = knn_Armadillo$nn_index[, -1]


# KFOLD
f=floor(nsample/n.cores)
f.1=f+1
mod=nsample%%n.cores
ll=list()
if(mod!=0){
   for(i in 1:mod){
     ll[[i]]=(1:f.1)+(i-1)*f.1
   }
 }
 for(i in (mod+1):n.cores){
   ll[[i]]=(1:f)+(i-1)*f+mod
 }

print("Calculation of dissimilarity matrix...")

 pb <- txtProgressBar(min = 1, max = n.cores, style = 1)
  

 res_parallel <- foreach(k = 1:n.cores, 
                  .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dopar%
{

  library("KODAMA")
knn_nn_index=knn_Armadillo$nn_index[ll[[k]],]
knn_distances=knn_Armadillo$distances[ll[[k]],]
                                       
#  for (i_tsne in 1:nrow(data)) {
  for (i_tsne in 1:length(ll[[k]])) {
    for (j_tsne in 1:neighbors) {
      h_tsne=ll[[k]][i_tsne]
      kod_tsne = mean(res[, h_tsne] == res[, knn_nn_index[i_tsne, j_tsne]], na.rm = TRUE)
      knn_distances[i_tsne, j_tsne] = knn_distances[i_tsne,  j_tsne]/kod_tsne
        
    }
    oo_tsne = order(knn_distances[i_tsne, ])
    knn_distances[i_tsne, ] = knn_distances[i_tsne, oo_tsne]
    knn_nn_index[i_tsne, ] = knn_nn_index[i_tsne, oo_tsne]
  }


   list(knn_distances=knn_distances,knn_nn_index=knn_nn_index)

  }
  
  parallel::stopCluster(cl = my.cluster)
  
  for(k in 1:n.cores){
      
    knn_Armadillo$nn_index[ll[[k]],]=res_parallel[[k]]$knn_nn_index
    knn_Armadillo$distances[ll[[k]],] =res_parallel[[k]]$knn_distances

  }
  close(pb)


                                       

    
  knn_Armadillo$neighbors = neighbors
  return(list(acc = accu, 
              v = vect_acc, res = res, 
              knn_Armadillo = knn_Armadillo, 
              data = data))

#  return(list(dissimilarity = dissimilarity, acc = accu, proximity = ma, 
#              v = vect_acc, res = res, 
#              knn_Armadillo = knn_Armadillo, 
#              data = data))
}


