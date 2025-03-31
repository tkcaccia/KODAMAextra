RunKODAMAmatrix <- function(...) {
    UseMethod("RunKODAMAmatrix")
}

RunKODAMAvisualization <- function(...) {
    UseMethod("RunKODAMAvisualization")
}


refinecluster <- function(...) {
    UseMethod("refinecluster")
}

multi_SPARKX <- function(...) {
    UseMethod("multi_SPARKX")
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
    data= data[ , 1:dims]
  }
  kk <- KODAMA.matrix.parallel(data = data, ...)
#  object@int_colData@listData[["reducedDims"]]@listData[["KODAMA"]] <- kk
 mcols(object@int_colData@listData$reducedDims)$info = list(KODAMA=kk)
  return(object)
}

RunKODAMAmatrix.SpatialExperiment = function(object, reduction= "PCA", dims=50, ...) {
  if (!is(object, "SpatialExperiment")) {
    stop("object is not a spatialExperiment object")
  }
  ##############################################
    #CHECK
    data <- SingleCellExperiment::reducedDim(object, reduction)
    nc <- ncol(data)
    if(nc < dims){
      dims = nc
      message("dims is set higher than the number of features")
    }
    
    data= data[ , 1:dims]
  ##############################################
      
  #spat_coord = NULL
  # In this dataset, the names of the assays are "counts" and "logcounts"
 

  spat_coord <- as.matrix(SpatialExperiment::spatialCoords(object))
  samples <- colData(object)$sample_id
      
  kk=KODAMA.matrix.parallel(data = data, spatial = spat_coord, samples = samples, ...)

   
  
  #So I assigned KODAMA object manually
#  object@int_colData@listData[["reducedDims"]]@listData[["KODAMA"]] <- kk
  mcols(object@int_colData@listData$reducedDims)$info = list(KODAMA=kk)

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
    data=data[,1:dims]
  }
    #expression_data= Seurat::GetAssayData(brain, assay = assay)
  spat_coord = NULL

   # Giotto pipeline deals only with spatial omics data
   spat_coord=getSpatialLocations(object,spat_unit = NULL,name = NULL,output = "data.table",copy_obj = TRUE,verbose = TRUE,set_defaults = TRUE)
   xy_names=spat_coord$cell_ID 
   spat_coord=as.matrix(spat_coord[,-ncol(as.matrix(spat_coord))])
    
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

#' Perform KODAMA.matrix on a Seurat object.
#' 
#' @method RunKODAMAmatrix Seurat
#' @export
#' @param assay Name of assay to retrieve the data if dimension = null.
#' @rdname RunKODAMAmatrix 
RunKODAMAmatrix.Seurat <- function (object, reduction = "pca", dims = 50, ...) 
{

  if (is.list(object)){ 
    #----------------------------------------------------------------------------#
    #                Running KODAMA on a list of Seurat objects                  #
    #----------------------------------------------------------------------------#
    for(i in seq_along(object)){
      if (!is(object[[i]], "Seurat")) {
        stop("object is not a Seurat object")
      }
    }
    for (i in seq_along(object)){
      data <- Seurat::Embeddings(object[[i]], reduction = reduction)
      nc = ncol(data)
      if (nc < dims) {
        dims = nc
        message("dims is set higher than the number of dimensions")
      }
      data = data[, 1:dims]
      spat_coord <- GetTissueCoordinates(object[[i]])
      samples <- object[[i]]@meta.data$orig.ident
    
      kk <- KODAMAextra::KODAMA.matrix.parallel(data = data, 
                                                spatial = spat_coord, samples = samples, ...)
      KODAMA = CreateDimReducObject(embeddings = data[ , 1:2], 
                                    key = "Dimensions_", assay = "RNA", misc = kk)
      object[[i]]@reductions[["KODAMA"]] <- KODAMA
      
    }
  }else{
    #----------------------------------------------------------------------------#
    #                Extract from integrated  or merged Seurat object            #
    #----------------------------------------------------------------------------#
    n_slide=length(object@images)
    if (!is(object, "Seurat")) { # 1- extract data
      stop("object is not a Seurat object")
    }
    data <- Seurat::Embeddings(object, reduction = reduction)
    nc = ncol(data)
    if (nc < dims) {
      dims = nc
      message("dims is set higher than the number of dimensions")
    }
    data = data[, 1:dims]
    
    shift=c(0,0)
    spat_coord=NULL
    for (f in seq_along(object@images)){ # 2- extract coordinates
      new_slide=as.matrix(GetTissueCoordinates(object@images[[f]])[,c("x","y")])
      slide <- t(t(new_slide)+shift)
      shift=c(round(max(new_slide[,1])*1.2),0)
      spat_coord <- rbind(spat_coord, slide)
    }
######################################################################################
      
    samples <- object@meta.data$orig.ident
    

#######################################################################################

      
      
    kk = KODAMAextra::KODAMA.matrix.parallel(data = data, spatial = spat_coord, samples = samples, ...)
    KODAMA = CreateDimReducObject(embeddings = data[ , 1:2],  # should we choose larger number of dims
                                  key = "Dimensions_", assay = "RNA", misc = kk)
    object@reductions$KODAMA = KODAMA
  }
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
 # reducedDims_KODAMA <- object@int_colData@listData[["reducedDims"]]@listData[["KODAMA"]]
reducedDims_KODAMA = mcols(object@int_colData@listData$reducedDims)$info[["KODAMA"]]
    
  vis <- KODAMA.visualization(reducedDims_KODAMA, ...)
 # object@int_colData@listData[["reducedDims"]]@listData[["KODAMA"]] <- vis
    SingleCellExperiment::reducedDim(object, "KODAMA") <- vis
  return(object)
}

RunKODAMAvisualization.SpatialExperiment = function(object, ...) {
  if (!is(object, "SpatialExperiment")) {
    stop("object is not a SpatialExperiment object")
  }
  #reducedDims_KODAMA <- object@int_colData@listData[["reducedDims"]]@listData[["KODAMA"]]
  reducedDims_KODAMA = mcols(object@int_colData@listData$reducedDims)$info[["KODAMA"]]  
  vis <- KODAMA.visualization(reducedDims_KODAMA, ...)
  SingleCellExperiment::reducedDim(object, "KODAMA") <- vis
  return(object)
}

RunKODAMAvisualization.Seurat = function(object, ...) {

  if (is.list(object)){ 
    for(i in seq_along(object)){
      if (!is(object[[i]], "Seurat")) {
        stop("object is not a Seurat object")
      }
    }
    for (i in seq_along(object)){
      vis=KODAMA.visualization(object[[i]]@reductions$KODAMA@misc, ...)
      KODAMA=CreateDimReducObject(
        embeddings = vis,
        key = "Dimensions_",
        assay = "RNA",
        misc=object[[i]]@reductions$KODAMA@misc
      )
      object[[i]]@reductions$KODAMA=KODAMA
      
    }
  }else{
    n_slide=length(object@images)
    if (!is(object, "Seurat")) { # 1- extract data
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

  }
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
    

#refine_SVM =
#function (xy, labels, samples=NULL, ...) 
#{
#  if(is.null(samples)){
#    samples=rep(1,length(labels))
#  }
#  samples = as.factor(samples)
#  labels = as.factor(labels)
#  sa = levels(samples)
#  refine = rep(NA, length(labels))
#  for (s in sa) {
#    print(s)
#    sel = samples == s
#    xr = xy[sel, ]
#    yr <- as.factor(as.vector(labels[sel]))
#    if (length(levels(yr)) > 1) {
#      model <- svm(x = xr, y = yr, ...)
#      refine[sel] <- as.vector(fitted(model))
#    }
#    else {
#      refine[sel] <- yr
#    }
#  }
#  refine = factor(refine, levels = levels(labels))
#  refine
#}

  
refine_SVM =
  function (xy, labels, samples = NULL, newdata=NULL,newsamples=NULL,tiles=NULL, ...) 
  {
    nam=rownames(xy)
    names(labels)=nam
    dimensions=length(dim(xy))
    if(!is.null(tiles)){
      ex=list()
      for(k in 1:dimensions){
        ex=c(ex,list(1:tiles[k]))
      }
      grid_tiles=expand.grid(ex)
    }
    if (is.null(samples)) {
      samples = rep(1, length(labels))
    }
    samples = as.factor(samples)
    labels = as.factor(labels)
    sa = levels(samples)
    
    if(is.null(newdata)){
      refine = rep(NA, length(labels))
      names(refine) = rownames(xy)
    }else{
      refine = rep(NA, nrow(newdata))
      
      names(refine) = rownames(newdata)
    }
    
    for (s in sa) {
      print(s)
      sel = samples == s
      sel_newsamples = newsamples == s
      xr = xy[sel, ,drop=FALSE]
      yr <- as.factor(as.vector(labels[sel]))
      names(yr)=names(labels[sel])
      if(!is.null(tiles)){
        ll_tiles=list()
        for(ii in 1:dimensions){
          vect=xr[,ii]
          mi=min(vect)
          ma=max(vect)
          delta=(ma-mi)/tiles[ii]
          low= mi+(0:(tiles[ii]-1))*delta
          high= mi+(1:tiles[ii])*delta
          high[tiles[ii]]=high[tiles[ii]]+0.1
          ll_tiles[[ii]]=list()
          ll_tiles[[ii]]$low=low
          ll_tiles[[ii]]$high=high
          
        }       
        if(dimensions==2){
          for(jj in 1:nrow(grid_tiles)){
            xr_x=xr[,1]
            xr_y=xr[,2]
            sel_x=ll_tiles[[1]]$low[grid_tiles[jj,1]]<xr_x & xr_x<ll_tiles[[1]]$high[grid_tiles[jj,1]]
            sel_y=ll_tiles[[2]]$low[grid_tiles[jj,2]]<xr_y & xr_y<ll_tiles[[2]]$high[grid_tiles[jj,2]]
            xr_tiles=xr[sel_x & sel_y,]
            sel_name=rownames(xr_tiles)
            ### newdata per samplpes
            if(!is.null(newdata)){
              newdata_x=newdata[sel_newsamples,1]
              newdata_y=newdata[sel_newsamples,2]
              
              sel_newdata_x=ll_tiles[[1]]$low[grid_tiles[jj,1]]<newdata_x & newdata_x<ll_tiles[[1]]$high[grid_tiles[jj,1]]
              sel_newdata_y=ll_tiles[[2]]$low[grid_tiles[jj,2]]<newdata_y & newdata_y<ll_tiles[[2]]$high[grid_tiles[jj,2]]
              xr_newdata_tiles=newdata[sel_newsamples,][sel_newdata_x & sel_newdata_y,]
              sel_name_newdata=rownames(xr_newdata_tiles)
              
            }
            
            yr_tiles=yr[sel_name]
            
            if (length(unique(yr_tiles)) > 1) {
              model <- svm(x = xr_tiles, y = yr_tiles, ...)
              if(is.null(newdata)){
                
                refine[sel_name] <- as.vector(fitted(model))
              }else{
                refine[sel_name_newdata] <- as.vector(predict(model,newdata = xr_newdata_tiles))
              }
            }            else {
              if(is.null(newdata)){
                refine[sel_name] <- yr_tiles
              }else{
                refine[sel_name_newdata] <- unique(yr_tiles)
              }
              
            }
          }
          
        }
        
      }else{
        
        
        
        if (length(unique(yr)) > 1) {
          model <- svm(x = xr, y = yr, ...)
          if(is.null(newdata)){
            refine[sel] <- as.vector(fitted(model))
          }else{
            refine[sel_newsamples] <- as.vector(predict(model,newdata = newdata))
          }
        }
        else {
          refine[sel] <- yr
        }
      }
    }
    refine = factor(refine, levels = levels(labels))
    refine
    }
    




new_trajectory =
function (input, plotting = FALSE, n = 20, data = NULL, draw = TRUE, 
          knn = 10, trace = NULL, ...) 
{
  if (plotting) {
    plot(input, ...)
  }
  x = input[, 1]
  y = input[, 2]
  if (is.null(trace)) {
    ii = identify(x, y, order = TRUE, plot = FALSE)
    ii = ii$ind[order(ii$order)]
    xx = x[ii]
    yy = y[ii]
    dd = xspline(xx, yy, shape = c(0, rep(-1, 10 - 2), 0), 
                 draw = FALSE)
    ll = length(dd$x)
    new_x = NULL
    new_y = NULL
    for (i in 1:(ll - 1)) {
      xxi = dd$x[(0:1) + i]
      yyi = dd$y[(0:1) + i]
      d2 = data.frame(x = xxi, y = yyi)
      xp = seq(xxi[1], xxi[2], length.out = n)
      yp = predict(lm(y ~ x, data = d2), data.frame(x = xp))
      new_x = c(new_x, xp)
      new_y = c(new_y, yp)
    }
    di = NULL
    di[1] = 0
    for (i in 2:length(new_x)) {
      di[i] = di[i - 1] + dist(cbind(new_x[((-1):0) + i], 
                                     new_y[((-1):0) + i]))
    }
    section = seq(0, di[length(new_x)], length.out = n)
    temp = Rnanoflann::nn(as.matrix(di), as.matrix(section),1)$indices
    dd$x = new_x[temp]
    dd$y = new_y[temp]
  }  else {
    dd = trace
  }
  if (draw) {
    xspline(dd, lwd = 5, border = "white")
    xspline(dd, lwd = 3, border = "red")
  }
  xy = cbind(dd$x, dd$y)
  kk = Rnanoflann::nn(input, as.matrix(xy), knn)$indices
  if (!is.null(data)) {
    trajectory = t(apply(kk, 1, function(x) colMeans(as.matrix(data)[x, ])))
  }  else {
    trajectory = NULL
  }
  points(dd$x[1],dd$y[1], col = 2, bg = "#eeeeee", lwd = 2, pch = 22,cex=2)
  points(dd$x[-1],dd$y[-1], col = 2, bg = "#eeeeee", lwd = 2, pch = 21)
  list(xy = dd, trajectory = trajectory, kk = kk, settings = list(x = x, 
                                                                  y = y, n = n, data = data))
}




# plot(as.raster(imgData(spe)$data[[2]]))
# axis(1)
# axis(2)

# sel=colData(spe)$sample_id=="Acinar_Cell_Carcinoma"
# xy_sel=xy_original[sel,]
# xy_sel=xy_sel*0.02184042
# xy_sel[,2]=550-xy_sel[,2]
# points(xy_sel,cex=0.5,type="n")

# spe@int_metadata$imgData$data[[2]]
# data=t(logcounts(spe)[,sel])
# trace=nn$xy
# nn=new_trajectory (xy_sel,data = data,trace=trace)
# ma=multi_analysis(nn$trajectory,1:20,FUN="correlation.test",method="spearman")


                       


#plot(xy,col=rainbow(20)[dd$kk],pch=20)
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
    selection=Rnanoflann::nn(xy_total,xy,k = dd$settings$knn)$indices
    trajectory=apply(selection,1,function(z) apply(data[z,],2,FUN))
  }
  points(branch,col=3,bg="#eeeeee",lwd=2,pch=21)
  
}




KODAMA.matrix.parallel =
  function (data,                       # Dataset
            spatial = NULL,             # In spatial are conteined the spatial coordinates of each entries
            samples = NULL, 
            M = 100, Tcycle = 20, 
            FUN = c("fastpls","simpls"), 
            ncomp = 50,
            W = NULL, metrics="euclidean",
            constrain = NULL, fix = NULL, epsilon = 0.05, landmarks = 10000,  
            splitting = 50, spatial.resolution = 0.3, n.cores = 1, seed=1234) 
  {
    set.seed(seed)
    f.par.pls = ncomp
    neighbors = min(c(landmarks, nrow(data)-1),500) 
    if (sum(is.na(data)) > 0) {
      stop("Missing values are present")
    } 
    data = as.matrix(data)
    nsample = nrow(data)
    nvariable = ncol(data)
    nsample_spatial= nrow(spatial)


    writeLines("Calculating Network")
    
    knn_Rnanoflann = Rnanoflann::nn(data, data, neighbors + 1,method=metrics)
    knn_Rnanoflann$distances = knn_Rnanoflann$distances[, -1]
    knn_Rnanoflann$indices = knn_Rnanoflann$indices[, -1]  


      
        
    if (is.null(spatial)) {
      spatial_flag = FALSE
    }  else {
      spatial_flag = TRUE
      
        
      writeLines("\nCalculating Network spatial")
      knn_Rnanoflann_spatial = Rnanoflann::nn(spatial, spatial, neighbors,method="euclidean" )
      
      aa=colMeans(abs(spatial[knn_Rnanoflann_spatial$indices[,1],]-spatial[knn_Rnanoflann_spatial$indices[,20],]))*3
      
      # A horizontalization of the spatial information is done
      # Each different sample will be placed side by side
      
      if(!is.null(samples)){
        samples_names=names(table(samples))          
        if(length(samples_names)>1){
          ma=0
          for (j in 1:length(samples_names)) {
            sel <- samples_names[j] == samples
            spatial[sel, 1]=spatial[sel, 1]+ma
            ran=range(spatial[sel, 1])
            ma=ran[2]+ dist(ran)[1]*0.5
          }
        }
      }
      ##  aa=sqrt(apply(spatial,2,function(x) dist(range(x))))
      
    }
    if (is.null(fix)) 
      fix = rep(FALSE, nsample)
    if (is.null(constrain)) 
      constrain = 1:nsample
    is.na.constrain=is.na(constrain)
    if(any(is.na.constrain)){
      constrain=as.numeric(as.factor(constrain))
      constrain[is.na.constrain]=max(constrain,na.rm = TRUE)+(1:length(constrain[is.na.constrain]))
    }
 
    
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
                       f.par.pls = f.par.pls)
    matchFUN=QC$matchFUN
    
    f.par.pls=QC$f.par.pls
    
    
    
    res = matrix(nrow = M, ncol = nsample)
    res_constrain = matrix(nrow = M, ncol = nsample)
    
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
                            .options.snow = list(progress = function(n) setTxtProgressBar(pb, n)),
                            .packages = c('KODAMA','Rnanoflann')) %dopar%
      {
        
        set.seed(seed+k)
        
        
        
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
        Xspatial = spatial[landpoints, , drop = FALSE]
        
        
        if (spatial_flag) {
          delta=as.numeric(unlist(tapply(aa,1:length(aa),function(x) runif(nsample,-x,x))))
          spatialclusters=as.numeric(kmeans(spatial+delta, nspatialclusters)$cluster)
          ta_const=table(spatialclusters)
          ta_const=ta_const[ta_const>1]
          sel_cluster_1=spatialclusters %in% as.numeric(names(ta_const))
          if(sum(!sel_cluster_1)>0){
            spatialclusters[!sel_cluster_1]=spatialclusters[sel_cluster_1][Rnanoflann::nn(spatial[sel_cluster_1,],spatial[!sel_cluster_1,,drop=FALSE],1)$indices]
          }        
          
          
          #######      tab = apply(table(spatialclusters, constrain), 2,which.max)
          #######      constrain_clean = tab[as.character(constrain)]
          #######     Does not work for big number
          
          constrain_clean=NULL
          for(ic in 1:max(constrain)){
            sel_ic=ic==constrain
            constrain_clean[sel_ic]=as.numeric(names(which.max(table(spatialclusters[sel_ic]))))
          }
              
        }else{
          constrain_clean=constrain
          
        }
        Xconstrain = as.numeric(as.factor(constrain_clean[landpoints]))
        if(!is.null(W)){
          SV_startingvector = W[landpoints]
          unw = unique(SV_startingvector)
          unw = unw[-which(is.na(unw))]
          ghg = is.na(SV_startingvector)
          SV_startingvector[ghg] = as.numeric(as.factor(SV_startingvector[ghg])) + length(unw)
          
          #################  tab = apply(table(SV_startingvector,Xconstrain), 2,  which.max)
          ################   XW = as.numeric(as.factor(tab[as.character(Xconstrain)]))
          XW=NULL
          for(ic in 1:max(Xconstrain)){
            XW[ic==Xconstrain]=as.numeric(names(which.max(table(SV_startingvector[ic==Xconstrain]))))
          }
          
          
        }else{
          if (landmarks<200) {
            XW = Xconstrain
          } else {
            clust = as.numeric(kmeans(Xdata, splitting)$cluster)
            #############   tab = apply(table(clust, Xconstrain), 2, which.max)
            #############   XW = as.numeric(as.factor(tab[as.character(Xconstrain)]))
            
            XW=NULL
            for(ic in 1:max(Xconstrain)){
              XW[ic==Xconstrain]=as.numeric(names(which.max(table(clust[ic==Xconstrain]))))
            }
            
            
          }
        }
        
        
        clbest = XW
        options(warn = -1)
        yatta = 0
        attr(yatta, "class") = "try-error"
        while (!is.null(attr(yatta, "class"))) {
          yatta = try(core_cpp(Xdata, Tdata, clbest, Tcycle, FUN, 
                               f.par.pls,
                               Xconstrain, Xfix), silent = FALSE)
          
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
          
          if(!is.null(W))
            yatta$vect_proj[Tfix] = W[-landpoints][Tfix]
          
          temp=rep(NA,nsample)
          res_k[landpoints] = clbest
          res_k[-landpoints] = yatta$vect_proj
          
          
          ###########        tab = apply(table(res_k, constrain_clean), 2, which.max)
          ###########        res_k = as.numeric(as.factor(tab[as.character(constrain_clean)]))
          
          res_k_temp=NULL
          for(ic in 1:max(constrain_clean)){
            res_k_temp[ic==constrain_clean]=as.numeric(names(which.max(table(res_k[ic==constrain_clean]))))
          }
          res_k=res_k_temp 
          
          
        }
        
        
        
        list(res_k=res_k,constrain_k=constrain_clean)
      }
    
    writeLines("\nFinished parallel computation")
    
    close(pb)
    for(k in 1:M){
      res[k,] = res_parallel[[k]]$res_k
      res_constrain[k,]=res_parallel[[k]]$constrain_k
    }
    
    
    
    
    
    
    print("Calculation of dissimilarity matrix...")
    
    pb <- txtProgressBar(min = 1, max = nrow(data), style = 1)
    
    
    res_parallel <- foreach(k = 1:nrow(data), 
                            .options.snow = list(progress = function(n) setTxtProgressBar(pb, n)),
                            .packages = c('KODAMA')) %dopar%
      {
        
        
        
        
        knn_indices=knn_Rnanoflann$indices[k,]
        knn_distances=knn_Rnanoflann$distances[k,]
        
        
        
        mean_knn_distances=mean(knn_distances)                             
        for (j_tsne in 1:neighbors) {
          
          kod_tsne = mean(res[, k] == res[, knn_indices[j_tsne]], na.rm = TRUE)
          knn_distances[j_tsne] = (1+knn_distances[j_tsne])/(kod_tsne^2)
          
        }
        
        
        oo_tsne = order(knn_distances)
        knn_distances = knn_distances[oo_tsne]
        knn_indices = knn_indices[oo_tsne]
        
        list(knn_distances=knn_distances,knn_indices=knn_indices)
      }
    
    parallel::stopCluster(cl = my.cluster)
    
    
    for(k in 1:nrow(data)){
      
      knn_Rnanoflann$indices[k,]=res_parallel[[k]]$knn_indices
      knn_Rnanoflann$distances[k,] =res_parallel[[k]]$knn_distances
      
    }
    
    close(pb)
    
    
    knn_Rnanoflann$neighbors = neighbors
    return(list(acc = accu, 
                v = vect_acc, res = res, 
                knn_Rnanoflann = knn_Rnanoflann, 
                data = data,
                res_constrain=res_constrain))
    
  }





                     


                                                 
photo=function(vis,xy,range=0.05,n_pixels=25){
  
  image=matrix(0,nrow(xy),ncol=2*n_pixels^2)
  for(j in 1:nrow(xy)){
    print(j)
    tile_sel=xy[,1]>(xy[j,1]-range) & xy[,1]<(xy[j,1]+range) &
      xy[,2]>(xy[j,2]-range) & xy[,2]<(xy[j,2]+range) 
    tile_j=t(t(xy[tile_sel,,drop=FALSE])-xy[j,])+range
    
    tile_KODAMA_1=vis[tile_sel,1]
    tile_KODAMA_2=vis[tile_sel,2]
    
    tile_j=ceiling(tile_j*10*n_pixels)
    
    ma1=matrix(0,ncol=n_pixels,nrow=n_pixels)
    for(i in 1:nrow(tile_j)){
      ma1[tile_j[i,]]=tile_KODAMA_1[i]
    }
    ma2=matrix(0,ncol=n_pixels,nrow=n_pixels)
    for(i in 1:nrow(tile_j)){
      ma2[tile_j[i,]]=tile_KODAMA_2[i]
    }
    image[j,]=c(as.numeric(ma1),as.numeric(ma2))
    
  }
  image
}



volume_rendering <- function(xyz,  tissue_segments,selection=NULL, alpha=NULL, colors=NULL,cells=c(20, 20, 20), level=exp(-3)) {
  if(!is.factor(tissue_segments)){
    stop("tissue_segments is not a factor")
  }
  option_tissue=levels(tissue_segments)
  if(is.null(colors)){
    colors=rainbow(length(option_tissue))
  }else{
    if(length(option_tissue)!=length(alpha)){
      stop("The number of color does not match")
    }
  }
  if(is.null(alpha)){ 
    alpha=rep(1,length(option_tissue))
  }else{
    if(length(option_tissue)!=length(alpha)){
      stop("The number of alpha does not match")
    }
  }
  if(!is.null(selection)){
    option_tissue=selection
  }
    ww=which(levels(tissue_segments) %in% option_tissue)
    colors=colors[ww]
    alpha=alpha[ww]
    
    

  cells[1]=min(cells[1],length(unique(xyz[,1])))
  cells[2]=min(cells[2],length(unique(xyz[,2])))
  cells[3]=min(cells[3],length(unique(xyz[,3])))
  sel.alpha=alpha>0
  option_tissue=option_tissue[sel.alpha]
  alpha=alpha[sel.alpha]
  colors=colors[sel.alpha]
  for (i in 1:length(option_tissue) ){
    segment <- option_tissue[i]
    sel <- tissue_segments == segment
    d <- kde3d(xyz[sel, 1], xyz[sel, 2], xyz[sel, 3], n = cells)
    e=array(0,dim=cells+2)
    e[2:(cells[1]+1),2:(cells[2]+1),2:(cells[3]+1)]=d$d
    dx=c(d$x[1]-d$x[2],d$x,d$x[cells[1]]+d$x[2]-d$x[1])
    dy=c(d$y[1]-d$y[2],d$y,d$y[cells[2]]+d$y[2]-d$y[1])
    dz=c(d$z[1]-d$z[2],d$z,d$z[cells[3]]+d$z[2]-d$z[1])
    contour3d(e, level=level, dx, dy, dz,
              color = colors[i],  scale = FALSE,
              engine = "rgl", draw = TRUE, alpha = alpha[i], add = (i != 1))
  }
  rglwidget()
}


passing.message = 
  function (data, spatial, number_knn = 15) 
  {
    data=as.matrix(data)
    spatial=as.matrix(spatial)
    # Get dimensions of the input data
    nspots = nrow(data)
    nvariables = ncol(data)

    # Initialize result matrix
    data.2 <- matrix(0, nrow = nspots, ncol = nvariables)
    knn=Rnanoflann::nn(spatial,spatial,number_knn)
    for(h in 1:nspots){       
      # Initialize a temporary vector for computations
      temp = rep(0, nvariables)
      RNA.temp = data[knn$indices[h,], ]
      knn$distances=knn$distances/max(knn$distances)
      # Compute weighted sum for each variable
      for (i in 1:number_knn) {
        temp = temp + RNA.temp[i, ]/exp(knn$distances[h,i])
      }
      # Return the computed row
      data.2[h, ] <-temp
    }

    rownames(data.2) = rownames(data)
    colnames(data.2) = colnames(data)
    data.2
  }



  multi_SPARKX.default = function (data,spatial,samples=NULL,n.cores=1){
    nsamples=nrow(data)
    nfeatures=ncol(data)
    
    if(is.null(samples)){
      samples=as.factor(rep("",nrow(data)))
    }
    sample_names=levels(samples)
    
    pvalue_mat=matrix(nrow=nfeatures,ncol=length(sample_names))
    rownames(pvalue_mat)=colnames(data)
    
    for(i in 1:length(sample_names)){
      sel=samples==sample_names[i]
      data_sub= data[sel,]
      spatial_sub= spatial[sel,]
      options(warn=-1)
      sparkX <- sparkx(t(data_sub),spatial_sub,numCores=n.cores,option="mixture")
      options(warn=0)
      pvalue_mat[rownames(sparkX$res_mtest),i]=sparkX$res_mtest$combinedPval
      print(sample_names[i])
    }
    pvalue_mat=-log(pvalue_mat)
    zeros=apply(data,2,function(y) any(tapply(y,samples,function(x) all(x==0))))
    pvalue_mat[zeros,]=0
    #  oo=order(apply(pvalue_mat,1,function(x) median(-log(x),na.rm=TRUE)),decreasing = TRUE)
    oo=order(apply(pvalue_mat,1,function(x) max(x,na.rm=TRUE)),decreasing = TRUE)
    top=rownames(pvalue_mat)[oo]
    top
  }
  
  
  
  
  
  multi_SPARKX.SpatialExperiment = function(object, ...) {
    if (!is(object, "SpatialExperiment")) {
      stop("object is not a spatialExperiment object")
    }
    data.2=as.matrix(logcounts(object))
    data=t(data.2)
    spatial=as.matrix(spatialCoords(object))
    samples=as.factor(colData(object)$sample_id)
    top=multi_SPARKX(data,spatial,samples, ...)
    return(top)
  }
  





plot_slide = 
function (xy, slide, labels, col = NULL, nrow = 1,scales = "free") 
{
  if (is.null(col)) {
    labels = as.factor(labels)
    nn = length(levels(labels))
    col = rainbow(nn)
  }
  df <- data.frame(xy, slide, labels)
  colnames(df) = c("x", "y", "slide", "labels")
  df$slide = as.factor(df$slide)
  df$labels = as.factor(df$labels)
  ggplot(df, aes(x, y, color = labels)) + geom_point(size = 1) + 
    facet_wrap(~slide, nrow = nrow,scales = scales) + theme_bw() + theme(legend.position = "bottom", 
                                                         axis.title = element_blank(), axis.text = element_blank(), 
                                                         axis.ticks = element_blank(), panel.grid = element_blank()) + 
    scale_color_manual("Domain", values = col) + guides(color = guide_legend(nrow = 1, 
                                                                             override.aes = list(size = 2)))
}




                   
 

Lscore = function(data,l,knn=10){
  require(Rfast)
  score=NULL
  for(j in 1:7){
    uni=unique(l)
    sel=l==uni[j]
    data.sel=data[sel,]
    l.sel=l[sel]
    nr=nrow(data.sel)
    if(nr>knn){
    dod=matrix(NA,ncol=nr,nrow=nr)
    di=Rnanoflann::nn(data.sel,data.sel,1+knn)
    for(i in 1:nr){
      dod[i,di$indices[i,]]=di$distances[i,]
    }
    ff=Rfast::floyd(dod)
     score[j]=max(ff,na.rm = TRUE)/sqrt((dist(range(data.sel[,1]))*dist(range(data.sel[,2]))))
  #  score[j]=max(ff,na.rm = TRUE)/nr
    }
  }
  score

}
                            

















#     library(ggplot2)
#     library(sf)
#     library(dplyr)

# Sample data for scatter plot
#     set.seed(123)
#     x <- rnorm(500)*100+6000
#     y <- rnorm(500)*100+6000

# Create a data frame from the scatter plot data
#     data <- data.frame(x = x, y = y)

# Create a base scatter plot with contours using ggplot2
#     ggplot(data, aes(x = x, y = y)) +
#       geom_point() +
#       stat_density_2d(aes(fill = ..level..), geom = "polygon") +
#       geom_density_2d()

#     density_est <- MASS::kde2d(x, y, n = 100)

#     contours <- contourLines(density_est$x, density_est$y, density_est$z)
#     contours=list(contours[[1]])

# Convert the contours to sf LINESTRING objects
#     sf_contours <- lapply(contours, function(contour) {
#       st_linestring(cbind(contour$x, contour$y))
#     })

# Combine all linestrings into an sf object
#     sf_contours <- st_sfc(sf_contours)

#     layer_name <- "Necrosis"
#     layer_color <- list(c(50, 50, 50))  # RGB color

# Create a simple feature (sf) object with attributes
#     sf_contours <- st_sf(
#       level = sapply(contours, `[[`, "level"),
#       geometry = sf_contours,
#       classification = "{\"name\":\"Necrosis\",\"color\":[50, 50, 50]}"    
#     )

# Optionally, add some attribute data (like the contour level)
#     sf_contours <- st_sf(level = sapply(contours, `[[`, "level"), geometry = sf_contours)

# Save the contour lines as a GeoJSON file
#     st_write(sf_contours, "contour_plot4.geojson", driver = "GeoJSON",delete_dsn=TRUE) #append=TRUE





read_annotations = function(address){
  rr <- read.csv(address)
ss=strsplit(rr[,2],":")
ss=unlist(lapply(ss, function(x) x[2]))
ss=strsplit(ss,",")
ss=unlist(lapply(ss, function(x) x[1]))
ss=gsub("\"","",ss)
names(ss)=rr[,1]
n=ave(1:length(rr[,1]), rr[,1], FUN = seq_along)
ss=ss[n==1]
# Remove the first character
ss <- substr(ss, 2, nchar(ss))
ss
}






leiden = function(g,ncluster,init=0,delta=0.2){
  ##### Leiden algorithm
  res=init
  t=0
  res=init-delta
  while(ncluster>t){
    res=res+delta
    clu=cluster_leiden(g,resolution_parameter=res)
    t=clu$nb_clusters
    print(c(res,t))
  }

  if(t!=ncluster){
    res=res-delta
    t=0
  }
  while(ncluster!=t){
    delta=delta/2
    if(t>ncluster){
      res=res-delta
    }else{
      res=res+delta
    }
    clu=cluster_leiden(g,resolution_parameter=res)
    t=clu$nb_clusters
    print(c(res,t))
  }
  clu
}

louvain = function(g,ncluster,init=0,delta=0.2){
  ##### Louvain algorithm
  res=init
  t=0
  res=init-delta
  while(ncluster>t){
    res=res+delta
    clu=cluster_louvain(g,resolution=res)
    t=length(unique(clu$membership))
    print(c(res,t))
  }

  if(t!=ncluster){
    res=res-delta
    t=0
  }
  while(ncluster!=t){
    delta=delta/2
    if(t>ncluster){
      res=res-delta
    }else{
      res=res+delta
    }
    clu=cluster_louvain(g,resolution=res)
    t=length(unique(clu$membership))
    print(c(res,t))
  }
  clu
}




                   

