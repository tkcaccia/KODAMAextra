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
    data= data[ , 1:dims]
  }
  kk <- KODAMA.matrix.parallel(data = data, ...)
  object@int_colData@listData[["reducedDims"]]@listData[["KODAMA"]] <- kk
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
  samples <- colData(spe)$sample_id
      
  kk=KODAMA.matrix.parallel(data = data, spatial = spat_coord, samples = samples, ...)

   
  
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
    

refine_SVM =
function (xy, labels, samples=NULL, ...) 
{
  if(is.null(samples)){
    samples=rep(1,length(labels))
  }
  samples = as.factor(samples)
  labels = as.factor(labels)
  sa = levels(samples)
  refine = rep(NA, length(labels))
  for (s in sa) {
    print(s)
    sel = samples == s
    xr = xy[sel, ]
    yr <- as.factor(as.vector(labels[sel]))
    if (length(levels(yr)) > 1) {
      model <- svm(x = xr, y = yr, ...)
      refine[sel] <- as.vector(fitted(model))
    }
    else {
      refine[sel] <- yr
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
    temp = knn_Armadillo(as.matrix(di), as.matrix(section), 
                         1)$nn_index
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
  kk = knn_Armadillo(input, as.matrix(xy), knn)$nn_index
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
    selection=knn_Armadillo(xy_total,xy,k = dd$settings$knn)$nn_index
    trajectory=apply(selection,1,function(z) apply(data[z,],2,FUN))
  }
  points(branch,col=3,bg="#eeeeee",lwd=2,pch=21)
  
}



KODAMA.matrix.parallel =
  function (data,                       # Dataset
            spatial = NULL,             # In spatial are conteined the spatial coordinates of each entries
            samples = NULL, 
            M = 100, Tcycle = 20, 
            FUN = c("PLS","PK","KNN"), 
            f.par.knn = 5, f.par.pls = 5,
            W = NULL, 
            constrain = NULL, fix = NULL, epsilon = 0.05, landmarks = 10000,  
            splitting = 50, spatial.resolution = 0.3, n.cores = 1, seed=1234) 
  {
    set.seed(seed)
    neighbors = min(c(landmarks, nrow(data)-1),1000) 
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
      snr=sqrt(nsample)
      aa=apply(spatial,2,function(x) dist(range(x))/snr)

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
    
    
    
#    res = matrix(nrow = M, ncol = nsample)
#    res_constrain = matrix(nrow = M, ncol = nsample)
    
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
    
    
    res <- big.matrix(M, nsample, type = "double", backingfile = "big_matrix.bin.res", descriptorfile = "big_matrix.desc.res")
    res_constrain <- big.matrix(M, nsample, type = "double", backingfile = "big_matrix.bin.res_constrain", descriptorfile = "big_matrix.desc.res_constrain")
    
    
    
    res_parallel <- foreach(k = 1:M, 
                            .options.snow = list(progress = function(n) setTxtProgressBar(pb, n)),
                            .packages = c('bigmemory','KODAMA')) %dopar%
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
          #    clu=sample(nsample,nspatialclusters)
          #    spatialclusters=knn_Armadillo(spatial[clu,],spatial,1)$nn_index
          #    tab = apply(table(spatialclusters, constrain), 2,which.max)
          #    constrain_clean = tab[as.character(constrain)]
          delta=as.numeric(unlist(tapply(aa,1:length(aa),function(x) runif(nsample,-x,x))))
          
          
          spatialclusters=as.numeric(kmeans(spatial+delta, nspatialclusters)$cluster)
          ta_const=table(spatialclusters)
          ta_const=ta_const[ta_const>1]
          sel_cluster_1=spatialclusters %in% as.numeric(names(ta_const))
          if(sum(!sel_cluster_1)>0){
            spatialclusters[!sel_cluster_1]=spatialclusters[sel_cluster_1][knn_Armadillo(spatial[sel_cluster_1,],spatial[!sel_cluster_1,,drop=FALSE],1)$nn_index]
          }        

                                         
    #######      tab = apply(table(spatialclusters, constrain), 2,which.max)
    #######      constrain_clean = tab[as.character(constrain)]
    ####### Does not work for big number

            constrain_clean=NULL
         for(ic in 1:max(constrain)){
           constrain_clean[ic==constrain]=as.numeric(names(which.max(table(spatialclusters[ic==constrain]))))
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
                               f.par.knn,f.par.pls,
                               Xconstrain, Xfix, shake), silent = FALSE)
          
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
        
        
        res <- attach.big.matrix("big_matrix.desc.res")  # Attach the big.matrix in the worker
        res_constrain <- attach.big.matrix("big_matrix.desc.res_constrain")  # Attach the big.matrix in the worker <- attach.big.matrix(desc_path)  # Attach the big.matrix in the worker
        
        res[k, ] <- res_k  # Store the result in the big.matrix
        res_constrain[k, ] <- constrain_clean  # Store the result in the big.matrix
        
        NULL  # Return NULL to avoid accumulating results in memory
      }
    
    writeLines("\nFinished parallel computation")

    close(pb)
    
    
    ######## knn_Armadillo = knn_Armadillo(data, data, neighbors + 1)
    ######## knn_Armadillo$distances = knn_Armadillo$distances[, -1]
    ######## knn_Armadillo$nn_index = knn_Armadillo$nn_index[, -1]
    
    
    KKnn_index <- big.matrix(nsample,neighbors, type = "double", backingfile = "big_matrix.bin.KKnn_index", descriptorfile = "big_matrix.desc.KKnn_index")
    KKdistances <- big.matrix(nsample, neighbors, type = "double", backingfile = "big_matrix.bin.KKdistances", descriptorfile = "big_matrix.desc.KKdistances")
    
    
    print("Calculation of dissimilarity matrix...")
    
    pb <- txtProgressBar(min = 1, max = nrow(data), style = 1)
    
    
    res_parallel <- foreach(k = 1:nrow(data), 
                            .options.snow = list(progress = function(n) setTxtProgressBar(pb, n)),
                            .packages = c('bigmemory','KODAMA')) %dopar%
      {
        
        res <- attach.big.matrix("big_matrix.desc.res")  # Attach the big.matrix in the worker
        res_constrain <- attach.big.matrix("big_matrix.desc.res_constrain")  # Attach the big.matrix in the worker <- attach.big.matrix(desc_path)  # Attach the big.matrix in the worker

        KKnn_index <- attach.big.matrix("big_matrix.desc.KKnn_index")  # Attach the big.matrix in the worker
        KKdistances <- attach.big.matrix("big_matrix.desc.KKdistances")  # Attach the big.matrix in the worker <- attach.big.matrix(desc_path)  # Attach the big.matrix in the worker
        
        
        
    
        ########### knn_nn_index=knn_Armadillo$nn_index[k,]
        ########### knn_distances=knn_Armadillo$distances[k,]
        
        knn_Armadillo = knn_Armadillo(data,data[k,,drop=FALSE], neighbors + 1)
        knn_distances = knn_Armadillo$distances[, -1]
        knn_nn_index = knn_Armadillo$nn_index[, -1]  
          
        mean_knn_distances=mean(knn_distances)                             
        for (j_tsne in 1:neighbors) {
          
          kod_tsne = mean(res[, k] == res[, knn_nn_index[j_tsne]], na.rm = TRUE)
          knn_distances[j_tsne] = (1+knn_distances[j_tsne])/(kod_tsne^2)
          
        }
        
        
        oo_tsne = order(knn_distances)
        knn_distances = knn_distances[oo_tsne]
        knn_nn_index = knn_nn_index[oo_tsne]
        
        KKnn_index[k, ] = knn_nn_index
        KKdistances[k, ] = knn_distances
          
    #    list(knn_distances=knn_distances,knn_nn_index=knn_nn_index)
     NULL   
      }
    
    parallel::stopCluster(cl = my.cluster)
    
    knn_Armadillo=list()
    knn_Armadillo$nn_index=KKnn_index
    knn_Armadillo$distances=KKdistances
                                                                
                            
   # for(k in 1:nrow(data)){
      
   #   knn_Armadillo$nn_index[k,]=res_parallel[[k]]$knn_nn_index
   #   knn_Armadillo$distances[k,] =res_parallel[[k]]$knn_distances
      
  #  }
    close(pb)
    
    
    
    
    
    knn_Armadillo$neighbors = neighbors
    return(list(acc = accu, 
                v = vect_acc, res = res, 
                knn_Armadillo = knn_Armadillo, 
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




passing.message <- function(data, spatial, number_knn = 10, quantile = 0.5, n.cores = 1) {
  
  # Get dimensions of the input data
  nspots <- nrow(data)
  nvariables <- ncol(data)
  
  # Set up parallel backend
  my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
  doParallel::registerDoParallel(cl = my.cluster)
  pb <- txtProgressBar(min = 1, max = nspots, style = 1)
  
  # Initialize result matrix
  data.2 <- matrix(0, nrow = nspots, ncol = nvariables)
  
  # Parallel processing
  res_parallel <- foreach(h = 1:nspots, .options.snow = list(progress = function(n) setTxtProgressBar(pb, n)),
                          .packages = c("KODAMA")) %dopar% {
                            
                            # Find nearest neighbors using spatial information
                            knn <- knn_Armadillo(spatial, spatial[h, , drop = FALSE], number_knn)
                            
                            # Initialize a temporary vector for computations
                            temp <- numeric(nvariables)
                            RNA.temp <- data[knn$nn_index, ]
                            
                            # Compute knn_gene (optional if needed for further calculations)
                            knn_gene <- knn_Armadillo(RNA.temp, RNA.temp[1, , drop = FALSE], 
                                                      round(number_knn * quantile))$nn_index
                            
                            # Compute weighted sum for each variable
                            for (i in 1:number_knn) {
                              temp <- temp + RNA.temp[i, ] / (1 + knn$distances[, i]^2)
                            }
                            
                            # Return the computed row
                            temp
                          }
  
  # Combine results back into data.2 matrix
  for (h in 1:nspots) {
    data.2[h, ] <- res_parallel[[h]]
  }
  
  # Stop the cluster
  parallel::stopCluster(cl = my.cluster)
  
  # Set row and column names to match the input data
  rownames(data.2) <- rownames(data)
  colnames(data.2) <- colnames(data)
  
  # Return the result matrix
  data.2
}




#spe_sub=spe[,colData(spe)$sample_id=="151507"]
#data=as.matrix(t(assay(spe_sub,"counts")[1:2000,]))
#xy=spatialCoords(spe_sub)
#spe=passing.message(data,xy,n.cores=4)

                            
                                                 



