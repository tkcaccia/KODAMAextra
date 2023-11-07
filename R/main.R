
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

#  pb <- txtProgressBar(min = 1, max = M, style = 1)
#  for (k in 1:M) {
#    setTxtProgressBar(pb, k)

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

  
#  res_parallel=foreach(k=1:M) %dopar%   {
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
 # }
#  close(pb)

    list(res_k=res_k)

  }
  
  parallel::stopCluster(cl = my.cluster)
  
  for(k in 1:M){
    res[k,] = res_parallel[[k]]$res_k
  }
  close(pb)

  

    
  dissimilarity=NULL
  ma=NULL
  if(simm_dissimilarity_matrix){
    ma = matrix(0, ncol = nsample, nrow = nsample)
    for(k in 1:M){
      uni = unique(res[k,])
      nun = length(uni)
      res_k=res[k,]
      for (ii in 1:nun) 
        ma[res[k,] == uni[ii], res_k ==  uni[ii]] = ma[res_k == uni[ii], res_k == uni[ii]] + 1
    }
    ma = ma/M
    Edist = as.matrix(dist(data))
    ma[ma < epsilon] = 0

#  Entropy calculation
#    y = ma
#    diag(y) = NA
#    yy = as.numeric(y)
#    yy = yy[!is.na(yy)]
#    yy = yy/sum(yy)
#    H = -sum(ifelse(yy > 0, yy * log(yy), 0))
    
    mam = (1/ma) * Edist
    mam[is.na(mam)] <- .Machine$double.xmax
    mam[is.infinite(mam) & mam > 0] <- .Machine$double.xmax
    mam = floyd(mam)
    mam[mam == .Machine$double.xmax] <- NA
    prox = Edist/mam
    diag(prox) = 1
    prox[is.na(prox)] = 0
    maxvalue = max(mam, na.rm = TRUE)
    mam[is.na(mam)] = maxvalue

    dissimilarity = mam
  }

  knn_Armadillo = knn_Armadillo(data, data, neighbors + 1)
  knn_Armadillo$distances = knn_Armadillo$distances[, -1]
  knn_Armadillo$nn_index = knn_Armadillo$nn_index[, -1]
  for (i_tsne in 1:nrow(data)) {
    for (j_tsne in 1:neighbors) {
      kod_tsne = mean(res[, i_tsne] == res[, knn_Armadillo$nn_index[i_tsne, j_tsne]], na.rm = TRUE)
      knn_Armadillo$distances[i_tsne, j_tsne] = knn_Armadillo$distances[i_tsne,  j_tsne]/kod_tsne
    }
    oo_tsne = order(knn_Armadillo$distance[i_tsne, ])
    knn_Armadillo$distances[i_tsne, ] = knn_Armadillo$distances[i_tsne, oo_tsne]
    knn_Armadillo$nn_index[i_tsne, ] = knn_Armadillo$nn_index[i_tsne, oo_tsne]
  }

  knn_Armadillo$neighbors = neighbors
  return(list(dissimilarity = dissimilarity, acc = accu, proximity = ma, 
              v = vect_acc, res = res, 
              knn_Armadillo = knn_Armadillo, 
              data = data))
}


