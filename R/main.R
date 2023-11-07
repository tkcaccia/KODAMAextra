
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
function (data, M = 100, Tcycle = 20, FUN_VAR = function(x) {
  ceiling(ncol(x))
}, FUN_SAM = function(x) {
  ceiling(nrow(x) * 0.75)
}, bagging = FALSE, FUN = c("PLS-DA", "KNN"), f.par = 5, W = NULL, 
constrain = NULL, fix = NULL, epsilon = 0.05, dims = 2, landmarks = 10000, 
neighbors = min(c(landmarks, nrow(data)/3)) + 1, spatial = NULL, 
spatial.knn = 10, splitting = 50, clust_contrain = FALSE, n.cores=1) 
{
  if (is.null(spatial)) {
    spatial = data
    spatial_flag = TRUE
  }
  else {
    spatial_flag = FALSE
  }
  if (sum(is.na(data)) > 0) {
    stop("Missing values are present")
  }
  if (is.null(fix)) 
    fix = rep(FALSE, nrow(data))
  if (is.null(constrain)) 
    constrain = 1:nrow(data)
  data = as.matrix(data)
  shake = FALSE
  nsample = nrow(data)
  landpoints = NULL
  nlandmarks = landmarks
  if (length(landmarks) > 1) {
    if (max(landmarks) > nsample) {
      stop("A selected landmark exceed the number of entries")
    }
    if (length(table(table(landmarks))) > 1) {
      stop("Repeated landmarks are not allowed")
    }
    if (length(landmarks) > nsample) {
      stop("The number of landmarks exceed the number of entries")
    }
    nlandmarks = length(landmarks)
  }
  LMARK = (nsample > nlandmarks)
  if (LMARK) {
    if (length(landmarks) > 1) {
      landpoints = landmarks
    }
    else {
      landpoints = sort(sample(nrow(data), landmarks))
      clust = as.numeric(kmeans(data, landmarks)$cluster)
      landpoints = NULL
      for (ii in 1:landmarks) {
        www = which(clust == ii)
        if (length(www) == 1) {
          landpoints = c(landpoints, www)
        }
        else {
          landpoints = c(landpoints, sample(www)[1])
        }
      }
    }
    Tdata = data[-landpoints, , drop = FALSE]
    Xdata = data[landpoints, , drop = FALSE]
    Xdata_landpoints = Xdata
    Tfix = fix[-landpoints]
    Xfix = fix[landpoints]
    Tconstrain = constrain[-landpoints]
    Xconstrain = constrain[landpoints]
    vect_proj = matrix(NA, nrow = M, ncol = nrow(Tdata))
    Xspatial = spatial[landpoints, , drop = FALSE]
    Tspatial = spatial[-landpoints, , drop = FALSE]
  }
  else {
    Xdata = data
    Xdata_landpoints = Xdata
    Xfix = fix
    Xconstrain = constrain
    landpoints = 1:nsample
    Xspatial = spatial
    Tspatial = NULL
  }
  nva = ncol(Xdata)
  nsa = nrow(Xdata)
  res = matrix(nrow = M, ncol = nsa)
  ma = matrix(0, ncol = nsa, nrow = nsa)
  normalization = matrix(0, ncol = nsa, nrow = nsa)
  FUN_VAR = FUN_VAR(Xdata)
  FUN_SAM = FUN_SAM(Xdata)
  if (f.par > FUN_VAR & FUN[1] == "PLS-DA") {
    message("The number of components selected for PLS-DA is too high and it will be automatically reduced to ", 
            FUN_VAR)
    f.par = FUN_VAR
  }
  if (f.par > FUN_VAR & FUN[1] == "KNNPLS-DA") {
    message("The number of components selected for PLS-DA is too high and it will be automatically reduced to ", 
            FUN_VAR)
    f.par = FUN_VAR
  }
  vect_acc = matrix(NA, nrow = M, ncol = Tcycle)
  accu = NULL
  whF = which(!Xfix)
  whT = which(Xfix)
  FUN_SAM = FUN_SAM - length(whT)
#  pb <- txtProgressBar(min = 1, max = M, style = 1)
  
  
  my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "PSOCK"
  )
  
  print(my.cluster)
  doParallel::registerDoParallel(cl = my.cluster)

  foreach::getDoParRegistered()
  foreach::getDoParWorkers()
  
  res_parallel=foreach(k=1:M) %dopar%   {
    library("KODAMA")
  #  setTxtProgressBar(pb, k)
    sva = sample(nva, FUN_VAR, FALSE, NULL)
    ssa = c(whT, sample(whF, FUN_SAM, bagging, NULL))
    if (LMARK) {
      xTdata = Tdata[, sva]
      if (spatial_flag) {
        Tspatial_ssa = Tspatial[, sva]
        Xspatial_ssa = Xspatial[ssa, sva]
      }
      else {
        Tspatial_ssa = Tspatial
        Xspatial_ssa = Xspatial[ssa, ]
      }
    }
    else {
      xTdata = NULL
      Xspatial_ssa = Xspatial[ssa, ]
      Tspatial_ssa = NULL
      if (spatial_flag) {
        Xspatial_ssa = Xspatial[, sva]
      }
    }
    x = Xdata[ssa, sva]
    xva = ncol(x)
    xsa = nrow(x)
    Xconstrain_ssa = as.numeric(as.factor(Xconstrain[ssa]))
    Xconstrain_ssa_previous = Xconstrain[ssa]
    Xfix_ssa = Xfix[ssa]
    del_n = rep(NA, nrow(x))
    for (ik in 1:(nrow(x) - 1)) {
      if (is.na(del_n[ik])) {
        del_n[ik] = ik
        for (ij in 2:nrow(x)) {
          if (all(x[ik, ] == x[ij, ])) 
            del_n[ij] = ik
        }
      }
    }
    if (is.na(del_n[nrow(x)])) 
      del_n[nrow(x)] = nrow(x)
    xsa_same_point = length(unique(del_n))
    if (is.null(W)) {
      if (xsa_same_point <= 200 || length(unique(x)) < 
          50) {
        XW = Xconstrain_ssa
      }
      else {
        clust = as.numeric(kmeans(Xspatial_ssa, splitting)$cluster)
        tab = apply(table(clust, Xconstrain_ssa), 2, 
                    which.max)
        XW = as.numeric(as.factor(tab[as.character(Xconstrain_ssa)]))
        if (clust_contrain == TRUE) {
          XW = clust
          Xconstrain_ssa = clust
        }
      }
    }
    else {
      XW = W[landpoints][ssa]
      if (any(is.na(XW))) {
        if (xsa_same_point <= 200 || length(unique(x)) < 
            50) {
          unw = unique(XW)
          unw = unw[-which(is.na(unw))]
          ghg = is.na(XW)
          nnew = length(unique(Xconstrain_ssa[ghg]))
          XW[ghg] = as.numeric(as.factor(Xconstrain_ssa[ghg])) + 
            length(unw)
        }
        else {
          clust = as.numeric(kmeans(Xspatial_ssa, splitting)$cluster)
          tab = apply(table(clust, Xconstrain_ssa), 2, 
                      which.max)
          constrain_temp = as.numeric(as.factor(tab[as.character(Xconstrain_ssa)]))
          unw = unique(XW)
          unw = unw[-which(is.na(unw))]
          ghg = is.na(XW)
          nnew = length(unique(constrain_temp[ghg]))
          XW[ghg] = as.numeric(as.factor(constrain_temp[ghg])) + 
            length(unw)
          if (clust_contrain == TRUE) {
            XW = clust
            Xconstrain_ssa = clust
          }
        }
      }
    }
    clbest = XW
    options(warn = -1)
    yatta = 0
    attr(yatta, "class") = "try-error"
    while (!is.null(attr(yatta, "class"))) {
      yatta = try(core_cpp(x, xTdata, clbest, Tcycle, FUN, 
                           f.par, Xconstrain_ssa, Xfix_ssa, shake, Xspatial_ssa, 
                           Tspatial_ssa, spatial.knn), silent = FALSE)
      if (!is.null(attr(yatta, "class"))) {
        save(yatta, x, xTdata, clbest, Tcycle, FUN, f.par, 
             Xconstrain_ssa, Xfix_ssa, shake, Xspatial_ssa, 
             Tspatial_ssa, spatial.knn, file = "/Users/stefano/Desktop/Chepalle2.RData")
      }
    }
    options(warn = 0)
    if (is.list(yatta)) {
      clbest = as.vector(yatta$clbest)
      accu = yatta$accbest
      yatta$vect_acc = as.vector(yatta$vect_acc)
      yatta$vect_acc[yatta$vect_acc == -1] = NA
      vect_acc[k, ] = yatta$vect_acc
      if (LMARK) {
        yatta$vect_proj = as.vector(yatta$vect_proj)
        yatta$vect_proj[Tfix] = W[-landpoints][Tfix]
        vect_proj[k, ] = yatta$vect_proj
      }

  #    normalization[ssa, ssa] = normalization[ssa, ssa] +  1
  #    res[k, ssa] = clbest
    }
    list(ssa=ssa,clbest=clbest)

  }
  parallel::stopCluster(cl = my.cluster)
  
  for(k in 1:M){
    ssa_i=res_parallel[[k]]$ssa
    clbest_i=res_parallel[[k]]$clbest
    
    uni = unique(clbest_i)
    nun = length(uni)
    for (ii in 1:nun) ma[ssa_i[clbest_i == uni[ii]], ssa_i[clbest_i == 
                                                       uni[ii]]] = ma[ssa_i[clbest_i == uni[ii]], ssa_i[clbest_i == 
                                                                                                    uni[ii]]] + 1
    

    
    normalization[ssa_i, ssa_i] = normalization[ssa_i, ssa_i] +  1
    res[k, ssa_i] = clbest_i
  }
#  close(pb)
  ma = ma/normalization
  Edist = as.matrix(dist(Xdata_landpoints))
  ma[ma < epsilon] = 0
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
  y = ma
  diag(y) = NA
  yy = as.numeric(y)
  yy = yy[!is.na(yy)]
  yy = yy/sum(yy)
  H = -sum(ifelse(yy > 0, yy * log(yy), 0))
  dissimilarity = mam
  if (LMARK) {
    total_res = matrix(nrow = M, ncol = nsample)
    total_res[, landpoints] = res
    total_res[, -landpoints] = vect_proj
    knn_Armadillo = knn_Armadillo(data, data, neighbors + 
                                    1)
    knn_Armadillo$distances = knn_Armadillo$distances[, -1]
    knn_Armadillo$nn_index = knn_Armadillo$nn_index[, -1]
    for (i_tsne in 1:nrow(data)) {
      for (j_tsne in 1:neighbors) {
        kod_tsne = mean(total_res[, i_tsne] == total_res[, 
                                                         knn_Armadillo$nn_index[i_tsne, j_tsne]], na.rm = TRUE)
        knn_Armadillo$distances[i_tsne, j_tsne] = knn_Armadillo$distances[i_tsne, 
                                                                          j_tsne]/kod_tsne
      }
      oo_tsne = order(knn_Armadillo$distance[i_tsne, ])
      knn_Armadillo$distances[i_tsne, ] = knn_Armadillo$distances[i_tsne, 
                                                                  oo_tsne]
      knn_Armadillo$nn_index[i_tsne, ] = knn_Armadillo$nn_index[i_tsne, 
                                                                oo_tsne]
    }
  }
  else {
    knn_Armadillo = list()
    knn_Armadillo$nn_index = matrix(ncol = ncol(mam), nrow = nrow(mam))
    for (i_tsne in 1:nrow(data)) {
      oo_tsne = order(mam[i_tsne, ])
      mam[i_tsne, ] = mam[i_tsne, oo_tsne]
      knn_Armadillo$nn_index[i_tsne, ] = oo_tsne
    }
    knn_Armadillo$nn_index = knn_Armadillo$nn_index[, -1][, 
                                                          1:neighbors]
    knn_Armadillo$distances = mam[, -1][, 1:neighbors]
    total_res = res
  }
  knn_Armadillo$neighbors = neighbors
  return(list(dissimilarity = dissimilarity, acc = accu, proximity = ma, 
              v = vect_acc, res = total_res, f.par = f.par, entropy = H, 
              landpoints = landpoints, knn_Armadillo = knn_Armadillo, 
              data = data))
}















                              
KODAMA.matrix.parallel =
function (data,                       # Dataset
          spatial = NULL,             # In spatial are conteined the spatial coordinates of each entries
          M = 100, Tcycle = 20, 
          FUN = c("PLS","PK","KNN"), 
          f.par.knn = 5, f.par.pls = 5,
          W = NULL, 
          constrain = NULL, fix = NULL, epsilon = 0.05, landmarks = 10000,  
          splitting = 50, spatial.resolution = 0.3, ncore = 1) 
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
  res_parallel=foreach(k=1:M) %dopar%   {
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
    res_k=NULL
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
#  close(pb)

  

    
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
              landpoints = landpoints, knn_Armadillo = knn_Armadillo, 
              data = data))
}


