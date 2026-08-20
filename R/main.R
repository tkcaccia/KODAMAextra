

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



                   
plot_slide =  function (xy, slide, labels, col = NULL, nrow = 1, scales = "free",size.dot = 3,
                        size.legend.text=10,size.legend.title=10,size.legend.dot=5,size.strip.text=10)
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
  ggplot(df, aes(x, y, color = labels)) + geom_point(size = size.dot) +
    facet_wrap(~slide, nrow = nrow, scales = scales) + theme_bw() +
    theme(legend.position = "bottom", axis.title = element_blank(),
          legend.text = element_text(size = size.legend.text),  
          legend.title = element_text(size = size.legend.title), 
          strip.text = element_text(size = size.strip.text, face = "bold"), 
          axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank()) +
    scale_color_manual("Domain",  values = col, drop = FALSE) + guides(color = guide_legend(nrow = 1,  override.aes = list(size = size.legend.dot)))  # Legend
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





