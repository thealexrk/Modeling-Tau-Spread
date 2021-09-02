source('code/misc/miscfxns.R')

##################################
### Atlas conversion functions ###
##################################

map.ABA.to.CNDR <- function(X.ABA,CNDR.names,ABA.to.CNDR.key){
  #### NEED TO ENSURE ALL REGION NAMES ARE CORRECT ####
  # Input:
  # X.ABA: named matrix of values associated with ABA regions, must have ABA names as rows, and any number of columns
  # CNDR.names: character vector of CNDR names 
  # ABA.to.CNDR.key: list whose element names are CNDR names and whose elements contain vectors of corresponding ABA names
  
  # Output:
  # X.CNDR: named matrix of values associated with CNDR regions, has CNDR names for rows
  
  X.ABA <- as.matrix(X.ABA)
  X.CNDR <- matrix(NA,nrow=length(CNDR.names),ncol=ncol(X.ABA)) # make new vector to hold data, default is NA
  colnames(X.CNDR) <- colnames(X.ABA)
  rownames(X.CNDR) <- CNDR.names
  for(CNDR.region in CNDR.names){
    ABA.region.match <- ABA.to.CNDR.key[[CNDR.region]] # get ABA region(s) corresponding to that CNDR region, if any
    if(length(ABA.region.match) > 0){      # print out names in key that aren't in ABA connectome
      display.missing.items(names(X.ABA),ABA.region.match,'Designation names not in ABA connectome names:')
      X.CNDR[CNDR.region,] <- colMeans(X.ABA[rownames(X.ABA) %in% ABA.region.match,,drop=FALSE]) # average together ABA regions in each CNDR region, for all columns
    } else if(length(ABA.region.match) ==0){# otherwise if no match, leave as NA    
      print(paste('no match for CNDR region',CNDR.region))
    }
  }  
  return(X.CNDR)
}

display.missing.items <- function(ref,items,msg=''){
  # find and print the elements of items that are not in ref
  # msg: extra text to display before printing missing items
  if(length(items) > 0){
    missingitems <- !(items %in% ref) 
    if(sum(missingitems) > 0){ # if there are any missing items
      print(msg)
      print(items[missingitems])
    }
  } else {print('WARNING: items is empty')}
}

map.CNDR.to.ABA <- function(X.CNDR,ABA.names,CNDR.to.ABA.key){
  # Input:
  # X.CNDR: named vector of values associated with CNDR regions, must have CNDR names
  # ABA.names: character vector of ABA names
  # CNDR.to.ABA.key: list whose element names are CNDR names and whose elements contain vectors of corresponding ABA names
  
  # Output:
  # X.ABA: named vector of values associated with ABA regions, has ABA names
  if(is.matrix(X.CNDR)){ # convert X.CNDR to named 1D vector if it's a named matrix
    if(nrow(X.CNDR) == 1){X.CNDR <- setNames(X.CNDR,colnames(X.CNDR))}
    if(ncol(X.CNDR) == 1){X.CNDR <- setNames(X.CNDR,rownames(X.CNDR))}
  }

  X.ABA <- matrix(NA,nrow=length(ABA.names)) # make new vector to hold data, default is NA
  names(X.ABA) <- ABA.names
  for(ABA.region in ABA.names){
    CNDR.region.match <- CNDR.to.ABA.key[[ABA.region]] # get ABA region(s) corresponding to that CNDR region, if any
    if(length(CNDR.region.match) > 0){      # print out names in key that aren't in ABA connectome
      display.missing.items(names(X.CNDR),CNDR.region.match,'Designation names not in CNDR path data names:')
      X.ABA[ABA.region] <- mean(X.CNDR[CNDR.region.match])
    } else if(length(CNDR.region.match) ==0){# otherwise if no match, leave as NA    
      print(paste('no match for ABA region',ABA.region))
    }
  }  
  return(X.ABA)
}

bootstrap.path.tps <- function(Mice,nboot=1000){
  # INPUTS:
  # Mice: list of dataframes that is mice x regions whose values contain pathology
  # nboot: number of bootstraps
  #
  # OUTPUTS:
  # List with nboot elements, each of which contains as many elements as Mice. These sub elements contained bootstrapped
  # log10 mean pathology for each element of Mice
  
  n.mice <- lapply(Mice,nrow)
  X.boot <- lapply(1:nboot, function(n) # for nboot bootstraps, make a list of 3 single bootstraps for each time point
    lapply(1:length(Mice), function(t) log(colMeans(Mice[[t]][sample(1:n.mice[[t]],replace=T),],na.rm = T),base=10)))
  return(X.boot)
}

###############################
### fitting diffusion model ###
###############################

get.Xo <- function(region.names,ROIs){
  # generate initial "pathology seed" vector, which is all 0's except
  # a 1 in all regions listed in character vector ROI
  # ROIs should correspond to region.names
  n.regions <- length(region.names)
  Xo <- matrix(0, nrow=n.regions,ncol=1)
  Xo[which(region.names %in% ROIs)] <- 1
  rownames(Xo) <- region.names
  return(Xo)
}

get.Lout <- function(W,S,ant.ret='retro'){
  # INPUTS:
  # W: NxN adjacency matrix, may be asymmetric (if so, connections are from row i to column j)  
  # S: vector of weights to apply to each row of W (here, synuclein expression)
  # ant.ret: character specifying anterograde or retrograde spread
  #
  # The connection matrix's default directionallity here is opposite the convention for matrix multiplication
  # so we are capturing retrograde spread along connections
  # compute out-degree laplacian for retrograde spread
  #
  # OUTPUTS:
  # return specified laplacian

  n.regions <- nrow(W)
  W <- W * !diag(n.regions) # get rid of diagonal
  mag.orig <- sum(W) # original sum of weights
  W <- diag(as.numeric(S)) %*% W
  #W <- W * mag.orig/sum(W) # scale new W to have the same global strength... lets you fit with similar c range
  #W <- W / (max(Re(eigen(W)$values))) # scale to max eigenvalue

  # Where i is row element and j is column element
  # Wij is a connection from region i to region j
  # convention is the opposite, so without transposing W
  # I am capturing "retrograde" connectivity

  in.deg <- colSums(W)
  out.deg <- rowSums(W)
  if(ant.ret=='retro'){L.out <- diag(x = out.deg) - W} # outdegree laplacian
  if(ant.ret=='antero'){L.out <- diag(x = in.deg) - t(W)} #indegree laplacian and transposed matrix for anterograde
  return(L.out)
}

predict.Lout <- function(L,Xo,c,t=1,fxn=reticulate::import('scipy.linalg')$expm){
  # Generate prediction for given time points using
  # L: Laplacian of adjacency matrix
  # Xo: initial conditions
  # c: time constant
  # t: time points
  # fxn: function for matrix exponential. default is python's scipy.linalg.expm 
  # using reticulate for faster performance (not as fast as matlab or python native though)
  # but R would have you use expm::expm() which is v slow

  Xt <- do.call('cbind', lapply(t, function(t.i) fxn(-L*c*t.i)%*%Xo))
  rownames(Xt) <- rownames(Xo)
  return(Xt)
}

c.fit <- function(log.path,L.out,Xo,c.rng){
  # INPUTS:
  # log.path: vector of log-10 transformed pathology scores. Time constant c is fit to predict these
  # L.out: out-degree graph laplacian of anatomical connectivity matrix
  # Xo: vector of initial pathology
  # c.rng: range of time constants to test

  # OUTPUTS:
  # c: optimal time constant within c.rng for predicting log.path

  # exclusion mask... don't count regions with 0 path
  mask <- log.path != -Inf
  # compute fit at each time point for range of time
  Xt.sweep <- sapply(c.rng, function(c) # no time here, because in this project there's only 1 time pointo 
    cor(log.path[mask],log(predict.Lout(L.out,Xo,c),base=10)[mask]))
  c <- c.rng[which.max(Xt.sweep)] # select c giving max correlation
  print(c)

  return(c)
}

c.CNDRspace.fit <- function(log.path,tps,L.out,Xo,c.rng,ABA.to.CNDR.key,excl.inj.CNDR){
  # fits time constant by modeling CNDR data
  # INPUTS:
  # log.path: list of vectors of log-10 transformed pathology scores *in CNDR space*  
  # for each time point specified in tps (below). Time constant c is fit to predict these
  # tps: vector of numeric time points post injection
  # L.out: out-degree graph laplacian of anatomical connectivity matrix
  # Xo: vector of initial pathology
  # c.rng: range of time constants to test
  # excl.inj.CNDR: vector of injection sites in CNDR space to exclude when computing fit

  # this function runs diffusion model in ABA space, 
  # but computes correlation with real data in CNDR annotation space to assess fit of time constant c
  # CNDR names are names of log.path, ABA names are names of Xo
  # OUTPUTS:
  # c: optimal time constant within c.rng for predicting log.path
  
  # exclusion mask... don't count regions with 0 path
  scipy.linalg <- reticulate::import('scipy.linalg') # import python function for matrix exponential because it's faster
  #log.path <- lapply(1:length(log.path), function(t) log.path[[t]][mask[[t]]])
  # compute fit at each time point for range of time
  Xt.sweep <- matrix(NA,nrow=length(c.rng),ncol=length(tps))
  
  if(!is.null(excl.inj.CNDR)){ # exclude injection sites from fit by setting them to -Inf in log.path so they'll be excluded by cor.mask below
    for(t. in 1:length(tps)){log.path[[t.]][excl.inj.CNDR] <- -Inf}
  }

  for(c.i in 1:length(c.rng)){
  print(paste0('c ',c.i,' out of ',length(c.rng)))
    c.val <- c.rng[c.i]
    #ptm <- proc.time()
    Xt.c <- do.call('cbind',lapply(tps,function(t) predict.Lout(L.out,Xo,c.val,t,fxn=scipy.linalg$expm))) # predict path using linear diff model into matrix that is region-by-time
    #print(paste0('predict: ',(proc.time()-ptm)['elapsed']))
    #ptm <- proc.time()
    Xt.c <- quiet(map.ABA.to.CNDR(Xt.c,names(log.path[[1]]),ABA.to.CNDR.key)) # convert matrix to CNDR space
    #print(paste0('map: ',(proc.time()-ptm)['elapsed']))
    Xt.sweep[c.i,] <- sapply(1:length(tps), function(t) cor.mask(Xt.c[,t],log.path[[t]])) # this replaces the two lines below
    #mask <- lapply(1:length(tps), function(t) is.na(Xt.c[,t]) | mask[[t]])
    #lapply(1:length(tps), function(t) print(paste('missing CNDR region',names(Xt.c[[t]])[!mask[[t]]],'from ABA connectome')))
    #Xt.sweep[c.i,] <- sapply(1:length(tps), function(t) cor(log.path[[t]][!mask[[t]]],log(Xt.c[!mask[[t]],t],base=10)))
  }
  c.best <- c.rng[which.max(rowMeans(Xt.sweep))] # select c giving max correlation
  print(c.best)

  return(list(c.best=c.best,Xt.sweep=Xt.sweep))
}

cor.mask <- function(Xt,y){
  # INPUTS:
  # Xt: predicted pathology
  # y: observed log pathology
  #
  # OUTPUTS:
  # compute pearson correlation between log predicted (computed) and log observed (input)
  # excluding elements that were originally 0 and thus log(0) = -Inf
  Xt <- log(Xt,base=10)

  mask <- y == -Inf | Xt == -Inf | is.na(Xt)
  return(cor(y[!mask],Xt[!mask]))
}

c.ABAspace.fit <- function(log.path,L.out,Xo,c.rng,CNDR.to.ABA.key){
  # fits time constant by modeling CNDR data
  # INPUTS:
  # log.path: vector of log-10 transformed pathology scores *in CNDR space*. Time constant c is fit to predict these
  # L.out: out-degree graph laplacian of anatomical connectivity matrix
  # Xo: vector of initial pathology
  # c.rng: range of time constants to test

  # this function runs diffusion model in ABA space, 
  # abd computes correlation with real data in ABA annotation space to assess fit of time constant c

  # OUTPUTS:
  # c: optimal time constant within c.rng for predicting log.path
  
  log.path <- quiet(log(map.CNDR.to.ABA(10^(log.path),rownames(Xo),CNDR.to.ABA.key),base=10)) # undo log, convert to ABA space, retake log
  # exclusion mask... don't count regions with 0 path
  mask <- log.path != -Inf & !is.na(log.path)
  # compute fit at each time point for range of time
  Xt.sweep <- sapply(c.rng, function(c) # no time here, because in this project there's only 1 time pointo 
    cor(log.path[mask],log(predict.Lout(L.out,Xo,c),base=10)[mask]))
  c <- c.rng[which.max(Xt.sweep)] # select c giving max correlation
  print(c)

  return(c)
}

