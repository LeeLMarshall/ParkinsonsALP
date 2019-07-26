# library(methylKit)

computeSVs <- function(Y, model, model0, n.sv = NULL, devel = FALSE, controls = NULL) {
  require(sva)

  # filter out loci where variance is 0
  v <- matrixStats::rowVars(Y, na.rm=TRUE)
  lidx <- v > 0
  Y <- Y[lidx, ]
  v <- v[lidx]

  if (devel) {
    message("DEVEL: computeSVs chosing most variable 5k loci")
    i <- order(v, decreasing = TRUE)[1:5000]
    Y <- Y[i,]
  }
  if (!is.numeric(n.sv)) {
    n.sv <- num.sv(na.omit(Y), model)
    message("computeSVs: using num.sv to guess number of SVs:", n.sv)
  } else {
    message("computeSVs: fixed number of n.sv:", n.sv)
  }
  if (is.null(controls)) {
    svobj <- sva(na.omit(Y), model, model0, n.sv = n.sv)  
  } else {
    controls <- controls[lidx]
    svobj <- sva(Y, model, model0, n.sv = n.sv, controls = controls)
  }
  
  res <- svobj$sv
  n <- ncol(res)
  colnames(res) <- paste("SVA", 1:n, sep="")
  res
}

permuteModelData <- function(modelData, permute, i = NULL) 
{
  modelData <- as.data.frame(modelData)
  if (is.null(i)) {
    i <- sample(nrow(modelData))
  }
  if (is.logical(permute)) {
    if (permute == TRUE) {
      message("permuteModelData: permuting all variables")
      modelData <- modelData[i,, drop=FALSE]
    } else {
      message("permuteModelData: not permuting")
    }
  } else {
    permute <- permute[ permute %in% colnames(modelData)]
    message("permuteModelData: shuffling variables ", 
      paste(permute, collapse=", "))
    modelData[, permute] <- modelData[i, permute, drop=FALSE]
  }
  return(modelData)
}



getBeta <- function(dataObj) {
  message("Processing count matrix")
  counts <- as(dataObj$counts, "methylRawList")
  counts <- unite(counts, destrand=FALSE, 
    mc.cores = params$workers, min.per.group = 1L)
  counts <- S3Part(counts, strictS3 = TRUE)

  # select T and C columns
  Tcols <- seq(7,ncol(counts),by=3)
  Ccols <- Tcols-1
  ids <- paste(counts[,1], counts[,2], sep="_")
  counts <- as.matrix(counts[,c(Ccols,Tcols)])
  rownames(counts) <- ids

  # get the beta value matrix
  y <- counts[, 1:(ncol(counts)/2)]
  w <- y+counts[, ((ncol(counts)/2)+1):ncol(counts)]
  beta <- y/w
  rm(w, y, Tcols, Ccols, ids)
  return(beta)  
}


# computeFitsLM <- function(dataObj, useData, useSVA,
#   formula, formula0, modelData, permute, 
#   comBatVariable = NULL,
#   ...) 
# {
#   require(sva)
#   require(limma)

#   modelData <- as.data.frame(modelData)

#   # Handle the permute
#   modelData <- permuteModelData(modelData, permute)

#   # Get the data object
#   stopifnot(useData %in% c("b", "m"))
#   Y <- dataObj$padlock[[useData]]

#   # Fix batch effect using ComBat
#   if (!is.null(comBatVariable))
#     Y <- sva::ComBat(Y, comBatVariable)

#   # Build the models
#   model <- model.matrix(as.formula(formula), modelData)
#   model0<- model.matrix(as.formula(formula0), modelData)

#   # Add SVA
#   if (useSVA) {
#     svs <- computeSVs(Y, model, model0, n.sv=useSVA, ...)
#     model <- cbind(model, svs)
#     model0<- cbind(model0, svs)
#   }

#   # Compute the p values
#   f.pvalue(Y, model, model0)
# }




# estimatePhi<-function(counts, modelMat, useBayes = FALSE){
#   i <- is.na(counts[1:nrow(modelMat)])
#   modelMat <- modelMat[!i,]
#   counts <- counts[!is.na(counts)]
  
#   n=counts[1:(length(counts)/2)]+counts[((length(counts)/2)+1):length(counts)]
#   y=counts[1:(length(counts)/2)]
#   prop=y/n

#   if (useBayes) {
#     require(arm)
#     glmfit <- bayesglm.fit(modelMat, prop, weights=n, family=binomial(link=logit), 
#       prior.df = 10)
#   } else {
#     glmfit <- glm.fit(modelMat, prop, weights=n, family=binomial(link=logit))
#   }
  
#   # fit glm
#   mu <- fitted(glmfit)
#   nprm=length(glmfit$coef) # number of parameters fitted
  
#   # calculate and record results
#   resids <- (y-n*mu)/sqrt(mu*(n-n*mu))

#   # get phi correction coefficients 
#   phi <- sum( resids^2 )/(length(n)-nprm)
#   c(phi,(length(n)-nprm))
# }

# estimateShrinkageMN<-function(cntlist, modelData, formula, 
#                               sample.size=100000, workers = 40,
#                               ...)
# {  
#   # get formula and construct model matrix
#   modelMat<-model.matrix(as.formula(formula), modelData)
  
#   # check number of samples
#   sample.size <- ifelse(length(cntlist)<=sample.size,length(cntlist),sample.size)
  
#   # calculate phis up to the first 100000 sites (stabilizing point)
#   if (workers > 1) {
#     estimation <- mclapply(cntlist[1:sample.size], estimatePhi,
#       modelMat=modelMat, mc.cores = workers, ...)
#   } else {
#     estimation <- lapply(cntlist[1:sample.size], estimatePhi,
#       modelMat=modelMat, ...)
#   }
#   estimation=simplify2array(estimation)

#   # for each phi, take the correct df (depending on number of model parameters)
#   phis<-estimation[1,]
#   df<-estimation[2,]
  
#   # squeeze sample variances 
#   shr <- limma::squeezeVar(phis,df)
  
#   # output prior df and variances (to be used later as input for parShrinkMN)
#   res <- list(df.prior=shr$df.prior,var.prior=shr$var.prior)
#   message("Shrinkage params:", res)
#   return(res)
# }

# myglm.fit <- function(model, prop, w, useBayes = FALSE) {
#   if (!useBayes) {
#     glm.fit(model, prop, 
#         weights=as.integer(w), family=binomial(link=logit))
#   } else {
#     require(arm)
#     bayesglm.fit(model, prop, 
#         weights=as.integer(w), 
#         family=binomial(link=logit), prior.df = 10)
#   }
# }

# # apply logistic regression
# logReg<-function(counts, formula, formula0 = NULL, modelData,
#                  overdispersion=c("none","MN","shrinkMN"),
#                  effect=FALSE, 
#                  parShrinkMN=list(), 
#                  test=c("F","Chisq"),
#                  useBayes = FALSE)
# {
#   tryCatch({
#     i <- is.na(counts[1:nrow(modelData)])
#     modelData <- modelData[!i,,drop=FALSE]
#     counts <- counts[!is.na(counts)]


#     w=counts[1:(length(counts)/2)]+counts[((length(counts)/2)+1):length(counts)]
#     y=counts[1:(length(counts)/2)]
#     prop=y/w

#     # Null model
#     if(!is.null(formula0)){
#       modelCov<-model.matrix(as.formula(formula0), modelData)
#       objCov <- myglm.fit(modelCov, prop, w, useBayes=useBayes)
#     }

#     # full model with all variables
#     modelMat<-model.matrix( as.formula(formula), as.data.frame(modelData) )
#     obj <- myglm.fit(modelMat, prop, w, useBayes=useBayes)

#     mu=fitted(obj)
#     nprm=length(obj$coef) # number of parameters fitted

#     #get dispersion
#     overdispersion <- match.arg(overdispersion)
#     phi=switch(overdispersion,
#                none=1,
#                MN={
#                  mu=fitted(obj)
#                  uresids <- (y-w*mu)/sqrt(mu*(w-w*mu)) # pearson residuals
#                  phi=sum( uresids^2 )/(length(w)-nprm) # correction factor  
#                  ifelse(phi>1,phi,1)
#                },
#                shrinkMN={
#                  mu=fitted(obj)
#                  uresids <- (y-w*mu)/sqrt(mu*(w-w*mu)) # pearson residuals
#                  phi=sum( uresids^2 )/(length(w)-nprm) # correction factor
                 
#                  # change phi with parameters from parShrinkMN
#                  df.prior=parShrinkMN$df.prior
#                  var.prior=parShrinkMN$var.prior
#                  df.total=(length(w)-nprm)+df.prior
#                  phi=((length(w)-nprm)*phi + df.prior*var.prior)/df.total
#                  ifelse(phi>1,phi,1)
#                })

#     if(!is.null(formula0)){
#       deviance <- objCov$deviance - obj$deviance
#       ddf=objCov$df.residual-obj$df.residual
#     }else{
#       deviance <- obj$null.deviance - obj$deviance
#       ddf=obj$df.null-obj$df.residual # difference in degrees of freedom
#     }

#     test=match.arg(test)

#     # do F-test when overdispersion >1 given
#     test=ifelse(test=="F" & phi>1,"F","Chisq")

#     p.value=switch(test,
#                    F={
#                      pf(deviance/phi, ddf, (length(w)-nprm), lower.tail = FALSE)     
#                    },
#                    Chisq={
#                      pchisq(deviance/phi, 1, lower.tail = FALSE)
#                    })
#     if (effect == TRUE) {
#       vars <- modelMat[, !colnames(modelMat) %in% colnames(modelCov)[-1], drop=FALSE]
#       # classMeans <- c(`(Intercept)` = mean(mu[apply(!vars, 1, all)]),
#       #                 apply(vars, 2, function(X) mean(mu[X])))
#       effect <-
#         lapply(1:ncol(vars), function(i) ifelse(vars[,i] == 1, colnames(vars)[i], "")) %>% 
#         do.call("cbind", .) %>% 
#         apply(., 1, paste, collapse="+") %>% 
#         data.table(mu, G=.) %>%
#           .[, mean(mu), by=G ] %>%
#         .[, G := gsub("[+]+", "+", G)] %>%
#         .[, G := gsub("[+]$", "", G)] %>%
#         .[order(G)]
#       classMeans <- as.numeric(effect$V1)
#       names(classMeans) <- effect$G
#       return( c(classMeans, p.value=p.value) )
#     } else {
#       return(p.value)
#     }
#   },
#   error = function(e) {
#     message(e)
#     return(NA)
#   })
# }

# computeFitsMethylKit <- function(dataObj, useSVA, 
#   formula, formula0, modelData, 
#   permute, effect = FALSE, devel=FALSE,
#   workers = 40, filter = TRUE, ...) 
# {
#   require(methylKit)

#   # Alternatively
#   modelData <- as.data.frame(modelData)

#   # Handle the permute
#   modelData <- permuteModelData(modelData, permute)

#   # Get counts
#   if (class(dataObj) == "list") {
#     message("Processing count matrix")
#     counts <- as(dataObj$counts, "methylRawList")
#     counts <- unite(counts, destrand=FALSE, 
#       mc.cores = workers, min.per.group = 1L)    
#     counts <- S3Part(counts, strictS3 = TRUE)
#   } else if (class(dataObj) == "methylBase") {
#     counts <- S3Part(dataObj, strictS3 = TRUE)
#   } else {
#     stop("Unknown type of dataObj.")
#   }
  

#   # select T and C columns
#   Tcols <- seq(7,ncol(counts),by=3)
#   Ccols <- Tcols-1
#   ids <- paste(counts[,1], counts[,2], sep="_")
#   counts <- as.matrix(counts[,c(Ccols,Tcols)])
#   rownames(counts) <- ids

#   # get the beta value matrix
#   y <- counts[, 1:(ncol(counts)/2)]
#   w <- y+counts[, ((ncol(counts)/2)+1):ncol(counts)]
#   beta <- y/w
#   rm(w, y, Tcols, Ccols, ids)

#   if (is.logical(filter) && filter == TRUE) {
#     message("Filtering loci.")
#     # filter out loci where variance is 0 or 
#     # NAs, 0s and 100s constitute more than 70% of samples
#     n <- ncol(beta)
#     v <- matrixStats::rowVars(beta, na.rm=TRUE)
#     zero <- min(beta, na.rm=TRUE)
#     one <- max(beta, na.rm=TRUE)
#     nas <- apply(beta, 1, function(X) {
#       sum(is.na(X) | X == zero | X == one)
#     })
#     lidx <- v == 0 | nas / n >= 0.7
#     beta <- beta[!lidx, ]
#     counts <- counts[!lidx, ]
#   }

#   # Compute SVA vectors
#   if (useSVA) {
#     model <- model.matrix(as.formula(formula), modelData)
#     model0<- model.matrix(as.formula(formula0), modelData)    
#     svs <- computeSVs(beta, model, model0, n.sv=useSVA, devel = devel)
#     modelData <- cbind(modelData, svs)
#   }


#   # Apply MethylKit
#   cntlist <- split(counts, 1:nrow(counts))
#   names(cntlist) <- rownames(counts)
#   parShrinkMN <- estimateShrinkageMN(cntlist, modelData, 
#     formula = formula, 
#     sample.size = ifelse(devel==TRUE, 1000, 100000),
#     workers = workers,
#     ...)
#   if (is.numeric(filter)) {
#     message("Reducing the computed loci to specified ones.")
#     cntlist <- cntlist[filter]    
#   }
#   message("Computing logistic regression.")
#   res <- mclapply(cntlist, logReg,
#     formula = formula, 
#     formula0 = formula0, 
#     modelData = modelData, 
#     overdispersion = "shrinkMN", parShrinkMN = parShrinkMN,
#     test = "Chisq", effect = effect, mc.cores = workers,
#     ...) 
#   names(res) <- names(cntlist)
#   return(simplify2array(res))
# }


# computeFits <- function(
#   dataObj, useData, useSVA, formula, formula0, modelData, permute, workers, ...
#   ) {
#   if (useData %in% c("m", "b")) {
#     return(computeFitsLM(dataObj, 
#       useData, useSVA, formula, formula0, modelData, permute, ...))
#   } else {
#     return(computeFitsMethylKit(dataObj, 
#       useSVA, formula, formula0, modelData,
#       permute, devel=FALSE, workers=workers, ...))
#   }
# }


# computeFitsWithPermutation <- function(
#   dataObj,
#   useData = c("b", "m", "counts"),
#   useSVA = FALSE,
#   nPerm = 40,
#   workers = 41,
#   formula, formula0, modelData,
#   shuffleVariables,
#   ...
#   ) 
# {
#   require(doFuture)
#   registerDoFuture()
#   if (workers > 1) {
#     options(future.globals.maxSize= 2^30)
#     plan(future::cluster, workers=workers)
#     message("Cluster future initialized.")
#   } else {
#     plan(future::sequential)
#     message("Sequential future initialized.")
#   }

#   useData <- match.arg(useData)

#   res <- foreach(i = 0:nPerm, .combine=cbind) %dopar% {
#     message(i)
#     if (i == 0)
#       computeFits(dataObj, useData, useSVA, formula, formula0, modelData,
#                   permute = FALSE, workers = 1, ...)
#     else 
#       computeFits(dataObj, useData, useSVA, formula, formula0, modelData,
#                   permute = shuffleVariables, workers = 1, ...)
#   }
#   return(res)
# }



withWarnings <- function(expr) {
  myWarnings <- NULL
  myErrors <- NULL
  wHandler <- function(w) {
    myWarnings <<- c(myWarnings, list(w))
    invokeRestart("muffleWarning")
  }
  eHandler <- function(e) {
    myErrors <<- c(myErrors, list(e))
    return(NA)   
  }
  val <- withCallingHandlers(try(expr, silent=TRUE), 
    warning = wHandler, error = eHandler)
  # list(value = val, warnings = myWarnings, errors = myErrors)
  val$warnings <- myWarnings
  val$errors <- myErrors
  return(val)
}


estimatePhi <- function(fit) {
  mu <- fitted(fit)
  nprm <- length(fit$coef) # number of parameters fitted
  n <- fit$prior.weights # total coverage
  y <- fit$y * n # beta values
  
  # calculate and record results
  resids <- (y-n*mu)/sqrt(mu*(n-n*mu))

  # get phi correction coefficients 
  sum( resids^2 )/(length(n)-nprm)
}

# This should only return useful information
singleGlmFit <- function(b, w, model, estimatePhi=FALSE) {
  res <- withWarnings({
    i <- is.na(b) | is.na(w)
    fit <- glm.fit(model[!i,], b[!i], weights = w[!i], family = binomial(link=logit))
    list(
      phi=ifelse(estimatePhi==TRUE, estimatePhi(fit), NA), 
      deviance=fit$deviance, 
      df.residual=fit$df.residual, 
      nprm=length(fit$coef),
      mu = fitted(fit)
    )
  })
  return(res)
}


fitLocus <- function(b, w, model, model0) 
{
  fit <-  singleGlmFit(b, w, model, estimatePhi = TRUE)
  fit0 <- singleGlmFit(b, w, model0, estimatePhi = FALSE)
  return(list(fit=fit, fit0=fit0, b=b, w=w))
}

# define operation on a chunk of beta and counts
multipleGlmFit <- function(X, model, model0, ...) {
  require(foreach)
  m <- ncol(X)/2
  n <- nrow(X)
  foreach (i = 1:n) %do% {
    fitLocus(
      b = X[i, 1:m],
      w = X[i, (m+1):(2*m)],
      model = model, 
      model0 = model0,
      ...
    )
  }
}

# estimate shrinkage parameters from phi of each full model
estimateShrinkageParameters <- function(res) {
  estimation <- lapply(res, function(R) c(R$fit$phi, R$fit$df.residual))
  estimation <- simplify2array(estimation)
  phis <- estimation[1,]
  df <- estimation[2,]
  shr <- limma::squeezeVar(phis, df)
  df.prior <- shr$df.prior
  var.prior <- shr$var.prior
  return(list(df.prior=df.prior, var.prior=var.prior))
}

computePValue <- function(beta, counts, fit, fit0, shrinkParams = NULL) {
  y <- na.omit(beta * counts)
  w <- na.omit(counts)
  mu <- fit$mu
  nprm <- fit$nprm
  uresids <- (y-w*mu)/sqrt(mu*(w-w*mu)) # pearson residuals
  phi <- sum( uresids^2 )/(length(w)-nprm) # correction factor

  if (!is.null(shrinkParams)) {
    df.prior <- shrinkParams[[1]]
    var.prior <- shrinkParams[[2]]
    # change phi with parameters from parShrinkMN
    df.total <- (length(w)-nprm)+df.prior
    phi <- ((length(w)-nprm)*phi + df.prior*var.prior)/df.total
    phi <- ifelse(phi>1,phi,1)
  } else {
    phi <- 1
  }
  deviance <- fit0$deviance - fit$deviance
  ddf <- fit0$df.residual-fit$df.residual
  p.value <- pchisq(deviance/phi, 1, lower.tail = FALSE)
  return(p.value)
}

computeEffect <- function(fit, model, model0) {
  require(dplyr)
  require(data.table)

  vars <- model[, !colnames(model) %in% colnames(model0)[-1], drop=FALSE]
  effect <-
    lapply(1:ncol(vars), function(i) ifelse(vars[,i] == 1, colnames(vars)[i], "")) %>% 
    do.call("cbind", .) %>% 
    apply(., 1, paste, collapse="+") %>% 
    data.table(mu=fit$mu, G=.) %>%
      .[, mean(mu), by=G ] %>%
    .[, G := gsub("[+]+", "+", G)] %>%
    .[, G := gsub("[+]$", "", G)] %>%
    .[order(G)]
  classMeans <- as.numeric(effect$V1)
  names(classMeans) <- effect$G
  classMeans
}


fastLogReg <- function(beta, counts, formula, formula0, modelData,
  effect = FALSE,
  shrinkMN = TRUE,
  chunkSize = 20000) {
  require(itertools)
  require(foreach)
  require(bigmemory)

  modelData <- as.data.frame(modelData)
  model <- model.matrix(formula %>% as.formula, modelData)
  model0 <- model.matrix(formula0 %>% as.formula, modelData)
  X <- cbind(beta, counts)
  export  <- c("multipleGlmFit", "fitLocus", "singleGlmFit", "withWarnings", 
    "estimatePhi", "computePValue", "computeEffect")

  sharedData <- as.big.matrix(X)
  desc <- describe(sharedData)
  iter <- isplitIndices(nrow(sharedData), chunkSize=chunkSize)

  res <- foreach(i=iter, .packages=c("bigmemory"), .export=export) %dopar% 
  { 
          sharedData <- attach.big.matrix(desc)
          multipleGlmFit(sharedData[i,], model=model, model0=model0)
  }
  res <- unlist(res, recursive = FALSE)
  rm(sharedData)

  if (shrinkMN) {
    message("Estimating shrinkage parameters.")
    shrinkParams <- estimateShrinkageParameters(res)    
  } else {
    shrinkParams <- NULL
  }

  message("Computing p values.")
  p.value <- foreach(R = isplitVector(res, chunkSize=chunkSize), 
    .export=export) %dopar% {
    lapply(R, function(r) {
      p <- computePValue(r$b, r$w, r$fit, r$fit0, 
        shrinkParams = shrinkParams)
      e <- NULL
      if (effect == TRUE)
        e <- computeEffect(r$fit, model, model0)
      list(
        p.value = p,
        effect = e,
        fW = r$fit$warnings,
        f0W = r$fit0$warnings)
    })
  }
  p.value <- unlist(p.value, recursive=FALSE)
  names(p.value) <- rownames(beta)

  return(p.value)
}

