fit.set <- function(obj, ...){

  regr.vinfl <- estim.regrvinfl(attr(obj, "likdat")$anclik.dat)
  
  out <- list()
  attr(out, "country") <- attr(obj, "country")
  attr(out, "region") <- attr(obj, "region")


  print("rtrend none")
  out$rtrend.none <- fitmod(obj, eppmod="rtrend", v.infl=0, ...)
  print("rtrend regression")
  out$rtrend.regr <- do.call(fitmod, list(obj=obj, eppmod="rtrend", v.infl=regr.vinfl, ...))
  print("rtrend mle")
  out$rtrend.mle <- fitmod(obj, eppmod="rtrend", vinfl.method='mle', ...)
  print("rtrend unbiased")
  out$rtrend.unbiased <- fitmod(obj, eppmod="rtrend", vinfl.method='unbiased', ...)
  print("rtrend fit")
  out$rtrend.fit <- fitmod(obj, eppmod="rtrend", vinfl.prior.rate=1/0.015, ...)
  print("rtrend fit 1")
  out$rtrend.fit1 <- fitmod(obj, eppmod="rtrend", vinfl.prior.rate=1/0.1, ...)
  print("rtrend fit 2")
  out$rtrend.fit2 <- fitmod(obj, eppmod="rtrend", vinfl.prior.rate=1, ...)
  

  print("rspline none")
  out$rspline.none <- fitmod(obj, eppmod="rspline", equil.rprior=TRUE, v.infl=0, ...)
  print("rspline regression")
  out$rspline.regr <- do.call(fitmod, list(obj=obj, eppmod="rspline", equil.rprior=TRUE, v.infl=regr.vinfl, ...))
  print("rspline mle")
  out$rspline.mle <- fitmod(obj, eppmod="rspline", equil.rprior=TRUE, vinfl.method='mle', ...)
  print("rspline unbiased")
  out$rspline.unbiased <- fitmod(obj, eppmod="rspline", equil.rprior=TRUE, vinfl.method='unbiased', ...)
  print("rspline fit")
  out$rspline.fit <- fitmod(obj, eppmod="rspline", equil.rprior=TRUE, vinfl.prior.rate=1/0.015, ...)
  print("rspline fit 1")
  out$rspline.fit1 <- fitmod(obj, eppmod="rspline", equil.rprior=TRUE, vinfl.prior.rate=1/0.1, ...)
  print("rspline fit 2")
  out$rspline.fit2 <- fitmod(obj, eppmod="rspline", equil.rprior=TRUE, vinfl.prior.rate=1, ...)

  
  print("rsplinenoeq none")
  out$rsplinenoeq.none <- fitmod(obj, eppmod="rspline", v.infl=0, ...)
  print("rsplinenoeq regression")
  out$rsplinenoeq.regr <- do.call(fitmod, list(obj=obj, eppmod="rspline", v.infl=regr.vinfl, ...))
  print("rsplinenoeq mle")
  out$rsplinenoeq.mle <- fitmod(obj, eppmod="rspline", vinfl.method='mle', ...)
  print("rsplinenoeq unbiased")
  out$rsplinenoeq.unbiased <- fitmod(obj, eppmod="rspline", vinfl.method='unbiased', ...)
  print("rsplinenoeq fit")
  out$rsplinenoeq.fit <- fitmod(obj, eppmod="rspline", vinfl.prior.rate=1/0.015, ...)
  print("rsplinenoeq fit 1")
  out$rsplinenoeq.fit1 <- fitmod(obj, eppmod="rspline", vinfl.prior.rate=1/0.1, ...)
  print("rsplinenoeq fit 2")
  out$rsplinenoeq.fit2 <- fitmod(obj, eppmod="rspline", vinfl.prior.rate=1, ...)

  return(out)
}



####################################
####  Bao regression estimator  ####
####################################

estim.regrvinfl <- function(ancdat){
  dat <- data.frame(site=rep(names(ancdat$W.lst), sapply(ancdat$W.lst, length)),
                    year=factor(unlist(ancdat$anc.idx)),
                    Wst=unlist(ancdat$W.lst),
                    vst=unlist(ancdat$v.lst))
  if(length(levels(dat$site)) == 1)
    fit <- lm(Wst~year,data=dat)
  else
    fit <- lm(Wst~site+year,data=dat)
  v.infl <- max(mean(fit$residuals^2 - dat$vst), 0)
  return(v.infl)
}




## ###################################################################################
## ####  Function to randomly partition ANC data into training and test datasets  ####
## ###################################################################################

## version of sample does not interpret as 1:x if length(x)==1
sample2 <- function(x, size){
  x[sample.int(length(x), size)]
}

sample.anclikdat <- function(anclik.dat, minobs=1, proptest=0.1){
  ## minobs: minimum number of observations per site to retain
  ## proptest: proportion of ANC prevalence observations to withold for test data

  nsites <- length(anclik.dat$W.lst)
  nobs <- sapply(anclik.dat$W.lst, length)
  idx <- seq_len(sum(nobs))

  ## sample minobs observations per site to keep
  keep.idx <- unlist(tapply(idx, rep(1:nsites, nobs), function(x) sample2(x, min(minobs, length(x)))))

  ## sample indices for test data, omitting keep.idx to ensure retain minobs per site
  test.idx <- sample2(idx[-keep.idx], round(proptest*sum(nobs)))

  ## training dataset is complement of test dataset
  train.idx <- idx[-test.idx]

  ## convert to lists of indexes per site
  test.idx <- split(test.idx, cut(test.idx, c(0, cumsum(nobs))))
  test.idx <- mapply("-", test.idx, cumsum(nobs) - nobs, SIMPLIFY=FALSE)
  train.idx <- split(train.idx, cut(train.idx, c(0, cumsum(nobs))))
  train.idx <- mapply("-", train.idx, cumsum(nobs) - nobs, SIMPLIFY=FALSE)
  
  train.anclik.dat <- lapply(anclik.dat[1:4], function(lst) mapply("[", lst, train.idx, SIMPLIFY=FALSE))
  test.anclik.dat <- lapply(anclik.dat[1:4], function(lst) mapply("[", lst, test.idx, SIMPLIFY=FALSE))

  return(list(train.anclik.dat=train.anclik.dat,
              test.anclik.dat=test.anclik.dat))
}



## ##########################################################
## ####  Sample from posterior predictive distributions  ####
## ##########################################################

add.v.infl <- function(fit){
  if(!is.null(fit$fp$v.infl))  # if fixed v.infl supplied
    fit$v.infl <- rep(fit$fp$v.infl, length(fit$param))
  else if(!is.null(fit$fp$vinfl.prior.rate))  # v.infl estimated
    fit$v.infl <- sapply(fit$param, "[[", "v.infl")
  else { # MLE estimator
    qM.mat <- sweep(qnorm(fit$prev), 2, sapply(fit$param, "[[", "ancbias"), "+")
    fit$v.infl <- apply(qM.mat, 2, epp:::calc.v.infl, fit$likdat$anclik.dat, fit$fp)
  }
  return(fit)
}

calc.b.site <- function(fit){
  qM.mat <- sweep(qnorm(fit$prev), 2, sapply(fit$param, "[[", "ancbias"), "+")
  anclik.dat <- fit$likdat$anclik.dat

  ## function to add v.infl value to clinic level variance
  vinfl.likdat <- function(anclik.dat, v.infl){
    anclik.dat$v.lst <- lapply(anclik.dat$v.lst, "+", v.infl);
    anclik.dat$v.lst <- lapply(anclik.dat$v.lst, pmax, 0);
    return(anclik.dat)
  }

  b.site <- matrix(sapply(seq_len(ncol(qM.mat)), function(ii) anclik::sample.b.site(qM.mat[,ii], vinfl.likdat(anclik.dat, fit$v.infl[ii]))), ncol=ncol(qM.mat))
  return(b.site)
}

calc.pred.site <- function(fit, anc.preddat=fit$likdat$anclik.dat){
  qM.mat <- sweep(qnorm(fit$prev), 2, sapply(fit$param, "[[", "ancbias"), "+")
  pred.site <- lapply(seq(along=fit$param), function(ii) anclik::sample.pred.site(qM.mat[,ii], fit$b.site[,ii], anc.preddat, v.infl=fit$v.infl[ii]))
  return(pred.site)
}

pred.coverage <- function(anc.pred, anc.preddat){
  pred.quant <- apply(matrix(sapply(anc.pred, unlist), ncol=length(anc.pred)), 1,
                      quantile, c(0.025, 0.975))
  obs <- pnorm(unlist(anc.preddat$W.lst))
  return(mean(obs > pred.quant[1,] & obs < pred.quant[2,]))
}

calc.pred.quantile <- function(anc.pred, anc.preddat){
  pred.mat <- matrix(sapply(anc.pred, unlist), ncol=length(anc.pred))
  obs <- pnorm(unlist(anc.preddat$W.lst))
  pred.quant <- sapply(seq_along(obs), function(i) ecdf(pred.mat[i,])(obs[i]))
  pred.quant <- split(pred.quant, rep(names(anc.preddat$W.lst), sapply(anc.preddat$W.lst, length)))
  return(pred.quant)
}

add.sigma2 <- function(fit){
  fit$sigma2 <- anclik::sample.sigma2(fit$b.site)
  return(fit)
}


## #####################################################
## ####  Function to fit model and compile results  ####
## #####################################################


fitmod_outpred <- function(obj, test.anclik.dat, ...){

  fit <- fitmod(obj, ...)
  fit$test.anclikdat <- test.anclik.dat
  
  ##   ## posterior predictive distributions from resample
  fit <- simfit(fit, fit$fp$eppmod == "rspline")
  fit <- add.v.infl(fit)
  fit$b.site <- calc.b.site(fit)
  fit <- add.sigma2(fit)
  pred.site <- calc.pred.site(fit)
  fit$pred.quant <- calc.pred.quantile(pred.site, fit$likdat$anclik.dat)
  fit$fit.pred.cov <- pred.coverage(pred.site, fit$likdat$anclik.dat)
  fit$test.pred.site <- calc.pred.site(fit, fit$test.anclikdat)
  fit$test.pred.quant <- calc.pred.quantile(fit$test.pred.site, fit$test.anclikdat)
  fit$test.pred.cov <- pred.coverage(fit$test.pred.site, fit$test.anclikdat)

  fit$b.site <- NULL
  fit$rvec <- fit$rvec.spline <- NULL
  fit$prev <- fit$incid <- fit$popsize <- fit$pregprev <- NULL
  fit$param <- NULL
  fit$resample <- NULL
  fit$test.pred.site <- NULL

  return(fit)
}

sim_outpred <- function(obj, seed, minobs=1, proptest=0.1, ...){

  set.seed(seed)
  
  testdat <- sample.anclikdat(attr(obj, 'likdat')$anclik.dat, minobs, proptest)
  attr(obj, 'likdat')$anclik.dat <- testdat$train.anclik.dat
  test.anclik.dat <- testdat$test.anclik.dat

  regr.vinfl <- estim.regrvinfl(attr(obj, "likdat")$anclik.dat)

  out <- list()
  attr(out, "seed") <- seed
  attr(out, "country") <- attr(obj, "country")
  attr(out, "region") <- attr(obj, "region")

  print("rtrend none")
  out$rtrend.none <- fitmod_outpred(obj, test.anclik.dat, eppmod="rtrend", v.infl=0, ...)
  print("rtrend regression")
  out$rtrend.regr <- do.call(fitmod_outpred, list(obj=obj, test.anclik.dat=test.anclik.dat, eppmod="rtrend", v.infl=regr.vinfl, ...))
  print("rtrend mle")
  out$rtrend.mle <- fitmod_outpred(obj, test.anclik.dat, eppmod="rtrend", vinfl.method='mle', ...)
  print("rtrend unbiased")
  out$rtrend.unbiased <- fitmod_outpred(obj, test.anclik.dat, eppmod="rtrend", vinfl.method='unbiased', ...)
  print("rtrend fit")
  out$rtrend.fit <- fitmod_outpred(obj, test.anclik.dat, eppmod="rtrend", vinfl.prior.rate=1/0.015, ...)

  print("rspline none")
  out$rspline.none <- fitmod_outpred(obj, test.anclik.dat, eppmod="rspline", equil.rprior=TRUE, v.infl=0, ...)
  print("rspline regression")
  out$rspline.regr <- do.call(fitmod_outpred, list(obj=obj, test.anclik.dat=test.anclik.dat, eppmod="rspline", equil.rprior=TRUE, v.infl=regr.vinfl, ...))
  print("rspline mle")
  out$rspline.mle <- fitmod_outpred(obj, test.anclik.dat, eppmod="rspline", equil.rprior=TRUE, vinfl.method='mle', ...)
  print("rspline unbiased")
  out$rspline.unbiased <- fitmod_outpred(obj, test.anclik.dat, eppmod="rspline", equil.rprior=TRUE, vinfl.method='unbiased', ...)
  print("rspline fit")
  out$rspline.fit <- fitmod_outpred(obj, test.anclik.dat, eppmod="rspline", equil.rprior=TRUE, vinfl.prior.rate=1/0.015, ...)

  print("rsplinenoeq none")
  out$rsplinenoeq.none <- fitmod_outpred(obj, test.anclik.dat, eppmod="rspline", v.infl=0, ...)
  print("rsplinenoeq regression")
  out$rsplinenoeq.regr <- do.call(fitmod_outpred, list(obj=obj, test.anclik.dat=test.anclik.dat, eppmod="rspline", v.infl=regr.vinfl, ...))
  print("rsplinenoeq mle")
  out$rsplinenoeq.mle <- fitmod_outpred(obj, test.anclik.dat, eppmod="rspline", vinfl.method='mle', ...)
  print("rsplinenoeq unbiased")
  out$rsplinenoeq.unbiased <- fitmod_outpred(obj, test.anclik.dat, eppmod="rspline", vinfl.method='unbiased', ...)
  print("rsplinenoeq fit")
  out$rsplinenoeq.fit <- fitmod_outpred(obj, test.anclik.dat, eppmod="rspline", vinfl.prior.rate=1/0.015, ...)

  return(out)
}
