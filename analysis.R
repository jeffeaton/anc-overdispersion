setwd("~/Documents/Research/anc-overdispersion/paper/anc-overdispersion/")

devtools::install_github("jeffeaton/anclik/anclik@79904d4")  # tagged git commits used for analysis
devtools::install_github("jeffeaton/epp@8d17e41")

library(epp)
source("functions.R")


#########################
####  Load EPP info  ####
#########################

## Files publicly available from UNAIDS by request: http://apps.unaids.org/spectrum/

filepath <- paste0("~/Documents/Data/Spectrum files/2016 final/SSA/",
                   list("Botswana_ Final_15_04_ 2016 upd.PJNZ",
                        "Lesotho_ 2016_20May.PJNZ",
                        "Malawi_2016 4May old structure.PJNZ",
                        "TZ_Mainland_6May2016.PJNZ",
                        "Uganda  21st April 2016.PJNZ",
                        "SouthAfrica_27May2016.PJNZ",
                        "Zambia 2016_NATIONAL Final.PJNZ",
                        "Kenya-Central Apr 6.PJNZ",
                        "Kenya-Coast Apr 6.PJNZ",
                        "Kenya-Eastern Apr 6.PJNZ",
                        "Kenya-Nairobi Apr 6upd.PJNZ",
                        "Kenya-North Eastern Apr 6.PJNZ",
                        "Kenya-Nyanza Apr 6.PJNZ",
                        "Kenya-Rift Valley Apr 6.PJNZ",
                        "Kenya-Western Apr 6.PJNZ",
                        "Zimbabwe_Midlands_revised updrevTFR.PJNZ",
                        "Zimbabwe_Matabeleland South_revised updrevTFR.PJNZ",
                        "Zimbabwe_Matabeleland North_revised updrevTFR.PJNZ",
                        "Zimbabwe_Masvingo_Revised updrevTFR.PJNZ",
                        "Zimbabwe_Mashonaland West_revisedTFRrev.PJNZ",
                        "Zimbabwe_Mashonaland Central_Revised updrevTFR.PJNZ",
                        "Zimbabwe_Manicaland_revised updrevTFR.PJNZ",
                        "Zimbabwe_Harare_4May updrevTFR.PJNZ",
                        "Zimbabwe_Mashonaland East_revised updrevTFR.PJNZ",
                        "Zimbabwe_Bulawayo_4may2016 updrevTFR.PJNZ"))
names(filepath) <- sub(".*/([^_ ]+).*.", "\\1", filepath)

inputs <- lapply(filepath, prepare_epp_fit)
inputs <- do.call(c, inputs)

inputs[23:30] <- mapply("attr<-", inputs[23:30], "region", c("Central", "Coast", "Eastern", "Nairobi", "North Eastern", "Nyanza", "Rift Valley", "Western"))
inputs[31:40] <- mapply("attr<-", inputs[31:40], "region", c("Midlands", "Matabeleland South", "Matabeleland North", "Masvingo", "Mashonaland West", "Mashonaland Central", "Manicaland", "Harare", "Mashonaland East", "Bulawayo"))

names(inputs) <- mapply(paste, lapply(inputs, attr, "country"), lapply(inputs, attr, "region"))


#############################
####  Configure cluster  ####
#############################

## Note: this is specific to Imperial College Department of Infectious Disease Epidemiology computing cluster.

workdir <- "/Volumes/jeff/anc-overdispersion/"
## system(paste0("mkdir ", workdir))
didewin::didewin_config_global(credentials="jwe08", cluster="fi--didemrchnb", workdir=workdir)

ctx <- context::context_save(workdir, packages="epp",
                             package_sources=context::package_sources(github=c("jeffeaton/anclik/anclik@79904d4/",
                                                                               "jeffeaton/epp@8d17e41")),
                             sources="functions.R")



mrcq <- didewin::queue_didewin(ctx, initialise=FALSE)

mrcq$sync_files()


###############################
####  Derivation of prior  ####
###############################

all.regrvinfl <- sapply(lapply(lapply(inputs, attr, "likdat"), "[[", "anclik.dat"), estim.regrvinfl)

mean(all.regrvinfl)
median(all.regrvinfl)
range(all.regrvinfl)


plot(density(all.regrvinfl, from=0))


######################
####  Fit models  ####
######################

library(queuer)

## Fit model variants to full data from 40 regions
qlapply(inputs, fit.set, mrcq, name="allfits")


## Out-of-sample prediction -- simulate 50 test/train datasets for each of 40 datasets
## (submitted in batches of 500 jobs at a time)
set.seed(461486902)
outpred.pars <- mapply(list, obj=rep(inputs, 50), seed=sample(1e9, 50*length(inputs)), SIMPLIFY=FALSE)

enqueue_bulk(mrcq, outpred.pars[1:500], sim_outpred, do.call=TRUE, name="outpred1")
enqueue_bulk(mrcq, outpred.pars[500+1:500], sim_outpred, do.call=TRUE, name="outpred2")
enqueue_bulk(mrcq, outpred.pars[1000+1:500], sim_outpred, do.call=TRUE, name="outpred3")
enqueue_bulk(mrcq, outpred.pars[1500+1:500], sim_outpred, do.call=TRUE, name="outpred4")


## Sensitivity analysis: fit EPP model to ANC SS data only, omitting household survey prevalence

inputs_nohhs <- lapply(inputs, function(obj){attr(obj, "likdat")$hhslik.dat <- attr(obj, "likdat")$hhslik.dat[NULL,]; obj})

qlapply(inputs_nohhs, fitmod, mrcq, name="nohhs_rtrend",
        eppmod="rtrend", vinfl.prior.rate=1/0.015, ancbias.pr.mean=0.15, ancbias.pr.sd=0.15)
qlapply(inputs_nohhs, fitmod, mrcq, name="nohhs_rtrend",
        eppmod="rspline", equil.rprior=TRUE, vinfl.prior.rate=1/0.015, ancbias.pr.mean=0.15, ancbias.pr.sd=0.15)
qlapply(inputs_nohhs, fitmod, mrcq, name="nohhs_rsplinenoeq",
        eppmod="rspline", equil.rprior=FALSE, vinfl.prior.rate=1/0.015, ancbias.pr.mean=0.15, ancbias.pr.sd=0.15)



######################################
####  Load model fitting results  ####
######################################

allfit <- mrcq$task_bundle_get("allfits")$results()

outpred <- c(mrcq$task_bundle_get("outpred1")$results(),
             mrcq$task_bundle_get("outpred2")$results(),
             mrcq$task_bundle_get("outpred3")$results(),
             mrcq$task_bundle_get("outpred4")$results())

outpred <- tapply(outpred, paste(sapply(outpred, attr, "country"), sapply(outpred, attr, "region")), list)

nohhs <- list(rtrend = mrcq$task_bundle_get("nohhs_rtrend")$results(),
              rspline = mrcq$task_bundle_get("nohhs_rspline")$results(),
              rsplinenoeq = mrcq$task_bundle_get("nohhs_rsplinenoeq")$results())


## Simulate model outputs 
sim_outputs <- function(fit){

  fit <- simfit(fit, fit$fp$eppmod == "rspline")
  fit <- add.v.infl(fit)
  fit$b.site <- calc.b.site(fit)
  fit <- add.sigma2(fit)
  pred.site <- calc.pred.site(fit)
  fit$pred.quant <- calc.pred.quantile(pred.site, fit$likdat$anclik.dat)
  fit$fit.pred.cov <- pred.coverage(pred.site, fit$likdat$anclik.dat)

  fit$param <- NULL

  return(fit)
}

library(parallel)
options(mc.cores=6)

fitsim <- mclapply(allfit, function(set){ print(paste(attr(set, "country"), attr(set, "region"))); lapply(set, sim_outputs)})

nohhs <- mclapply(nohhs, lapply, sim_outputs)


###################
####  Table 1  ####
###################

## Estimates of vinfl
all.vinfl.est <- lapply(fitsim, lapply, "[[", "v.infl")
all.vinfl.mean <- data.frame(t(sapply(all.vinfl.est, sapply, mean)))

## sampling error
vst.mean <- sapply(lapply(lapply(lapply(lapply(inputs, "attr", "likdat"), "[[", "anclik.dat"), "[[", "v.lst"), unlist), mean)

tab1 <- function(x) sprintf("%.3f %.3f (%.3f, %.3f) (%.3f, %.3f)", median(x), mean(x), quantile(x, 0.25), quantile(x, 0.75), min(x), max(x))

rbind(tab1(all.vinfl.mean$rtrend.regr),
      tab1(all.vinfl.mean$rtrend.unbia),
      tab1(all.vinfl.mean$rtrend.fit),
      tab1(all.vinfl.mean$rspline.regr),
      tab1(all.vinfl.mean$rspline.unbia),
      tab1(all.vinfl.mean$rspline.fit),
      tab1(all.vinfl.mean$rsplinenoeq.regr),
      tab1(all.vinfl.mean$rsplinenoeq.unbia),
      tab1(all.vinfl.mean$rsplinenoeq.fit),
      tab1(vst.mean),
      tab1((all.vinfl.mean$rspline.fit+vst.mean) / vst.mean))



###################
####  Table 2  ####
###################


## Fit to national survey data.
## Calculate log posterior predictive density (LPPD) for household survey data.
## Calculation following Gelman, Hwang, Vehtari 2013, equation 5.

calc.lppd <- function(fit){
  hhsdat <- fit$likdat$hhslik.dat
  hhsll <- matrix(dnorm(hhsdat$W.hhs, qnorm(fit$prev[hhsdat$idx,]), hhsdat$sd.W.hhs, TRUE), length(hhsdat$W.hhs))
  sum(apply(hhsll, 1, pomp::logmeanexp))
}

hhs.lppd <- t(sapply(fitsim, sapply, calc.lppd))

hhs.change.lppd <- cbind(hhs.lppd[,paste0("rtrend.", c("regr", "mle", "unbiased", "fit"))] - hhs.lppd[,"rtrend.none"],
                         hhs.lppd[,paste0("rspline.", c("regr", "mle", "unbiased", "fit", "fit1", "fit2"))] - hhs.lppd[,"rspline.none"],
                         hhs.lppd[,paste0("rsplinenoeq.", c("regr", "mle", "unbiased", "fit", "fit1", "fit2"))] - hhs.lppd[,"rsplinenoeq.none"])


## Mean and IQR of number of IMIS iterations       
imis.iterations <- lapply(outpred, sapply, sapply, function(x) nrow(x$stat))


insample.fit <- rowMeans(sapply(fitsim, sapply, "[[", "fit.pred.cov"))

outpred.cov <- rowMeans(sapply(lapply(outpred, sapply, sapply, "[[", "test.pred.cov"), rowMeans)) # out-of-sample prediction

colMeans(hhs.lppd[,t(outer(c("rtrend.", "rspline.", "rsplinenoeq."), c("none", "regr", "mle", "unbiased", "fit"), paste0))]) # HHS LPPD
colMeans(hhs.change.lppd)


## Table 2
out <- t(outer(c("rtrend.", "rspline.", "rsplinenoeq."), c("none", "regr", "unbiased", "fit"), paste0))

rbind(round(100*insample.fit, 1)[out],
      round(colMeans(hhs.lppd)[out], 2),
      round(100*outpred.cov[out], 1),
      round(colMeans(hhs.change.lppd)[out], 2),
      round(rowMeans(sapply(imis.iterations, apply, 1, median)), 1)[out], # number of IMIS iterations
      round(rowMeans(sapply(imis.iterations, apply, 1, function(x) diff(quantile(x, c(0.25, 0.75))))), 1)[out]) # IQR of IMIS iterations
      

##########################################
####  Table 3: Effect on uncertainty  ####
##########################################


## Coefficient of variation

coefvar <- function(x, idx = c(21, 26, 31, 36)){ ## 1990, 1995, 2000, 2005
  apply(x[idx,,drop=FALSE], 1, sd)/rowMeans(x[idx,,drop=FALSE])
}

rtrend.none.prevcv <- sapply(lapply(lapply(fitsim, "[[", "rtrend.none"), "[[", "prev"), coefvar)
rtrend.fit.prevcv <- sapply(lapply(lapply(fitsim, "[[", "rtrend.fit"), "[[", "prev"), coefvar)

rspline.none.prevcv <- sapply(lapply(lapply(fitsim, "[[", "rspline.none"), "[[", "prev"), coefvar)
rspline.fit.prevcv <- sapply(lapply(lapply(fitsim, "[[", "rspline.fit"), "[[", "prev"), coefvar)

rsplinenoeq.none.prevcv <- sapply(lapply(lapply(fitsim, "[[", "rsplinenoeq.none"), "[[", "prev"), coefvar)
rsplinenoeq.fit.prevcv <- sapply(lapply(lapply(fitsim, "[[", "rsplinenoeq.fit"), "[[", "prev"), coefvar)


rtrend.none.incidcv <- sapply(lapply(lapply(fitsim, "[[", "rtrend.none"), "[[", "incid"), coefvar)
rtrend.fit.incidcv <- sapply(lapply(lapply(fitsim, "[[", "rtrend.fit"), "[[", "incid"), coefvar)

rspline.none.incidcv <- sapply(lapply(lapply(fitsim, "[[", "rspline.none"), "[[", "incid"), coefvar)
rspline.fit.incidcv <- sapply(lapply(lapply(fitsim, "[[", "rspline.fit"), "[[", "incid"), coefvar)

rsplinenoeq.none.incidcv <- sapply(lapply(lapply(fitsim, "[[", "rsplinenoeq.none"), "[[", "incid"), coefvar)
rsplinenoeq.fit.incidcv <- sapply(lapply(lapply(fitsim, "[[", "rsplinenoeq.fit"), "[[", "incid"), coefvar)


rtrend.none.cv <- rbind(rtrend.none.prevcv, rtrend.none.incidcv)
rtrend.fit.cv <- rbind(rtrend.fit.prevcv, rtrend.fit.incidcv)

rspline.none.cv <- rbind(rspline.none.prevcv, rspline.none.incidcv)
rspline.fit.cv <- rbind(rspline.fit.prevcv, rspline.fit.incidcv)

rsplinenoeq.none.cv <- rbind(rsplinenoeq.none.prevcv, rsplinenoeq.none.incidcv)
rsplinenoeq.fit.cv <- rbind(rsplinenoeq.fit.prevcv, rsplinenoeq.fit.incidcv)


mediqrfmt <- function(x) sprintf("%.2f [%.2f-%.2f]", median(x), quantile(x, 0.25), quantile(x, 0.75))

rbind(apply(rtrend.none.cv, 1, mediqrfmt),
      apply(rtrend.fit.cv, 1, mediqrfmt),
      apply(rtrend.fit.cv/rtrend.none.cv, 1, mediqrfmt),
      ##
      apply(rspline.none.cv, 1, mediqrfmt),
      apply(rspline.fit.cv, 1, mediqrfmt),
      apply(rspline.fit.cv/rspline.none.cv, 1, mediqrfmt),
      ##
      apply(rsplinenoeq.none.cv, 1, mediqrfmt),
      apply(rsplinenoeq.fit.cv, 1, mediqrfmt),
      apply(rsplinenoeq.fit.cv/rsplinenoeq.none.cv, 1, mediqrfmt))



##########################################################
####  Supplementary Appendix: Plot page for each fit  ####
##########################################################

## prevalence, incidence, r(t)
## rtrend, rspline w/ equil, rspline no equil
## current, regr, mle, estimated1

## QQ plot
## posterior distribution of vinfl


## Create plot functions
library(RColorBrewer)
dark2 <- brewer.pal(5, "Dark2")

cred.region <- function(x, y, ...)
  polygon(c(x, rev(x)), c(y[1,], rev(y[2,])), border=NA, ...)

transp <- function(col, alpha=0.5)
  return(apply(col2rgb(col), 2, function(c) rgb(c[1]/255, c[2]/255, c[3]/255, alpha)))

plot.prev <- function(fit, ..., ylim=NULL, xlim=c(1980, 2015), col="blue"){
  if(is.null(ylim))
    ylim <- c(0, 1.1*max(apply(fit$prev, 1, quantile, 0.975)))
  plot(1970:2015, rowMeans(fit$prev), type="n", ylim=ylim, xlim=xlim, ylab="prevalence", xlab="", yaxt="n", xaxt="n")
  axis(1, labels=TRUE)
  axis(2, labels=TRUE)
  dots <- list(...)
  for(ii in seq_along(dots))
    cred.region(1970:2015, apply(dots[[ii]]$prev, 1, quantile, c(0.025, 0.975)), col=transp(col[1+ii], 0.3))
  cred.region(1970:2015, apply(fit$prev, 1, quantile, c(0.025, 0.975)), col=transp(col[1], 0.3))
  for(ii in seq_along(dots))
    lines(1970:2015, rowMeans(dots[[ii]]$prev), col=col[1+ii])
  lines(1970:2015, rowMeans(fit$prev), col=col[1])
  ##
  points(fit$likdat$hhslik.dat$year, fit$likdat$hhslik.dat$prev, pch=20)
  segments(fit$likdat$hhslik.dat$year,
           y0=pnorm(fit$likdat$hhslik.dat$W.hhs - qnorm(0.975)*fit$likdat$hhslik.dat$sd.W.hhs),
           y1=pnorm(fit$likdat$hhslik.dat$W.hhs + qnorm(0.975)*fit$likdat$hhslik.dat$sd.W.hhs))
}

plot.incid <- function(fit, ..., ylim=NULL, xlim=c(1980, 2015), col="blue"){
  if(is.null(ylim))
    ylim <- c(0, 1.1*max(apply(fit$incid, 1, quantile, 0.975)))
  plot(1970:2015, rowMeans(fit$incid), type="n", ylim=ylim, xlim=xlim, ylab="incidence rate", xlab="", yaxt="n", xaxt="n")
  axis(1, labels=TRUE)
  axis(2, labels=TRUE)
  dots <- list(...)
  for(ii in seq_along(dots))
    cred.region(1970:2015, apply(dots[[ii]]$incid, 1, quantile, c(0.025, 0.975)), col=transp(col[1+ii], 0.3))
  cred.region(1970:2015, apply(fit$incid, 1, quantile, c(0.025, 0.975)), col=transp(col[1], 0.3))
  for(ii in seq_along(dots))
    lines(1970:2015, rowMeans(dots[[ii]]$incid), col=col[1+ii])
  lines(1970:2015, rowMeans(fit$incid), col=col)
}

plot.rvec <- function(fit, ..., ylim=NULL, xlim=c(1980, 2015), col="blue"){
  dots <- list(...)
  fit$rvec <- sapply(seq_len(ncol(fit$rvec)), function(i) replace(fit$rvec[,i], fit$fp$proj.steps < fit$param[[i]]$tsEpidemicStart, NA))
  for(ii in seq_along(dots))
    dots[[ii]]$rvec <- sapply(seq_len(ncol(dots[[ii]]$rvec)), function(i) replace(dots[[ii]]$rvec[,i], dots[[ii]]$fp$proj.steps < dots[[ii]]$param[[i]]$tsEpidemicStart, NA))
  if(is.null(ylim))
    ylim <- c(0, quantile(apply(fit$rvec, 1, quantile, 0.975, na.rm=TRUE), 0.975, na.rm=TRUE))
  plot(fit$fp$proj.steps, rowMeans(fit$rvec, na.rm=TRUE), type="n", ylim=ylim, xlim=xlim, ylab="r(t)", yaxt="n", xlab="")
  axis(1, labels=TRUE)
  axis(2, labels=TRUE)
  for(ii in seq_along(dots))
    cred.region(dots[[ii]]$fp$proj.steps, apply(dots[[ii]]$rvec, 1, quantile, c(0.025, 0.975), na.rm=TRUE), col=transp(col[1+ii], 0.3))
  cred.region(fit$fp$proj.steps, apply(fit$rvec, 1, quantile, c(0.025, 0.975), na.rm=TRUE), col=transp(col[1], 0.3))
  for(ii in seq_along(dots))
    lines(dots[[ii]]$fp$proj.steps, rowMeans(dots[[ii]]$rvec, na.rm=TRUE), col=col[1+ii])
  lines(fit$fp$proj.steps, rowMeans(fit$rvec, na.rm=TRUE), col=col[1])
  return(invisible())
}

plot.qqplot <- function(fits, col="blue"){
  predquant <- sapply(fits, function(fit) sort(unlist(fit$pred.quant)))
  matplot(seq(0, 1, length.out=nrow(predquant)), predquant,
          pch=20, cex=0.5, col=col,
          xlab="Theoretical quantiles",
          ylab="Observed quantiles")
  abline(a=0, b=1)
  return(invisible())
}

plot.vinflpost <- function(fits, col=col){
  plot(NA, NA, type="n", xlim=c(0, 0.06), ylim=c(0, 400), xlab=expression(sigma["infl"]^2), ylab="density")
  for(i in seq_along(fits))
    lines(density(fits[[i]]$v.infl, from=0), col=col[i])
  return(invisible())
}

plot.obj <- function(fit, name){
  par(oma=c(0, 5, 1, 5), cex=1)
  layout(rbind(0, rep(1:3, each=2), c(0, rep(4:5, each=2), 0),
               0, rep(6:8, each=2), c(0, rep(9:10, each=2), 0),
               0, rep(11:13, each=2), c(0, rep(14:15, each=2), 0), 0),
         h=c(0.1, 1, 0.8, 0.1, 1, 0.8, 0.1, 1, 0.8, 0.1))
  par(mgp=c(1.8, 0.5, 0), tcl=-0.25, mar=c(2.0, 3.3, 1.0, 0.5), cex=1.0)
  plot.prev(fit$rtrend.none, fit$rtrend.regr, fit$rtrend.unbiased, fit$rtrend.fit, col=dark2)
  mtext("r-trend", 3, 1, adj=-0.5, font=2, cex=1.2)
  mtext(name, 3, -1.5, TRUE, cex=1.5, font=2)
  plot.incid(fit$rtrend.none, fit$rtrend.regr, fit$rtrend.unbiased, fit$rtrend.fit, col=dark2)
  plot.rvec(fit$rtrend.none, fit$rtrend.regr, fit$rtrend.unbiased, fit$rtrend.fit, col=dark2)
  par(mar=c(3, 3, 0.5, 5))
  plot.qqplot(list(fit$rtrend.none, fit$rtrend.regr, fit$rtrend.unbiased, fit$rtrend.fit), col=dark2)
  plot.vinflpost(list(fit$rtrend.regr, fit$rtrend.unbiased, fit$rtrend.fit), dark2[-1])
  legend("right", c("none", "regr.", "unbias.", "estim."), col=dark2, lwd=2, inset=-1, xpd=NA)
  ##
  par(mgp=c(1.8, 0.5, 0), tcl=-0.25, mar=c(2.0, 3.3, 1.0, 0.5), cex=1.0)
  plot.prev(fit$rspline.none, fit$rspline.regr, fit$rspline.unbiased, fit$rspline.fit, col=dark2)
  mtext("r-spline", 3, 1, adj=-0.5, font=2, cex=1.2)
  mtext(name, 3, -1.5, TRUE, cex=1.5, font=2)
  plot.incid(fit$rspline.none, fit$rspline.regr, fit$rspline.unbiased, fit$rspline.fit, col=dark2)
  plot.rvec(fit$rspline.none, fit$rspline.regr, fit$rspline.unbiased, fit$rspline.fit, col=dark2)
  par(mar=c(3, 3, 0.5, 5))
  plot.qqplot(list(fit$rspline.none, fit$rspline.regr, fit$rspline.unbiased, fit$rspline.fit), col=dark2)
  plot.vinflpost(list(fit$rspline.regr, fit$rspline.unbiased, fit$rspline.fit), dark2[-1])
  legend("right", c("none", "regr.", "unbias.", "estim."), col=dark2, lwd=2, inset=-1, xpd=NA)
  ##
  par(mgp=c(1.8, 0.5, 0), tcl=-0.25, mar=c(2.0, 3.3, 1.0, 0.5), cex=1.0)
  plot.prev(fit$rsplinenoeq.none, fit$rsplinenoeq.regr, fit$rsplinenoeq.unbiased, fit$rsplinenoeq.fit, col=dark2)
  mtext("r-spline, no equil. prior", 3, 1, adj=-0.5, font=2, cex=1.2)
  mtext(name, 3, -1.5, TRUE, cex=1.5, font=2)
  plot.incid(fit$rsplinenoeq.none, fit$rsplinenoeq.regr, fit$rsplinenoeq.unbiased, fit$rsplinenoeq.fit, col=dark2)
  plot.rvec(fit$rsplinenoeq.none, fit$rsplinenoeq.regr, fit$rsplinenoeq.unbiased, fit$rsplinenoeq.fit, col=dark2)
  par(mar=c(3, 3, 0.5, 5))
  plot.qqplot(list(fit$rsplinenoeq.none, fit$rsplinenoeq.regr, fit$rsplinenoeq.unbiased, fit$rsplinenoeq.fit), col=dark2)
  plot.vinflpost(list(fit$rsplinenoeq.regr, fit$rsplinenoeq.unbiased, fit$rsplinenoeq.fit), dark2[-1])
  legend("right", c("none", "regr.", "unbias.", "estim."), col=dark2, lwd=2, inset=-1, xpd=NA)
  return(invisible())
}

pdf("figures/anc-overdispersion-summary.pdf", w=8.27, h=11.69, pointsize=8)
mapply(plot.obj, fitsim, names(fitsim))
dev.off()


      
####################
####  Figure 1  ####
####################

## Example of HIV prevalence and incidence forBotswana Urban and Zimbabwe Manicaland

pdf("figures/figure1.pdf", h=4.54, w=6.81, pointsize=9)
##
par(oma=c(0, 0, 1, 0))
layout(rbind(1:3, 4:6, 7), h=c(1, 1, 0))
par(mgp=c(1.8, 0.5, 0), tcl=-0.25, mar=c(3.0, 3.3, 1.0, 0.5), cex=1.0)
fit <- fitsim[["Botswana Urban"]]
plot.vinflpost(list(fit$rspline.regr, fit$rspline.unbiased, fit$rspline.fit), dark2[-1])
legend("topright", c("regr.", "unbias.", "estim."), col=dark2[-1], lwd=1.5, xpd=NA, cex=0.8)
plot.prev(fit$rspline.none, fit$rspline.fit, col=dark2[c(1,4)])
legend("topleft", c("none", "estim."), col=dark2[c(1,4)], lwd=1.5, xpd=NA, cex=0.8)
mtext("Botswana Urban", 3, 0.5, cex=1.3, font=2)
plot.incid(fit$rspline.none, fit$rspline.fit, col=dark2[c(1,4)])
legend("topright", c("none", "estim."), col=dark2[c(1,4)], lwd=1.5, xpd=NA, cex=0.8)
##
par(mgp=c(1.8, 0.5, 0), tcl=-0.25, mar=c(3.0, 3.3, 1.0, 0.5), cex=1.0)
fit <- fitsim[["Zimbabwe Manicaland"]]
plot.vinflpost(list(fit$rspline.regr, fit$rspline.unbiased, fit$rspline.fit), dark2[-1])
legend("topright", c("regr.", "unbias.", "estim."), col=dark2[-1], lwd=1.5, xpd=NA, cex=0.8)
plot.prev(fit$rspline.none, fit$rspline.fit, col=dark2[c(1,4)])
legend("topleft", c("none", "estim."), col=dark2[c(1,4)], lwd=1.5, xpd=NA, cex=0.8)
mtext("Zimbabwe Manicaland", 3, 0.5, cex=1.3, font=2)
plot.incid(fit$rspline.none, fit$rspline.fit, col=dark2[c(1,4)])
legend("topright", c("none", "estim."), col=dark2[c(1,4)], lwd=1.5, xpd=NA, cex=0.8)
##
dev.off()


##########################################
####  Figure 2: sensitivity to prior  ####
##########################################

## effect of prior on v.infl

pdf("figures/figure2.pdf", h=3.0, w=3.23, pointsize=9)
par(mgp=c(1.8, 0.5, 0), tcl=-0.25, mar=c(3, 3.5, 0.5, 0.5), cex=1)
plot(NA, NA, type="n", ylim=c(0, 0.07), xlim=c(0, 0.07),
     xlab=expression("posterior mean"~sigma["infl"]^2~", "~v[0]==0.015^{-1}),
     ylab=expression("posterior mean"~sigma["infl"]^2~", diffuse prior"), main="")
abline(a=0, b=1, col="grey80")
points(all.vinfl.mean$rspline.fit, all.vinfl.mean$rspline.fit2, pch=2, cex=0.8)
points(all.vinfl.mean$rspline.fit, all.vinfl.mean$rspline.fit1, pch=1, cex=0.8)
legend("topleft", legend=c(expression(v[0]==0.1^{-1}), expression(v[0]==1)), pch=1:2)
dev.off()


####################################
####  Correlation of estimates  ####
####################################

with(all.vinfl.mean, cor(rspline.unbiased, rspline.fit))


#################################
####  Determinants of vinfl  ####
#################################

dat <- data.frame(country = sapply(allfit, attr, "country"),
                  region = sapply(allfit, attr, "region"),
                  nsite = sapply(allfit, function(x) length(x$rspline.fit$likdat$anclik.dat$W.lst)),
                  medobs = sapply(allfit, function(x) median(sapply(x$rspline.fit$likdat$anclik.dat$W.lst, length))),
                  vinfl = all.vinfl.mean$rspline.fit,
                  logvinfl = log(all.vinfl.mean$rspline.fit))

dat <- merge(dat, data.frame(country=c("Botswana", "Lesotho", "Malawi", "United Republic of Tanzania", "Uganda", "South Africa", "Zambia", "Kenya", "Zimbabwe"),
                             code=c("BW", "LS", "MW", "TZ", "UG", "ZA", "ZM", "KE", "ZW")))
dat <- merge(dat, aggregate(cbind(nstrat=code) ~ country, dat, length))

summary(lm(logvinfl ~ country, dat))
anova(lm(logvinfl ~ country, dat))

summary(lm(logvinfl ~ nstrat, dat))
summary(lm(logvinfl ~ nsite, dat))
summary(lm(logvinfl ~ medobs, dat))

t.test(subset(dat, region == "Urban")$logvinfl, subset(dat, region == "Rural")$logvinfl, paired=TRUE)


############################################################
####  Sensitivity analysis: fit model without HHS data  ####
############################################################

pdf("figures/vinfl-no-hhs.pdf", h=3.0, w=3.23, pointsize=9)
par(mgp=c(1.8, 0.5, 0), tcl=-0.25, mar=c(3, 3.5, 0.5, 0.5), cex=1)
plot(NA, NA, type="n", ylim=c(0, 0.07), xlim=c(0, 0.07),
     xlab=expression("posterior mean"~sigma["infl"]^2~", incl. HHS prev"),
     ylab=expression("posterior mean"~sigma["infl"]^2~", ANC SS only"), main="")
abline(a=0, b=1, col="grey80")
points(sapply(lapply(lapply(fitsim, "[[", "rspline.fit"), "[[", "v.infl"), mean),
       sapply(lapply(nohhs$rspline, "[[", "v.infl"), mean), pch=20)
dev.off()
