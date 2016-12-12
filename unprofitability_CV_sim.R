## EI: CV Super Learner
library(SuperLearner)
library(parallel)
library(doMC)
options(mc.cores=80)

# full simulated data set
load("/data/markscan_authorized_users/bergquist/EI/df_sim_full.Rdata")
data <- df_sim
# if downloading data from github, will be in csv form

# check the data
summary(data)

## took this out because performed in SAS
# check for duplicate rows 
# dupsdf <- data[duplicated(data),]
# dim(dupsdf)
# rm(dupsdf)

# min and max for all variables besides unprofits
max(data[,-1])
min(data[,-1])

# max and min for unprofits
max(data[,1])
min(data[,1])

set.seed(27)
# get rid of variables with no obs
newdat <-data[, colSums(data != 0) > 0] 

dim(newdat)
cSums<-colSums(newdat)

# "screener" for taking out therapeutic class variables
tgrp.fun <- function(X, ...){
  whichvars <- c(rep.int(TRUE, ncol(X)))
  names(whichvars) <- colnames(X)
  tclsvars <- grep("tcls", names(X), value=T)
  whichvars[tclsvars] <- FALSE 
  whichvars <- unname(whichvars)
  return(whichvars) 
}


# lasso screener that always retains classes for HIV and MS drugs
var.index <- c(which(colnames(newdat)=="tcls14"), which(colnames(newdat)=="tcls251"))

screen.glmnet10 <- function(Y, X, family, alpha = 1, minscreen = 2, nfolds = 10, nlambda = 100,fixed.var.index=var.index,...) {
  # .SL.require('glmnet')
  if(!is.matrix(X)) {
    X <- model.matrix(~ -1 + ., X)
  }
  fitCV <- glmnet::cv.glmnet(x = X, y = Y, lambda = NULL, type.measure = 'deviance', 
                             nfolds = nfolds, family = family$family, alpha = alpha, 
                             nlambda = nlambda, pmax=10, parallel=T)
  whichVariable <- (as.numeric(coef(fitCV$glmnet.fit, s = fitCV$lambda.min))[-1] != 0)
  # the [-1] removes the intercept; taking the coefs from the fit w/ lambda that gives minimum cvm
  if (sum(whichVariable) < minscreen) {
    warning("fewer than minscreen variables passed the glmnet screen, 
            increased lambda to allow minscreen variables")
    sumCoef <- apply(as.matrix(fitCV$glmnet.fit$beta), 2, function(x) sum((x != 0)))
    newCut <- which.max(sumCoef >= minscreen) 
    whichVariable <- (as.matrix(fitCV$glmnet.fit$beta)[, newCut] != 0)
  }
  whichVariable[c(var.index)] <- TRUE
  return(whichVariable)
}

# lasso
SL.glmnet1 <- function(Y, X, newX, family, obsWeights, id, alpha = 1, nfolds = 10, nlambda = 100, useMin = TRUE, ...) {
  #.SL.require('glmnet')
  # X must be a matrix, should we use model.matrix or as.matrix
  if(!is.matrix(X)) {
    X <- model.matrix(~ -1 + ., X)
    newX <- model.matrix(~ -1 + ., newX)
  }
  # now use CV to find lambda
  fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights, lambda = NULL, type.measure = 'deviance', nfolds = nfolds, family = family$family, alpha = alpha, nlambda = nlambda, parallel=T)
  # two options for lambda, fitCV$lambda.min and fitCV$lambda.1se
  pred <- predict(fitCV$glmnet.fit, newx = newX, s = ifelse(useMin, fitCV$lambda.min, fitCV$lambda.1se), type = 'response')
  fit <- list(object = fitCV, useMin = useMin)
  class(fit) <- 'SL.glmnet'
  out <- list(pred = pred, fit = fit)
  return(out)
}

# ridge
SL.glmnet0 <- function(Y, X, newX, family, obsWeights, id, alpha = 0, nfolds = 10, nlambda = 100, useMin = TRUE, ...) {
  #.SL.require('glmnet')
  # X must be a matrix, should we use model.matrix or as.matrix
  if(!is.matrix(X)) {
    X <- model.matrix(~ -1 + ., X)
    newX <- model.matrix(~ -1 + ., newX)
  }
  # now use CV to find lambda
  fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights, lambda = NULL, type.measure = 'deviance', nfolds = nfolds, family = family$family, alpha = alpha, nlambda = nlambda, parallel=T)
  # two options for lambda, fitCV$lambda.min and fitCV$lambda.1se
  pred <- predict(fitCV$glmnet.fit, newx = newX, s = ifelse(useMin, fitCV$lambda.min, fitCV$lambda.1se), type = 'response')
  fit <- list(object = fitCV, useMin = useMin)
  class(fit) <- 'SL.glmnet'
  out <- list(pred = pred, fit = fit)
  return(out)
}

alg_screen <- list(c("SL.nnet", "All", "tgrp.fun", "screen.glmnet10"), 
c("SL.glmnet1", "All", "tgrp.fun", "screen.glmnet10"),
c("SL.glm", "All", "tgrp.fun", "screen.glmnet10"),
c("SL.glmnet0", "All", "tgrp.fun", "screen.glmnet10"),
c("SL.rpart", "All", "tgrp.fun", "screen.glmnet10"))

set.seed(99)
start <- proc.time()
fitSL.data.CV <- CV.SuperLearner(Y=newdat[,1],X=newdat[,-1], V=10, SL.library=alg_screen, verbose = TRUE, method = "method.NNLS", family = gaussian(), parallel="multicore")
total <- proc.time()-start
print(total)

cvrisk.SL <- mean((newdat[,1]-fitSL.data.CV$SL.predict)^2) #CV risk for super learner
print(cvrisk.SL) 

save(fitSL.data.CV, file="/data/markscan_authorized_users/bergquist/EI/fullCVfit_sim.Rdata")
save(cvrisk.SL, file="/data/markscan_authorized_users/bergquist/EI/fullCVriskSL_sim.Rdata")