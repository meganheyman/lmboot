ANOVA.boot <-
function(formula, B=1000, type="residual", wild.dist="normal",  
                       seed=NULL, data = NULL, keep.boot.resp=FALSE){

  ############################################
  ##Check for function inputs               ##
  ############################################

  if(inherits(formula, "formula")==FALSE){
    stop("The input model must be a formula.\n")
  }
  
  full.model.frame <- model.frame(formula, data=data, na.action = na.pass) #get model variables in data.frame
  resp <- model.response(full.model.frame)                                 #get the response variable
  p <- dim(full.model.frame)[2]-1                                          #get the number of covariates
  n <- length(resp)                                                        #get the number of observations

  if(is.matrix(resp)!=TRUE && is.vector(resp)!=TRUE){
    stop("Response must be a vector or matrix.\n")
  }
  else if((dim(resp)[1]==0 || dim(resp)[2]==0) && length(resp)==0){
    stop("Response must have entries.\n")
  }
  else if(mode(resp)!="numeric"){
    stop("Response must be of type numeric.\n")
  }
  else if(anyNA(resp)==TRUE){
    stop("Response must not have any missing values.\n")
  }


  if (p > 2) {
    warning("This function has only been fully tested for one-way and two-way ANOVA.\n")
  }else if(p==0){
    stop("The model must have predictor variables. \n")
  }
  
  modelMat <- model.matrix(formula, data=data)[,]             #get the model matrix
  modelqr <- qr(modelMat)                                     #perform QR decomposition on model matrix for checks
  model.pivot <- modelqr$pivot[1:modelqr$rank]
  if (ncol(modelMat) > modelqr$rank) {
    warning("The design matrix isn't full column rank.\n")
  } else if (modelqr$rank == 1) {
    stop("Model is same as response~1. Halting ANOVA test.\n")
  }
  if(anyNA(modelMat)==TRUE){
    stop("Predictors must not have any missing values.\n")
  }
  
  if(type!="residual" & type!="wild"){
    stop("Only residual bootstrap or wild boostrap is allowed for type.\n")
  } 

  if(!( wild.dist %in% c("normal","uniform","exponential","laplace","lognormal",
                         "gumbel","t5","t8","t14") )){
    stop("Invalid value for wild.dist.\n")
  }
  
  if(mode(B)!="numeric"){
    stop("Number of bootstrap samples, B, must be of type numeric.\n")
  }
  else if(is.atomic(B)!=TRUE){
    stop("Number of bootstrap samples, B, must be a constant.\n")
  }
  else if(is.null(B)==TRUE){
    stop("Number of bootstrap samples, B cannot be NULL.\n")
  }
  else if( B < n){
    warning("Number of bootstrap samples is recommended to be more than the number of observations.\n")
  }
  
  if(is.null(seed)==TRUE){
    seed <- sample(seq(1,100000000), size=1)
  }
  else{
    if(mode(seed)!="numeric"){
      stop("The seed must be of type numeric.\n")
    }
    else if(is.atomic(seed)!=TRUE){
      stop("The seed must be a constant.\n")
    }
  }


  set.seed(seed)
  
  ##################################################
  ## ANOVA Model Fit                              ##
  ##################################################
  ### SOME OF THESE VARIABLES MAY NOT BE NEEDED???
  anovaMod <- aov(formula, data = data)               #ANOVA model (original data)
  SS.list <- summary(anovaMod)[[1]][,2]               #ANOVA sum of squares (original data)
  df.list <- summary(anovaMod)[[1]][,1]               #ANOVA degrees of freedom (original data)
  termNames <- row.names(summary(anovaMod)[[1]])      #names of SS terms (original data)
  df.len <- length(df.list)                           #number of SS terms (original data)
  obsFStats <- summary(anovaMod)[[1]][,4]             #get the F test stats from orig fit
  obsFStats <- obsFStats[-length(obsFStats)]          #remove NA next to MSE
  quadMat <- solve(t(modelMat) %*% modelMat) %*% t(modelMat) #projection matrix (X^TX)^-1 X^T (original data)
  y.hat <- modelMat %*% (quadMat %*% resp)                   #predicted values (original data)
  obsDataResid <- resp - y.hat                               #residuals (original data)

  modelMat.list <- list()                             #keep columns of model matrix assoc. with each factor, progressively
  modelMat.list[[1]] <- matrix(modelMat[, 1], ncol=1) #this is the intercept in the model matrix
  count <- 1
  for (i in 1:(df.len-1)) {
    count <- count + df.list[i]
    modelMat.list[[i+1]] <- modelMat[, model.pivot[1:count]]
  }





  ##################################################
  ## Bootstrap here                               ##
  ##################################################
  #keep bootstrap response if requested
  if (keep.boot.resp == TRUE) { bootResponseMatrix.list <- list() }
  
  bootFStats <- matrix(NA, nrow=B, ncol=(df.len-1))   #bootstrap F stats for each factor & potentially interaction
  bootSSEMat <- matrix(NA, nrow=B, ncol=1)            #bootstrap SSE for full model          
  bootSSTrMat <- matrix(NA, nrow=B, ncol=(df.len-1))  #bootstrap SSTr for each factor & potentially interaction
 

 
###CODE IS CHECKED TO HERE -- can we get bootstrap out of a loop?
### I'm not sure if the observed residuals need to be in here or if we need residuals associated with one of the inside models
  for(i in 1:(df.len-1)) {
    nullMat <- modelMat.list[[i]]                                #model matrix (X) for null model
    nullQuadMat <- solve(t(nullMat) %*% nullMat) %*% t(nullMat)  #projection matrix (X^T X)^{-1}X^T for null model
    y.hatnull <- nullMat %*% (nullQuadMat %*% resp)              #predicted values under null model
    modelMatTest <- modelMat.list[[i+1]]                         #model matrix (X) for factor being tested
    quadMatTest <- solve(t(modelMatTest) %*% modelMatTest) %*% t(modelMatTest)  #projection matrix (X^T X)^{-1}X^T for factor being tested
    y.hattest <- modelMatTest %*% (quadMatTest %*% resp)         #predicted values under model for factor being tested
    testResid <- resp - y.hattest                                #residuals associated with model for factor being tested
    df.Tr <- df.list[i]                                          #degrees of freedom for factor being tested
    df.E <- n - (1 + sum(df.list[1:i]))                    #degrees of freedom, error in test
    
    if (keep.boot.resp == TRUE) { bootResponseMatrix.list[[i]] <- matrix(NA, nrow=B, ncol=n) }
    


    for(j in 1:B){
      if (type == "residual") {
        bootResid <- matrix(sample(testResid, replace=TRUE), ncol=1)  #bootstrap residuals
        yb <- matrix(y.hatnull + bootResid, nrow=n, ncol=1)        #bootstrap response under H0
      } else {
        if(wild.dist=="normal"){
          bootResid <- matrix(testResid*rnorm(n, mean=0, sd=1), ncol=1)
        }else if(wild.dist=="uniform"){
          bootResid <- matrix(testResid*runif(n, min=-sqrt(3), max=sqrt(3)), ncol=1)
        }else if(wild.dist=="exponential"){
          bootResid <- matrix(testResid*(rexp(n)-1), ncol=1)
        }else if(wild.dist=="laplace"){
          uniforms <- runif(n, min=-1/2, max=1/2)
          bootResid <- matrix(testResid*(-1/sqrt(2)*sign(uniforms)*log(1-2*abs(uniforms))), ncol=1)
        }else if(wild.dist=="lognormal"){
          sig <- 1
          mu  <- (log(1/(exp(sig)-1))-sig)/2
          bootResid <- matrix(testResid*(rlnorm(n, meanlog=mu, sdlog=sig)-exp(mu+sig/2)), ncol=1)
        }else if(wild.dist=="gumbel"){
          bootResid <- matrix(testResid*(evd::rgumbel(n, scale=sqrt(6/pi^2), location=sqrt(6/pi^2)*digamma(1))), ncol=1)
        }else if(wild.dist=="t5"){
          bootResid <- matrix(testResid*rt(n, df=5)*sqrt(3/5), ncol=1)
        }else if(wild.dist=="t8"){
          bootResid <- matrix(testResid*rt(n, df=8)*sqrt(6/8), ncol=1)
        }else if(wild.dist=="t14"){
          bootResid <- matrix(testResid*rt(n, df=14)*sqrt(12/14), ncol=1)
        }
        yb <- matrix(y.hatnull + bootResid, nrow=n, ncol=1)      #bootstrap response under H0
      }


      if (keep.boot.resp == TRUE) { bootResponseMatrix.list[[i]][j, ] <- yb }

      y.hatnullb <- nullMat %*% (nullQuadMat %*% yb)                   #bootstrap predicted values under H0
      y.hattestb <- modelMatTest %*% (quadMatTest %*% yb)              #bootstrap predicted values for factor eing tested
      SSE <- t(yb - y.hattestb) %*% (yb - y.hattestb)                  #sum of squares, error
      SST <- t(yb - y.hatnullb) %*% (yb - y.hatnullb)                  #sum of squares, total
      SSTr <- SST - SSE                                                #sum of squares, treatment
      MSE <- SSE / df.E                                                #mean squares, error
      MSTr <- SSTr / df.Tr                                             #mean squares, treatment
      Ftest <- MSTr/MSE                                                #test statistic (under H0)
      bootSSTrMat[j, i] <- SSTr
      if (i == (df.len-1)){
        bootSSEMat[j, 1] <- SSE
      }
      bootFStats[j, i]<- Ftest
    }
 

    
  }
  
  pvalues <- rep(NA, df.len-1)
  for (i in 1:(df.len-1)) {
    pvalues[i] <- mean(bootFStats[,i] > obsFStats[i])
  }
  
##Potentially need to re-name or re-format output here
  if (keep.boot.resp == TRUE) {
    structure(invisible(list(terms = termNames,df = df.list, bootFStats=bootFStats, 
                             origSSE=SS.list[length(SS.list)], origSSTr=SS.list[1:(length(SS.list)-1)],
                             bootSSE=bootSSEMat, bootSSTr=bootSSTrMat,
                             origFStats=obsFStats, "p-values" = pvalues,
                             bootResponse=bootResponseMatrix.list)))
  } else {
    structure(invisible(list(terms = termNames, 
                             df= df.list, 
                             origFStats=obsFStats, 
                             origSSE=SS.list[length(SS.list)], 
                             origSSTr=SS.list[1:(length(SS.list)-1)],
                             bootFStats=bootFStats,
                             bootSSE=bootSSEMat, 
                             bootSSTr=bootSSTrMat,
                             "p-values" = pvalues)))
  }
}
