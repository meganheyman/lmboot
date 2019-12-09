paired.boot <-
function(formula, B=1000, seed=NULL, data=NULL){

  ###################################################
  ## Checks for function inputs                    ##
  ###################################################

  if(inherits(formula, "formula")==FALSE){
    stop("The input model must be a formula.\n")
  }

  full.model.frame <- model.frame(formula, data=data, na.action = na.pass) #get model variables in data.frame
  resp <- model.response(full.model.frame)                                 #get the response variable
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




  modelMat <- model.matrix(formula, data=data)                #get the model matrix

  if(dim(modelMat)[2] <= 0){
    stop("The model has no predictors or intercept.\n")
  }
  modelqr <- qr(modelMat)                                     #perform QR decomposition on model matrix for checks
  model.pivot <- modelqr$pivot[1:modelqr$rank]
  if (ncol(modelMat) > modelqr$rank) {
    warning("The design matrix isn't full column rank.\n")
  }

  if(dim(modelMat)[1]!=length(resp)){
    stop("Predictors must not have any missing values.\n")
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





  #######################################################
  ## Least Squares Fit                                 ##
  #######################################################
  obsDataregFit <- lm(formula, data=data)             #fit the linear model specified in formula input
  estParam <- matrix(obsDataregFit$coef, ncol=1)      #keep the param. estimates in a vector
  obsDataResid <- as.vector(residuals(obsDataregFit)) #keep the original residuals
  ParamNames <- names(obsDataregFit$coefficients)     #keep the coefficient name/association
  rownames(estParam) <- ParamNames                    #name the rows for the parameters so we know what they are
  modelMat <- model.matrix(obsDataregFit)                   #model matrix (X)
  hatMat <- solve(t(modelMat) %*% modelMat) %*% t(modelMat) #projection matrix (X^TX)^-1 X^T



  ######################################################
  ## Bootstrap                                        ##
  ######################################################
  ##Objects to keep Bootstrap Observations
  bootEstParam <- matrix(NA, nrow=B, ncol=dim(estParam)[1])  #bootstrap param. estimates
  colnames(bootEstParam) <- ParamNames

  for(i in 1:B){
    index.sel <- sample(1:n, size=n, replace=TRUE)                  #sample indices for observations in bootstrap
    bootResp <- matrix(resp[index.sel], ncol=1)                     #bootstrap response
    bootPred <- matrix(modelMat[index.sel, ], ncol=dim(modelMat)[2])#bootstrap predictors

    while(ncol(bootPred) > qr(bootPred)$rank){
      index.sel <- sample(1:n, size=n, replace=TRUE)         #attempt to ensure that X^TX is non-singular
      bootResp <- matrix(resp[index.sel], ncol=1)            #re-assign bootstrap response
      bootPred <- modelMat[index.sel, ]                      #re-assign bootstrap predictors
    }

    bootEstParam[i,]<- as.vector(solve(t(bootPred) %*% bootPred) %*% t(bootPred) %*% bootResp) #boot param est
  }


  #####################################################
  ## Returns
  #####################################################
  structure(invisible(list(bootEstParam=bootEstParam, 
                           origEstParam=estParam, seed=seed)))

}
