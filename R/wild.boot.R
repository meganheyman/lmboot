wild.boot <-
function(formula, B=1000, data=NULL, seed=NULL, bootDistn="normal"){

  ###################################################
  ## Checks for function inputs                    ##
  ###################################################

  if(inherits(formula, "formula")==FALSE){
    stop("The input model must be a formula. \n")
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


  

  if(!( bootDistn %in% c("normal","uniform","exponential","laplace","lognormal",
                         "gumbel","t5","t8","t14") )){
    stop("Invalid value for bootDistn.")
  }
  

  set.seed(seed)








  #######################################################
  ## Least Squares Fit                                 ##
  #######################################################
  obsDataregFit <- lm(formula, data=data)                   #fit the linear model specified in formula input
  estParam <- matrix(obsDataregFit$coef, ncol=1)            #keep the param. estimates in a vector
  obsDataResid <- as.vector(residuals(obsDataregFit))       #keep the original residuals
  ParamNames <- names(obsDataregFit$coefficients)           #keep the coefficient name/association
  rownames(estParam) <- ParamNames                          #name the rows for the parameters so we know what they are
  modelMat <- model.matrix(obsDataregFit)                   #model matrix (X)
  hatMat <- solve(t(modelMat) %*% modelMat) %*% t(modelMat) #projection matrix (X^TX)^-1 X^T



  ######################################################
  ## Bootstrap                                        ##
  ######################################################
  ##Objects to keep Bootstrap Observations
  bootEstParam <- matrix(NA, nrow=B, ncol=dim(estParam)[1])  #bootstrap param. estimates
  colnames(bootEstParam) <- ParamNames


  if(bootDistn=="normal"){
    bootWtMatrix <- matrix(rnorm(B*n, mean=0, sd=1), nrow=B, ncol=n)
  }else if(bootDistn=="uniform"){
    bootWtMatrix <- matrix(runif(B*n, min=-sqrt(3), max=sqrt(3)), nrow=B, ncol=n)
  }else if(bootDistn=="exponential"){
    bootWtMatrix <- matrix(rexp(B*n)-1, nrow=B, ncol=n)
  }else if(bootDistn=="laplace"){
    uniforms <- runif(B*n, min=-1/2, max=1/2)
    bootWtMatrix <- matrix(-1/sqrt(2)*sign(uniforms)*log(1-2*abs(uniforms)), nrow=B, ncol=n)
  }else if(bootDistn=="lognormal"){
    sig <- 1
    mu  <- (log(1/(exp(sig)-1))-sig)/2
    bootWtMatrix <- matrix(rlnorm(B*n, meanlog=mu, sdlog=sig)-exp(mu+sig/2), nrow=B, ncol=n)
  }else if(bootDistn=="gumbel"){
    bootWtMatrix <- matrix(evd::rgumbel(B*n, scale=sqrt(6/pi^2), location=sqrt(6/pi^2)*digamma(1)), nrow=B, ncol=n)
  }else if(bootDistn=="t5"){
    bootWtMatrix <- matrix(rt(B*n, df=5)*sqrt(3/5), nrow=B, ncol=n)
  }else if(bootDistn=="t8"){
    bootWtMatrix <- matrix(rt(B*n, df=8)*sqrt(6/8), nrow=B, ncol=n)
  }else if(bootDistn=="t14"){
    bootWtMatrix <- matrix(rt(B*n, df=14)*sqrt(12/14), nrow=B, ncol=n)
  }

  for(i in 1:B){
    bootResid <- matrix(obsDataResid*bootWtMatrix[i, ], ncol=1)      #bootstrap residuals
    bootEstParam[i,]<- as.vector( estParam + hatMat %*% bootResid )  #bootstrap parameter estimates
  }


  #####################################################
  ## Returns
  #####################################################
  structure(invisible(list(bootEstParam=bootEstParam, 
                           origEstParam=estParam, seed=seed, bootDistn=bootDistn)))

}
