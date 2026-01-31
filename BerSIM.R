Ber <-
function(x,y,tau=0.5, runs=2000, burn=1000, thin=1) {
    #x:    matrix of predictors.
    #y:    vector of dependent variable. 
    #tau:  expectile level.
    #runs: the length of the Markov chain.
    #burn: the length of burn-in.
    #thin: thinning parameter of MCMC draws
    x <- as.matrix(x)  
    if(ncol(x)==1) {x=x} else {
    x=x
    if (all(x[,2]==1)) x=x[,-2] }

      # Calculate some useful quantities
        n <- nrow(x)
        p <- ncol(x)
   
      # check input
        if (tau<=0 || tau>=1) stop ("invalid tau:  tau should be >= 0 and <= 1. 
               \nPlease respecify tau and call again.\n")
        if(n != length(y)) stop("length(y) not equal to nrow(x)")
        if(n == 0) return(list(coefficients=numeric(0),fitted.values=numeric(0),
              deviance=numeric(0)))
        if(!(all(is.finite(y)) || all(is.finite(x)))) stop(" All values must  be 
              finite and non-missing") 

      # Saving output matrices 
        betadraw  = matrix(nrow=runs, ncol=p)
        MuY  = matrix(nrow=runs, ncol=n)
        VarY  = matrix(nrow=runs, ncol=n)
        sigma2draw = matrix(nrow=runs, ncol=1)

      # Initial valus
        beta   = rep(0.99, p)
        v      = rep(1, n)
        sigma2  = 1

      # Parameters for Inverse-Gamma (on sigma^2)
        a0 <- 0.01
        b0 <- 0.01
        c0 <- 0.1
        d0 <- 0.1
        c1 <- 0.1
        d1 <- 0.1
        sigmaBeta=0.002
	  sigmaGamma=0.5

gamma= 0.5
lambda= 0.5
	acceptrate<-0

Cn<-matrix(gamma*exp(-(outer(x%*%beta,x%*%beta,FUN="-"))^(2)),ncol=n)
ETAn= matrix(mvrnorm(n=1,rep(0,n),Cn),ncol=1,nrow=n) 

      # Start the algorithm
        for (iter in 1: runs) {

     # Update latent weights (Data Augmentation)
       residuals <- y - ETAn

     # Asymmetric normal representation weights
       v <- numeric(n)
       for (i in 1:n) {
       if (residuals[i] >= 0) {
        v[i] <- tau
        } else {
        v[i] <- 1 - tau
        }
       }

 		E= sigma2*diag(c(v))
          
      # Draw sigma
        shape =   n/2 + p + a0
        rate  = sum((y - ETAn)^2/(2*v)) + lambda*sum(abs(beta)) + b0
        sigma2 = 1/rgamma(1, shape= shape, rate= rate)


      # Draw beta
       
        betatemp<-matrix(beta+sqrt(sigmaBeta)*matrix(rnorm(p),ncol=1),ncol=1)      #new beta#
        betatemp=c(betatemp)
        Cntemp<-matrix(gamma*exp(-(outer(x%*%betatemp,x%*%betatemp,FUN="-"))^(2)),ncol=n)
        alpha1<-exp(0.5*determinant(Cn+0.00001*diag(n),logarithm=TRUE)$modulus-0.5*determinant(Cntemp+0.00001*diag(n),
        logarithm=TRUE)$modulus-lambda*(sum(abs(betatemp))-sum(abs(beta)))/sigma2-0.5*
        t(ETAn)%*%(solve(Cntemp+0.001*diag(n))-solve(Cn+0.001*diag(n)))%*%(ETAn))

cat(alpha1," ")                        
		if (runif(1)<min(1,alpha1))
		{
			beta<-betatemp
			beta<-sign(beta[2])*beta/sqrt(sum(beta^2)) 
			Cn<-Cntemp
			acceptrate<-acceptrate+1
}else{
	      	beta <-sign(beta[2])*beta/sqrt(sum(beta^2))}

                  lambda<-rgamma(1,c1+p,rate=sum(abs(beta))/sigma2+d1)

            gammatemp<-exp(rnorm(1,mean=log(gamma),sd=sigmaGamma))       
	Cntemp<-(gammatemp/gamma)*Cn
		alpha2<-exp(0.5*determinant(Cn,logarithm=TRUE)$modulus-0.5*determinant(Cntemp,logarithm=TRUE)$modulus+
                    (a0+1+1)*(log(gamma)-log(gammatemp))-d0*(1/gammatemp-1/gamma)-
					+ 0.5*t(y)%*%(solve(E+Cntemp)-solve(E+Cn))%*%(y))
		
		
		if (runif(1)<min(1,alpha2))
		{ 
			gamma<-gammatemp
			Cn<-Cntemp
		}
		else {gamma<-gamma
                   Cn= Cn}





           
#print(betatemp)
#print(Cntemp)
#print(alpha1)

      # Sort beta and sigma
        betadraw[iter,]  = beta
        MuY[iter,  ]     = ETAn
        sigma2draw[iter,] = sigma2
}
        coefficients =apply(as.matrix(betadraw[-(1:burn), ]),2,mean)
        names(coefficients)=colnames(x)
        if (all(x[,1]==1))  names(coefficients)[1]= "Intercept"  

        result <- list(beta = betadraw[seq(burn, runs, thin),],
                       MuY = MuY[seq(burn, runs, thin),],
        sigma2  = sigma2draw[seq(burn, runs, thin),],
        y=y,    
        coefficients=coefficients)
    
      return(result)
      result
}

