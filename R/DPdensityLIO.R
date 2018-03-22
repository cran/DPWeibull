normwrapper <- function(y,y50=NULL,y95=NULL,mcmc=list(nburn = 1000, nsave = 1000, nskip = 10, ndisplay = 100),ngrid=1000,grid=NULL) {
    # Standardize data using median as center, scale derived from 95th percentile
    Y <- as.matrix(y)
    n <- nrow(Y)
    p <- ncol(Y)
  
    if(is.null(y50)) {m <- colMedians(Y)}
    else             {m <- y50}
    if(is.null(y95)) {c <- colQuantiles(Y,probs=0.95)}
    else             {c <- y95}
    
    s <- (c-m)/2
    M <- kronecker(rep(1,n),t(m))
    if(p==1)         {R <- s}
    else             {R <- diag(s)}
    Z <- (Y - M) %*% solve(R)
    
    if(!is.null(grid)) {
        mg <- kronecker(rep(1,nrow(grid)),t(m))
        if(p==1) {grid <- (grid - mg)%*%diag(s^(-1))}
        else {grid <- (grid - mg)%*%diag(s^(-1))}
    }
    
    # LIO prior specification:
    # Proportion of data coverage used for Chebyshev bound: P(|z|<=k) >= pi
    pi <- 0.99
    
    # Var(mu_i) = v^2*I
    v <- sqrt(p*(n-1)/(n*qchisq(pi,p)*(1-pi)))
    
    # DPdensity hyperparameters
    m1 <- 0
    tau1 <- 3
    tau2 <- v^2
    nu1 <- p+2
    nu2 <- p
    psiinv2 <- diag(rep(p,p))
    a0 <- 1
    b0 <- 1
    
    
    # Run DPdensity on standardized data
    if(p==1)  {prior <- list(a0=a0,b0=b0,m1=m1,nu1=nu1,nu2=nu2,psiinv2=psiinv2,tau1=tau1,tau2=tau2)}
    else {prior <- list(a0=a0,b0=b0,m2=m1,s2=10^(-9),nu1=nu1,nu2=nu2,psiinv2=psiinv2,tau1=tau1,tau2=tau2)}
    fit <- DPdensity(y=Z,prior=prior,mcmc=mcmc,status=TRUE,ngrid=ngrid,grid=grid)

    fit$y <- Y
    
    
    # Convert MCMC samples back
    if(p==1) {
      mu <- fit$save.state$randsave[,(2*1:n-1)]*s+m
      sigma <- sqrt(fit$save.state$randsave[,(2*1:n)])*s
      fit$save.state$randsave[,(2*1:n-1)] <- mu
      fit$save.state$randsave[,(2*1:n)] <- sigma^2
    }
    else {
      k <- p*(p+3)/2
      
      idx <- sort(kronecker(rep(1,p),seq(0,(n-1)*k,k)))
      muidx <- idx + kronecker(rep(1,n),1:p)
      sigmaidx <- setdiff(1:(k*n),muidx)
      r <- (s%*%t(s))[lower.tri(s%*%t(s),diag=TRUE)]
      
      mu <- fit$save.state$randsave[,muidx]%*%kronecker(diag(rep(1,n)),R) + kronecker(matrix(1,nrow=mcmc$nsave,ncol=n),t(m))
      fit$save.state$randsave[,muidx] <- mu
      sigma <- fit$save.state$randsave[,sigmaidx]%*%kronecker(diag(rep(1,n)),diag(r))
      fit$save.state$randsave[,sigmaidx] <- sigma
      
    }
    
    # Extract DPdensity parameter estimates for standardized data
    if(p<=2 & length(fit$grid1)>0) {
        if(p==1) {
          zgrid <- fit$x1
	  ygrid <- zgrid*s+m
	  fit$grid1 <- ygrid 
          fit$x1 <- fit$grid1 
	  densgrid <- fit$dens/s
	  fit$fun1 <- densgrid
          fit$dens <- fit$fun1
        }
      else {
          z1grid <- fit$x1
          z2grid <- fit$x2
	  ygrid1 <- z1grid*s[1]+m[1]
	  fit$grid1 <- ygrid1
          fit$x1 <- fit$grid1
	  ygrid2 <- z2grid*s[2]+m[2]
	  fit$grid2 <- ygrid2
          fit$x2 <- fit$grid2 
	  densgrid <- fit$dens/prod(s)
          fit$dens <-densgrid 
	  fit$fun1 <- fit$fun1/s[1]
          fit$fun2 <- fit$fun2/s[2]
      }
    }
    
    return(fit)
}



