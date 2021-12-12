source("R/basic.R")
source("R/dpm.R")
source("R/reg.R")
source("R/compdpm.R")
source("R/compreg.R")


dpweib<-function (formula,data, high.pct=NULL,predtime=NULL,comp=FALSE,alpha= 0.05, simultaneous=FALSE,burnin=8000, iteration=2000,
alpha00=1.354028, alpha0=0.03501257, lambda00=7.181247, alphaalpha=0.2, alphalambda=0.1, a=1, b=1, gamma0=1, gamma1=1,
thin=10, betasl=2.5, addgroup=2){   
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula","data"),names(mf),0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")

	if(identical(attr(y,"type"),"interval")){
		tl<-y[,1]
		tr<-y[,2]
		event<-y[,3]
		if(any(y[,3]==2)){
	  	stop("Sorry! We have not yet added the module for left-truncated data.")
		}
	}else{
    if(identical(attr(y,"model"),"competing.risks")){
      comp=TRUE
      if(attr(y,"cens.type") != "rightCensored"){
        stop("Sorry! We have not yet added the interval censoring module for competing risks data.")	
      }else{
      tevent<-y[,1]
      event<-ifelse(y[,2]==1,y[,3],y[,2])
      }
    }else{
	  tl<-y[,1]
	  tr<-y[,1]
	  event<-y[,2]
    }
	}
  if(is.null(high.pct)){
    censorratio<-sum(event==0)/length(event)
    if(censorratio>0.153){
      high.pct<-findhighpct(censorratio)*ifelse(identical(attr(y,"type"),"interval"),
                                                max(as.vector(y[,1:2]),na.rm=TRUE),
                                                max(y[,1],na.rm=TRUE))
    }else{
      high.pct<-ifelse(identical(attr(y,"type"),"interval"),
                       quantile(as.vector(y[,1:2]),0.95,na.rm=TRUE),
                       quantile(y[,1],0.95,na.rm=TRUE))
    }
  }
  usertime<-TRUE
  if(is.null(predtime)){
   usertime<-FALSE
	 predtime<-(1:41)/41*high.pct
  }
  x <- model.matrix(mt, mf)
  xlevels <- .getXlevels(mt, mf)
  Terms <- delete.response(mt)
  indicator<-apply(x,2,is.binary)

  if(comp==TRUE){
    if (ncol(x)==1) {
      z <- compdpm(tevent,event,high.pct,predtime, burnin,iteration,
            alpha00,alpha0,lambda00,
            alphaalpha,alphalambda,a,b,
            gamma0,gamma1,
            addgroup,thin)
        if(alpha!=0.05){
        	z<-dpmcompdiffalpha(alpha,z)
	      }
        if(simultaneous==TRUE){
        z$CIF1bandl<-confband(alpha,z$CIF1)[1,]
        z$CIF1bandu<-confband(alpha,z$CIF1)[2,]
        z$d1bandl<-confband(alpha,z$d1)[1,]
        z$d1bandu<-confband(alpha,z$d1)[2,]
        z$h1bandl<-confband(alpha,z$h1)[1,]
        z$h1bandu<-confband(alpha,z$h1)[2,]
        z$CIF2bandl<-confband(alpha,z$CIF2)[1,]
        z$CIF2bandu<-confband(alpha,z$CIF2)[2,]
        z$d2bandl<-confband(alpha,z$d2)[1,]
        z$d2bandu<-confband(alpha,z$d2)[2,]
        z$h2bandl<-confband(alpha,z$h2)[1,]
        z$h2bandu<-confband(alpha,z$h2)[2,]
	      }  
	      class(z)<-"dpmcomp"
    }else{
      x<-matrix(x[,-1],nrow=nrow(x))
	    z <- compreg(tevent,event,x,high.pct,predtime,indicator[-1],burnin,iteration,
                  alpha00, alpha0,lambda00,
                  alphaalpha,alphalambda,a,b,
                  gamma0,gamma1, addgroup,
                  thin, betasl)
        z$x<-x
        if(alpha!=0.05){
	      z<-ddpdiffalpha(alpha,z)
	      }
        if(simultaneous==TRUE){
        z$loghrbandl<-matrix(confband(alpha,z$loghr)[1,],byrow=TRUE,nrow=ncol(z$x))/z$xscale
      	z$loghrbandu<-matrix(confband(alpha,z$loghr)[2,],byrow=TRUE,nrow=ncol(z$x))/z$xscale
	  }
	class(z)<-"ddpcomp"
    }
}else{
    if (ncol(x)==1) {
        z <- dpm(tl,tr,event,high.pct,predtime,burnin,iteration,
                  alpha00, alpha0, lambda00,
                  alphaalpha,alphalambda,
                  a,b,addgroup,thin)
	class(z)<-"dpm"
        if(alpha!=0.05){
	z<-dpmdiffalpha(alpha,z)
	}
        if(simultaneous==TRUE){
        z$Sbandl<-confband(alpha,z$S)[1,]
        z$Sbandu<-confband(alpha,z$S)[2,]
        z$dbandl<-confband(alpha,z$d)[1,]
        z$dbandu<-confband(alpha,z$d)[2,]
        z$hbandl<-confband(alpha,z$h)[1,]
        z$hbandu<-confband(alpha,z$h)[2,]
	}    
	}
    else {
      x<-matrix(x[,-1],nrow=nrow(x))
        z <- reg(tl,tr,event,x,high.pct,predtime, indicator[-1],burnin,iteration,
                  alpha00,alpha0,lambda00,
                  alphaalpha,alphalambda,a,b,
                  addgroup,betasl,thin)
        z$x<-x
	class(z)<-"ddp"
        if(alpha!=0.05){
	z<-ddpdiffalpha(alpha,z)
	}
        if(simultaneous==TRUE){
        z$loghrbandl<-matrix(confband(alpha,z$loghr)[1,],byrow=TRUE,nrow=ncol(z$x))/z$xscale
      	z$loghrbandu<-matrix(confband(alpha,z$loghr)[2,],byrow=TRUE,nrow=ncol(z$x))/z$xscale
    	}
    }
   }
   z$usertime<-usertime
   z$alpha<-alpha
   z$simultaneous<-simultaneous
   z$high.pct<-high.pct
   z$predtime<-predtime
   z$Terms<-Terms
   z$xlevels<-xlevels
   z
}
source("R/continue.R")
source("R/pred.R")
source("R/summary.R")
source("R/plot.R")
