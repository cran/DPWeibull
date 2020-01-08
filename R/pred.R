predict.ddp<-function(object,xpred,alpha=0.05,tpred=NULL,...){
if(is.null(tpred)){
tpred<-object$predtime
}
tpred<-tpred/object$high.pct*10

xpredtemp<-rbind(xpred,object$x)
xpredframe<-data.frame(xpredtemp)
xpredtemp<-model.matrix(~xpredtemp,xpredframe,contrasts)

xpred<-matrix(xpredtemp[1:nrow(xpred),],nrow=nrow(xpred))
xpred<-matrix(xpred[,-1],nrow=nrow(xpred))
xpred<-(xpred-matrix(rep(object$xmean,times=nrow(xpred)),nrow=nrow(xpred), byrow=TRUE))/matrix(rep(2*object$xsd,times=nrow(xpred)),nrow=nrow(xpred), byrow=TRUE)	
predobject<-.Call('DPWeibull_predreg', PACKAGE = 'DPWeibull',
object$alpharec,object$lambdascaled,object$betarec,xpred,tpred,alpha)
class(predobject)<-"predddp"
predobject$tpred<-tpred
predobject$alpha<-alpha
predobject
}

predict.ddpcomp<-function(object,xpred,alpha=0.05,tpred=NULL,...){
if(is.null(tpred)){
tpred<-object$predtime
}
tpred<-tpred/object$high.pct*10

xpredtemp<-rbind(xpred,object$x)
xpredframe<-data.frame(xpredtemp)
xpredtemp<-model.matrix(~xpredtemp,xpredframe,contrasts)
xpredtemp<-xpredtemp[,-1]

xpred<-matrix(xpredtemp[1:nrow(xpred),],nrow=nrow(xpred))

xpred<-(xpred-matrix(rep(object$xmean,times=nrow(xpred)),nrow=nrow(xpred), byrow=TRUE))/matrix(rep(2*object$xsd,times=nrow(xpred)),nrow=nrow(xpred), byrow=TRUE)	
predobject<-.Call('DPWeibull_predcompreg', PACKAGE = 'DPWeibull',
object$alpharec1,object$lambdascaled1,object$betarec1,
object$alpharec2,object$lambdascaled2,object$betarec2,object$prec,
xpred,tpred,alpha)
class(predobject)<-"predddpcomp"
predobject$tpred<-tpred
predobject$alpha<-alpha
predobject
}
