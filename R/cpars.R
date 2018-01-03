
#' Forces correlation among operating model parameters for M, K, Linf and L50
#'
#' @description Uses typical correlations among estimated parameters to generate realistic samples for natural mortality rate (M), growth rate (K), maximum length (Linf) and length at 50% maturity (L50), these are placed in the cpars slot
#' @param OM An operating model object with M, growth, stock-recruitment and maturity parameters specified.
#' @param nsim The number of simulated values to create (note that OM@nsim will be used preferentially). 
#' @param plot Should the sampled parameters and distributions be plotted?
#' @return An object of class OM with a populated (or appended) cpars slot
#' @author T. Carruthers (Canadian DFO grant)
#' @importFrom mvtnorm rmvnorm
#' @export 
#' @examples
#' testOM<-ForceCor(testOM)
ForceCor<-function(OM,nsim=48,plot=T){
  
  if("nsim"%in%slotNames(OM))nsim<-OM@nsim
  if("seed"%in%slotNames(OM))set.seed(OM@seed)
  
  OM@nsim<-nsim
  
  colline=makeTransparent('blue',60)
  lwdline=4
  histcol='black'
  # Estimation cross correlation is typically -0.9 for Linf and K
  sigma<-matrix(c(1,-0.9,-0.9,1),nrow=2)
  
  # Other Parameters correlation  from meta analysis (technically incorrect in the sapce of other v
  #                  M     K    F_linf    Linf
  sigma=matrix(c(1.000,  0.721, -0.568, -0.721,
                 0.721,  1.000, -0.107, -0.910,
                -0.568, -0.107, 1.000,  0.407,
                -0.721, -0.910, 0.407,  1.000),nrow=4)
  
  
  means<-c(mean(OM@M),mean(OM@K),mean(OM@L50),mean(OM@Linf))
  
  sim<-as.data.frame(mvtnorm::rmvnorm(nsim,rep(0,4),sigma))
  sim<-sim-rep(apply(sim,2,mean),each=nsim)#mean 0
  sim<-sim/rep(apply(sim,2,sd),each=nsim)# cv=1
  cvs<-c(OM@M[2]-OM@M[1],OM@K[2]-OM@K[1],OM@L50[2]-OM@L50[1],OM@Linf[2]-OM@Linf[1])/(1.96*2)/means
  sim<-exp(sim*rep(cvs,each=nsim))
  sim<-sim/rep(apply(sim,2,mean),each=nsim)
  sim<-sim*rep(means,each=nsim)
  
  if(plot){ 
    par(mfrow=c(4,4),mai=c(0.3,0.3,0.4,0.05),omi=c(0.02,0.02,0.3,0.02))
    
    labs<-c("M","K","L50","Linf")
    bounds<-matrix(c(OM@M,OM@K,OM@L50,OM@Linf),nrow=2)
    
    for(i in 1:4){
      
      for(j in 1:4){
        
        if(i == j){
          
          if(i==1){
           
            hist(sim[,1],main="Natural mortality rate (M)",col=histcol,border='white',xlab="",axes=F)
            axis(1)
            abline(v=OM@M,col=colline,lwd=lwdline)
            
          }else if(i==2){
            
            hist(sim[,2],main="Growth rate (K)",col=histcol,border='white',xlab="",axes=F)
            axis(1)
            abline(v=OM@K,col=colline,lwd=lwdline)
            
            
          }else if(i==3){
            
            hist(sim[,3],main="Length at 50% maturity (L50)",col=histcol,border='white',xlab="",axes=F)
            axis(1)
            abline(v=OM@L50,col=colline,lwd=lwdline)
            
            
          }else{
            
            hist(sim[,4],main="Maximum length (Linf)",col=histcol,border='white',xlab="",axes=F)
            axis(1)
            abline(v=OM@Linf,col=colline,lwd=lwdline)
            
          }
          
        }else{ # not positive diagonal
          
          plot(sim[,j],sim[,i],axes=F,col="white")
          polygon(bounds[c(1,1,2,2),j],bounds[c(1,2,2,1),i],col=colline,border="white")
          points(sim[,j],sim[,i],pch=19)
          
          axis(1)
          axis(2)
          
        }
        
      }
    }
    
    mtext("Sampled parameters and cross-correlations",3,font=2,outer=T)
    
  } # if plot
  
  names(sim)<-c("M","K","L50","Linf")
  
  OM@cpars$M<-sim$M
  OM@cpars$K<-sim$K
  OM@cpars$L50<-sim$L50
  OM@cpars$Linf<-sim$Linf
  
  OM
  
}


