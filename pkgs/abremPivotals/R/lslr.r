lslr<-function(x, dist="weibull", npar=2, reg_method="XonY")  {				
	## a convergence limit is fixed here for 3rd parameter  convergence			
	## no longer an argument for the R function, but still an argument to C++ functions			
	limit<-1e-5			
				
	if(is.vector(x))  {			
		stop("use MRR functions for casual fitting, or pre-process with getPPP")		
	}else{			
		if(names(x)[1]=="time"&&names(x)[2]=="ppp")  {		
		## will handle the output from getPPP				
		}else{	
			if(length(x$ppp)<3)  {
				stop("insufficient failure points")
			}else{
				stop("input format not recognized")	
			}
		}		
	}			
				
	casenum<-0			
	if(reg_method=="YonX") casenum=casenum+1			
	if(npar==3) casenum=casenum+2			
	if(dist=="lnorm")casenum=casenum+4			
	if(dist=="gumbel") casenum=casenum+8			
				
				
	resultVec<-.Call("LSLR", x$time, x$ppp, limit, casenum , package="abremPivotals")			
				
	if(casenum < 4) {			
		if(length(resultVec)==3)  {	
			prr<-AbPval(dim(x)[1], resultVec[3])
			outVec<-c(Eta=resultVec[1],Beta=resultVec[2],Rsqr=resultVec[3], AbPval=prr[1])	
		}else{		
			outVec<-c(Eta=resultVec[1],Beta=resultVec[2], t0=resultVec[3],Rsqr=resultVec[4])
			if(resultVec[5]==1)  {			
				warn="3p optimization did not converge"
				attr(outvec,"warning")<-warn
			}			
		}		
	}else{			
		if(casenum < 8) {		
			if(length(resultVec)==3)  {	
				prr<-AbPval(length(x[,1]), resultVec[3],"lnorm")
				outVec<-c(Mulog=resultVec[1],Sigmalog=resultVec[2],Rsqr=resultVec[3], AbPval=prr[1])	
			}else{	
				outVec<-c(Mulog=resultVec[1],Sigmalog=resultVec[2], t0=resultVec[3],Rsqr=resultVec[4])
				if(resultVec[5]==1)  {
					warn="3p optimization did not converge"
					attr(outVec,"warning")<-warn
				}
			}	
		}else{		
			if(length(resultVec)==3)  {	
				outVec<-c(Etalog=resultVec[1],Betalog=resultVec[2],Rsqr=resultVec[3])
			}else{	
				outVec<-c(Etalog=resultVec[1],Betalog=resultVec[2], t0=resultVec[3],Rsqr=resultVec[4])
				if(resultVec[5]==1)  {
					warn="3p optimization did not converge"
					attr(outVec,"warning")<-warn
				}				
			}	
		}		
	}			
				
return(outVec)				
}				
