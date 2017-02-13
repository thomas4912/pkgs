abremLoglike<-function(x, par, dist="weibull" )  {				
## check basic format of x				
				
	if(class(x)!="data.frame") {stop("abremLoglike takes a structured dataframe input, use mleframe")}			
	if(ncol(x)!=3)  {stop("abremLoglike takes a structured dataframe input, use mleframe")}			
	xnames<-names(x)			
	if(xnames[1]!="left" || xnames[2]!="right"||xnames[3]!="qty")  {			
		 stop("abremLoglike takes a structured dataframe input, use mleframe")  }		
## test for any na's and stop, else testint below will be wrong				
				
				
## need this length information regardless of input object formation				
	testint<-x$right-x$left			
	failNDX<-which(testint==0)			
	suspNDX<-which(testint<0)			
	Nf<-length(failNDX)			
	Ns<-length(suspNDX)			
	discoveryNDX<-which(x$left==0)			
	Nd<-length(discoveryNDX)			
	intervalNDX<-which(testint>0)			
	interval<-x[intervalNDX,]			
	intervalsNDX<-which(interval$left>0)			
	Ni<-length(intervalsNDX)			
				
				
## need to stop if Nf<1?				
## or Nf+Ndi-Nd <3?				
				
## further validate the input arguments for non-frame.fsiq object				
	if(length(attributes(x)$fsiq)!=1)  {			
				
				
				
## stop if Nf+Ns+Ndi != nrow(x)				
	if( (Nf+Ns+Nd+Ni) != nrow(x))  {			
		stop("invalid input dataframe")		
	}			
				
## rebuild input vector from components, just to be sure				
	fsiq<-rbind(x[failNDX,], x[suspNDX,], x[discoveryNDX,], interval[intervalsNDX,])			
## end input validation code				
	}else{			
		fsiq<-x		
	}			
				
## now form the arguments for C++ call				
## no data limitation applies to getting a Loglikelihood value	
##	if((Nf+Ni)<3)  {stop("insufficient failure data")}	

	fsdi<-NULL
	if( (Nf+Ns)>0 )  {
		fsdi<-fsiq$left[1:(Nf + Ns)]
	}	
	if(Nd>0) {		
		fsdi<-c(fsdi,fsiq$right[(Nf + Ns + 1):(Nf +  Ns + Nd)])	
	}		
	if(Ni>0)  {		
		fsdi<-c(fsdi, fsiq$left[(Nf + Ns + Nd + 1):nrow(fsiq)], 	
			  fsiq$right[(Nf + Ns + Nd + 1):nrow(fsiq)])	
	}
	
	q<-fsiq$qty			
## third argument will be c(Nf,Ns,Nd,Ni)				
	N<-c(Nf,Ns,Nd,Ni)	

## establish distribution number
	if(tolower(dist)=="weibull"	)  {
	dist_num=1
	}else{
		if(tolower(dist)=="lognormal")  {
			dist_num=2
		}else{
			stop("distribution not resolved")
		}
	}

	MLEclassList<-list(fsdi=fsdi,q=q,N=N)
								
	outval<-.Call("MLEloglike",MLEclassList,par,dist_num, package="abremDebias")
				
				
				
			
outval			
}