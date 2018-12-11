###################################
###Caio Kenup######################
###Contact: caio.kenup@gmail.com###
###2018-12-11######################
###################################

#Functions to run visual diagnostics of model fit (for GLMs and GLMMs)

diagnostics.glmm<-function(mod,n){
	  
	  ###mod: a glmer model
	  ###n: number of simulations to run
	  
  	ri<-sort(resid(mod)) #residuals from model
    ri0<-ri[ri<0]        #only negative residuals
    ri1<-ri[ri>=0]       #only positive residuals
    
  	ft<-fitted(mod)      #fitted values from model
  	
  	#Create a 'null model' matrix, with probabilities equal to fitted values
    M<-matrix(rbinom(
  	  n=length(ft)*n,size=1,p=ft),
  	  nrow=length(ft), #rows equal to the number of observations
  	  ncol=n)          #columns equal to the number of simulations
    
  	resp.var<-colnames(mod@frame)[1] #get response variable as character
  	
  	#Create new formula (respose variable as a function of 'x', a 'null' variable)
  	new.formula<-as.formula(gsub(x=(mod@call%>%as.character)[2],resp.var,'x'))
  	
  	
  	#for each simulations, run a GLMM
  	X<-lapply(1:ncol(M),function(y){
	   x<-M[,y]
	   return(
	    glmer(formula = new.formula,
	          data = cbind(mod@frame,x=x),
	          family = binomial())
	  )})
	
  	#get residuals from each simulated GLMM
  	R<-sapply(X,function(x){sort(resid(x))}) 
  	
  	#Separate residual matrices into negative and positive residuals
  	R0<-matrix(
  	  apply(
  	    R,MARGIN=2,FUN=function(x){
  	      quantile(x[x<0],ppoints(ri0,a=1))}),ncol=n,nrow=length(ri0))
  	R1<-matrix(
  	  apply(
  	    R,MARGIN=2,FUN=function(x){
  	      quantile(x[x>=0],ppoints(ri1,a=1))}),ncol=n,nrow=length(ri1))
  	  	
  	#Get 50, 2.5 and 97.5 percentiles for negative and positive residuals
  	medianSim0<-apply(R0,1,median)
  	medianSim1<-apply(R1,1,median)
  	lclSim0<-apply(R0,1,quantile,.025)
  	uclSim0<-apply(R0,1,quantile,.975)
  	lclSim1<-apply(R1,1,quantile,.025)
  	uclSim1<-apply(R1,1,quantile,.975)  	
  	
  	#Created data.frame of fitted residuals and simulated residuals
  	df<-data.frame(
  	    Fitted_Residuals=c(ri0,ri1),
  	    Simulated_Residuals=c(medianSim0,medianSim1),
  	    UCLSim=c(uclSim0,uclSim1),LCLSim=c(lclSim0,lclSim1))
  	
  	#Plot this data.frame to assess model fit
  	g<-ggplot(df,aes(x=Simulated_Residuals,y=Fitted_Residuals))+
  	  geom_point()+
  	  geom_line(aes(x=Simulated_Residuals,y=LCLSim))+
  	  geom_line(aes(x=Simulated_Residuals,y=UCLSim))+
  	  geom_abline()
  	return(g)
	}
	
diagnostics.glm<-function(mod,n){
	  ri<-sort(resid(mod))
	  ri0<-ri[ri<0]
	  ri1<-ri[ri>=0]
	  
	  ft<-fitted(mod)
	  
	  M<-matrix(rbinom(
	    n=length(ft)*n,size=1,p=ft),nrow=length(ft),ncol=n)
	  form<-as.character(mod$formula)
	  resp.var<-form[2]
	  
	  new.formula<-paste('x',form[1],form[3:length(form)])%>%
	    as.formula
	  
	  X<-lapply(1:ncol(M),function(y){
	    x<-M[,y]
	    return(
	      glm(formula = new.formula,
	            data = cbind(mod$model,x=x),
	            family = binomial())
	    )})
	  
	  R<-sapply(X,function(x){sort(resid(x))})
	  
	  R0<-matrix(
	    apply(
	      R,MARGIN=2,FUN=function(x){
	        quantile(x[x<0],ppoints(ri0,a=1))}),ncol=n,nrow=length(ri0))
	  R1<-matrix(
	    apply(
	      R,MARGIN=2,FUN=function(x){
	        quantile(x[x>=0],ppoints(ri1,a=1))}),ncol=n,nrow=length(ri1))
	  
	  
	  medianSim0<-apply(R0,1,median)
	  medianSim1<-apply(R1,1,median)
	  
	  lclSim0<-apply(R0,1,quantile,.025)
	  uclSim0<-apply(R0,1,quantile,.975)
	  lclSim1<-apply(R1,1,quantile,.025)
	  uclSim1<-apply(R1,1,quantile,.975)  	
	  
	  df<-data.frame(
	    Fitted_Residuals=c(ri0,ri1),
	    Simulated_Residuals=c(medianSim0,medianSim1),
	    UCLSim=c(uclSim0,uclSim1),LCLSim=c(lclSim0,lclSim1))
	  
	  g<-ggplot(df,aes(x=Simulated_Residuals,y=Fitted_Residuals))+
	    geom_point()+
	    geom_line(aes(x=Simulated_Residuals,y=LCLSim))+
	    geom_line(aes(x=Simulated_Residuals,y=UCLSim))+
	    geom_abline()
	  return(g)
	}
	