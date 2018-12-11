###################################
###Caio Kenup######################
###Contact: caio.kenup@gmail.com###
###2018-12-11######################
###################################

#Function to generate predictions from GLMM models (and smooth them as well)

	smoothPredMer<-function(
	  #Model object
	  mod,
	  #Wether values are scaled
	  scaled=TRUE,
	  #List with mean and SD values for each covariable
	  scaleList=NULL,
	  #Number of simulations to run
	  nsims=1e3){
	  
	  #define which parameters are present in the model
		pars<-colnames(mod@frame)[-c(1,ncol(mod@frame))] 
		
		#define which is the random-effects variable
		random.var<-colnames(mod@frame)[ncol(mod@frame)]
		
		#If there are any explanatory variables
		if(length(pars)>0){
		  
		#define the response variable 
		resp.var<-colnames(mod@frame)[1] #
		
		#define the inverse of the link function of the model
		ilink<- mod@resp$family$linkinv 
		
		
		dat<-mod@frame #extract the underlying data
		dat<-dat[!is.na(dat[,resp.var]),] #remove NA entries from data
		
		dat[,resp.var]%<>%as.numeric #convert response variables to numeric data (T/F becomes 1/0)
		
		# Generate random sample of random random effects variable
		randomEffs<-sample(x = dat[,random.var],size = 100,replace=TRUE)
		
		###Simulating data
		dat.sim<-lapply(pars,function(x){ #for each parameters
			np<-pars[pars!=x] #define which are the other parameters, to be kept constant
			
		#Create temporary data
		temp.dat<-cbind(data.frame(
		    #100 values in the range of the current variable
				a = seq(min(dat[,x]), max(dat[,x]),length = 100)), 
			# bird,
			 #100 values with the mean of the fixed variables
		  	matrix(
		 	    apply(t(t(dat[,np])),2,mean),
				  nrow=100,ncol=length(np),byrow=TRUE))
			
			#Rename colnames accordingly
			colnames(temp.dat)<-c(x,np)
			
			#Add birds data
			temp.dat[,random.var]<-randomEffs
			
			####Return predicted values from bootstrap
			mySumm <- function(.) {
        predict(., newdata=temp.dat, re.form=NULL)
			}
			
      ####Collapse bootstrap results into median, 95% PI
      sumBoot <- function(merBoot) {
        return(
          data.frame(Fitted = ilink(apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.5, na.rm=TRUE)))),
                     Lower = ilink(apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE)))),
                     Upper = ilink(apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE))))
          )
        )
      }
			
      ##Bootstrap glmer model using simulated data
			boot1<-bootMer(mod, mySumm, nsim=nsims, use.u = FALSE, type="parametric",parallel = "multicore",ncpus = 4)
			##Summarise bootstrap
			PI.boot1<-sumBoot(boot1)
			##Added fited bootstrap data to temp.dat
			temp.dat%<>%cbind(PI.boot1)
			
			#Column names management
		  colnames(temp.dat)[1]<-'value'
			temp.dat$coef<-x
			
			#Return actual values, if the values were previously scaled
			if(scaled){
			  temp.dat$unscaled<-(temp.dat$value * scaleList$scale[x]) + scaleList$center[x]
			}else{
			temp.dat$unscaled<-temp.dat$value}
			
		return(temp.dat)		
		})
		
		#Name each data frame of simulated data
		names(dat.sim)<-pars
		
		#If values are scaled, transform them back to original values
		if(scaled){
		unsc<-sapply(par.cols,function(x){(dat[,x]*scaleList$scale[x])+scaleList$center[x]})
		colnames(unsc)<-par.cols
		dat[,par.cols]<-unsc}
		

		resu<-list()
		
		#For each set of simulated data (i.e., each covariable in the model)
		for (n in 1:length(dat.sim)){
		  
			#Base ggplot, to generate raw and smoothed plots
			gbase<-ggplot(dat,aes_string(
			  x=names(dat.sim)[n],
			  y=as.character(formula(mod))[2])) +
      scale_x_continuous(expand = c(0, 0)) + 
			geom_point(alpha=.4,size=3)
			
			
			#Raw plot, without smoothed means and intervals
			graw<-gbase +
			geom_ribbon(data = dat.sim[[n]],
			            aes(ymin = Lower, ymax = Upper, x = unscaled),
			                fill = 'black', alpha = .4, inherit.aes = FALSE) + 
			geom_line(data = dat.sim[[n]],color= 'black',aes(y = Fitted, x = unscaled))

			#Smoothed plot
			gsmooth0<-ggplot_build(gbase+
			    geom_smooth(data=dat.sim[[n]],
			      aes(x = unscaled,y = Fitted, colour = "mean"), method = "auto",se=F)+ 
			  geom_smooth(data=dat.sim[[n]],
			      aes(x = unscaled,y = Lower, colour = "ci"), method = "auto",se=F)+ 
			    geom_smooth(data=dat.sim[[n]],
			      aes(x = unscaled,y = Upper, colour = "ci"), method = "auto",se=F))
			# +

			#Creating a dataframe from the smoothed results
			smoothDF<-data.frame(x = gsmooth0$data[[2]]$x,
			            y    = gsmooth0$data[[2]]$y,
                  ymin = gsmooth0$data[[3]]$y,
                  ymax = gsmooth0$data[[4]]$y) 
			
			

			#Dataframe keeping only the explanatory and response variables
			dataTrue<-dat[,c(resp.var,names(dat.sim)[n])]
			colnames(dataTrue)<-c("X","Y")
			
			
			#Return, for each explanatory variable a list containg
			resu[[n]]<-list(plotRaw=graw, #the raw plot of the confidence intervals calculated
			                # plotSmooth=gsmooth,
			                respVar=resp.var,         #a character value expressing the response variable
			                expVar=names(dat.sim)[n], #a character value expressing the current explanatory variable
			                model=mod,                #the model used to generate predictions
			                dataSim=dat.sim[[n]],     #the simulated data for the current explanatory variable
			                dataTrue=dataTrue,        #the true, underlying data of the model
			                dataSmoothModel=smoothDF) #a data.frame containing the smoothed values for the predictions' means and confidence intervals
		}
		return(resu)}else{
		return(NULL)}
	}
	