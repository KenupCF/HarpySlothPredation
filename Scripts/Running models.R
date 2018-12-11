###################################
###Caio Kenup######################
###Contact: caio.kenup@gmail.com###
###2018-12-11######################
###################################

###Importing Required Packages###
	require(lubridate)
	require(stringr)
	require(plyr)
	require(dplyr)
	require(magrittr)
	require(MASS)
	require(car)
	require(lme4)
	require(lmerTest)
  require(merTools)
	require(AICcmodavg)
	require(ResourceSelection)
	require(devtools)
	require(grid)
	require(gridExtra)
	require(ggplot2)
	require(rgdal)
	require(rgeos)
	
  #Set working directory to the project folder on your machine
	setwd()
	
	options(stringsAsFactors=FALSE)

	###Import custom functions
	source(".\\Functions\\Model diagnostics.R")
	source(".\\Functions\\smoothPredMer.R")
	
	source_url("https://raw.githubusercontent.com/KenupCF/Kenup-Repo/master/Quick%20Functions.R")

	###Importng Spatial Objects###
	
	##Project4String object
	ll.wgs84<-'+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
	
	##Points of interest
	interests<-readOGR(dsn=".\\GeoData",layer="Points of Interest",p4s=ll.wgs84)
	##Soberania National Park Limits
	snp<-readOGR(dsn=".\\GeoData",layer="Soberania NP")
	##Fixes of Harpy Eagles with prey
	# preylocs<-readOGR(dsn=".\\GeoData",layer="Prey Locations")
	
	##Reproject to WGS-84
	preylocs.ll<-spTransform(preylocs,interests@proj4string)
	snp.ll<-spTransform(snp,interests@proj4string)
	
	###
  ##Importing data
	###
	
	##Predation data
	dat<-read.csv(".\\Data\\Full_Data_Environmental_Variables.csv")
	
	##Climatic data
	clim.dat<-read.csv(".\\Data\\Climatic_Data.csv")

	###Defining explanatory variables
	  exp.vars<-c(
		  "NDVI_1",    #Normalized difference vegetation over a 1-km buffer of detection
		  "moon.yest", #Illuminated portion of lunar disc on night before detection
		  "moon.max",  #Maximum illuminated portion of lunar disc over 3-day period
		  "moon.min",  #Minimum illuminated portion of lunar disc over 3-day period
		  "Rainfall",  #Precipitation on day of detection
		  "rfall.sum", #Total rainfall over 3-day period
		  "Max_temp",  #Maximum temperature on day of detection
		  "Min_temp",  #Minimum temperature on day of detection
		  "min_temp.min", #Minimum temperature over 3-day period
		  "max_temp.max") #Maximum temperature over 3-day period
	
	
	###Z-transforming variables
	  
  	##Finding missing values in explanatory variables
  	missing.values<-apply(dat[,exp.vars],1,function(x){any(is.na(x))})
  	##Creating scale data.frame object
  	dat.scale<-scale(dat[!missing.values,exp.vars])
  	##Creating list of center and scale values
  	scaleList <- list(scale = attr(dat.scale, "scaled:scale"),
  		center = attr(dat.scale, "scaled:center"))
  	##Adding scaled variables in back to main data set
  	dat.scale%<>%cbind(dat[!missing.values,!colnames(dat)%in%exp.vars])
	
	###
  ###Model specification	
	###
	
  #For all sloths
	#Response variable is a function of intercept and:
  	formula.sloth.int<-paste("Sloth ~ 1 +", 
		"NDVI_1 + ",   
		"moon.max + ", 
		"NDVI_1 : moon.max + ",
		"rfall.sum + ", 
		"rfall.sum : Min_temp + ",
		"moon.yest + ",
		"Max_temp + ",
		"Rainfall",
		sep="")
	
	#For nocturnal animals
	#Response variable is a function of intercept and:
	formula.night.int<-paste("Night_prey ~ 1 +",
		"NDVI_1 + ",
		"moon.max + ",
		"(NDVI_1 : moon.max) + ",
		"moon.yest + ",
		"moon.min",sep="")
	
	#For Bradypus
	#Same model as for all Sloths, only changing the response variable
	formula.brad.int<-gsub(x=formula.sloth.int,"Sloth","is.brad")
	
	#For Choloepus
	#Same model as for all Sloths, only changing the response variable
	formula.chol.int<-gsub(x=formula.sloth.int,"Sloth","is.chol")

	
	###Exploratory analysis
	
	##Variance inflation analysis
	
	vif(glm(data=dat.scale,formula=as.formula(formula.night.int),family=binomial()))
	vif(glm(data=dat.scale,formula=as.formula(formula.sloth.int),family=binomial()))
	vif(glm(data=dat.scale,formula=as.formula(formula.brad.int),family=binomial()))
	vif(glm(data=dat.scale,formula=as.formula(formula.chol.int),family=binomial()))

	####
	###Running stepwise analysis
	####
	
	sw.analysis<-list()
	
	sw.analysis$Sloth<-MASS::stepAIC(
	  ##Object: Full model to be reduced
		object=glm(data=dat.scale,formula=as.formula(formula.sloth.int),family=binomial()),
		##Scope: List containing minimum and maximum models
		scope=list(upper=as.formula(formula.sloth.int),lower=as.formula("Sloth ~ 1")),
		##Direction: Start from the most complex model, reducing it
		direction='backward')

	sw.analysis$Night<-MASS::stepAIC(
	  ##Object: Full model to be reduced
		object=glm(data=dat.scale,formula=as.formula(formula.night.int),family=binomial()),
		##Scope: List containing minimum and maximum models
		scope=list(upper=as.formula(formula.night.int),lower=as.formula("Night_prey ~ 1")),
		##Direction: Start from the most complex model, reducing it
		direction='backward')
	
	sw.analysis$Bradypus<-MASS::stepAIC(
	  ##Object: Full model to be reduced
		object=glm(data=dat.scale,formula=as.formula(formula.brad.int),family=binomial()),
		##Scope: List containing minimum and maximum models
		scope=list(upper=as.formula(formula.sloth.int),lower=as.formula("Sloth ~ 1")),
		##Direction: Start from the most complex model, reducing it
		direction='backward')
	
	sw.analysis$Choloepus<-MASS::stepAIC(
	  ##Object: Full model to be reduced
		object=glm(data=dat.scale,formula=as.formula(formula.chol.int),family=binomial()),
		##Scope: List containing minimum and maximum models
		scope=list(upper=as.formula(formula.sloth.int),lower=as.formula("Sloth ~ 1")),
		##Direction: Start from the most complex model, reducing it
		direction='backward')
	
		###GOF tests for GLMs (Hoslem tests)
	
  	hoslem.test(
  	  sw.analysis$Sloth$data$Sloth[!is.na(sw.analysis$Sloth$data$Sloth)], 
  	  fitted(sw.analysis$Sloth))
  	
  	hoslem.test(
  	  sw.analysis$Night$data$Night[!is.na(sw.analysis$Night$data$Night)], 
  	  fitted(sw.analysis$Night))
  	
  	hoslem.test(
  	  sw.analysis$Choloepus$data$is.chol[!is.na(sw.analysis$Choloepus$data$is.chol)],
  	  fitted(sw.analysis$Choloepus))
  	
  	hoslem.test(
  	  sw.analysis$Bradypus$data$is.brad[!is.na(sw.analysis$Bradypus$data$is.brad)],
  	  fitted(sw.analysis$Bradypus))
	
	
	
	####
	###Mixed model analysis
	####
  	
	##Creating mixed-model formulas from stepwise results
	glmm.forms<-lapply(sw.analysis,function(mod){
		form<-as.character(mod$formula)
		##Add Bird effect on the intercept
		form<-paste(form[2],form[1],form[3:length(form)],"+(1|Bird)") 
		return(form)
		})
	
  ##Run generalized linear mixed models using new formulas
	glmm.analysis<-lapply(1:length(sw.analysis),function(m){
		glmer(data=dat.scale,formula=as.formula(glmm.forms[[m]]),family=binomial())
		})
	
  ##Keep model names from stepwise analyses
	names(glmm.analysis)<-names(sw.analysis)
	

	###GLMM diagnostics
	##Run diagnostics for all glmms
	D<-lapply(glmm.analysis,diagnostics.glmm,n=1e1) 
	names(D)<-names(glmm.analysis)

	##Generating (smoothed) predictions from GLMMs
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
	
	
	#Generated smoothed prediction data for all models
	system.time({sloth.plots<-smoothPredMer(mod = glmm.analysis[[1]],scaled = TRUE,scaleList = scaleList,nsims=1e3)})
	system.time({night.plots<-smoothPredMer(mod = glmm.analysis[[2]],scaled = TRUE,scaleList = scaleList,nsims=1e3)})
	system.time({brad.plots<-smoothPredMer(mod = glmm.analysis[[3]],scaled = TRUE,scaleList = scaleList,nsims=1e3)})
	system.time({chol.plots<-smoothPredMer(mod = glmm.analysis[[4]],scaled = TRUE,scaleList = scaleList,nsims=1e3)})
	 
	  
	#Generating summary table for all models 
	glmm.table<-rbind.fill(
	  lapply(names(glmm.analysis),function(m){
	  mod<-glmm.analysis[[m]]
	  x<-(summary(mod)$coefficients)
	  y<-data.frame(Model=m,
	                Var=rownames(x),
	                Beta=x[,'Estimate'],
	                SE=x[,'Std. Error'],
	                p=x[,"Pr(>|z|)"])
	  return(y)}))%>%
	  filter(Var!="(Intercept)")
	
	
	######################	
	#Plotting predictions#
	######################
	
	all.preds<-c(sloth.plots,night.plots,chol.plots,brad.plots)
	
	#Extracting X and Y labels
	xlabs<-sapply(all.preds,function(x){
	  x$expVar
	  })
	ylabs<-sapply(all.preds,function(x){
	  x$respVar
	  })
	
	#Converting labels to more verbose versions
	ylabs[ylabs=="Sloth"]<-"Sloth predation probability"
	ylabs[ylabs=="is.chol"]<-"Two-toed sloth predation probability"
	ylabs[ylabs=="Night_prey"]<-"Nocturnal animal predation probability"
	ylabs[ylabs=="is.brad"]<-"Three-toed sloth predation probability"
		
	xlabs[xlabs=="Min_temp"]<-"Minimum temperature (ºC)"
	xlabs[xlabs=="moon.max"]<-"Max. illuminated portion of lunar disc (3-day, %)"
	xlabs[xlabs=="moon.min"]<-"Min. illuminated portion of lunar disc (3-day, %)"
	xlabs[xlabs=="moon.yest"]<-"Illuminated portion of lunar disc (%)"
	
	###Plotting parameters
	cex.<-.67
	
	#Generating one (smoothed) plot for each explanatory variable of each model
	plotsSmooth<-lapply(1:length(all.plots),function(n){
	    
	    #Get prediction object
	    x<-all.plots[[n]]
	    
	    #Create plot
	    ggplot() +
	    #Plot jittered actual data
	  	geom_jitter(data = x$dataTrue,
	                mapping = aes(x=Y,y=X),
	                height = 0,width=.4,alpha=.4,size=2*cex.)+
	    #Plot smoothed trend line
	    geom_line(data=x$dataSmoothModel,
	              mapping = aes(x=x,y=y))+
	    #Plot smoothed confidence intervals
	    geom_ribbon(data=x$dataSmoothModel,
	                mapping = aes(x=x,ymin=ymin,ymax=ymax),alpha=.4)+
	    #Label axes accordingly
	    xlab(xlabs[n]) +
			ylab(ylabs[n]) +
	    #Aesthethic changes to the plot
			theme(
			  text = element_text(size=10*cex.),
			  panel.grid.major = element_blank(), 
			  panel.grid.minor = element_blank(),
				panel.background = element_blank(), 
				axis.line = element_line(colour = "black")) + 
			scale_x_continuous(
			  # breaks=seq(from=18,to=30,by=1),
			  expand=c(0,0))
	  })

	#Generating one (raw) plot for each explanatory variable of each model
	plotsRaw<-lapply(1:length(all.plots),function(n){
	    
	    #Get prediction object
	    x<-all.plots[[n]]
	    
	    #Create plot
	    ggplot() +
	    #Plot jittered actual data
	  	geom_jitter(data = x$dataTrue,
	                mapping = aes(x=Y,y=X),
	                height = 0,width=.4,alpha=.4,size=2*cex.)+
	    #Plot smoothed trend line
	    geom_line(data=x$dataSim,
	              mapping = aes(y = Fitted, x = unscaled))+
	    #Plot smoothed confidence intervals
	    geom_ribbon(data=x$dataSim,
	                mapping = aes(ymin = Lower, ymax = Upper, x = unscaled),alpha=.4)+
	    #Label axes accordingly
	    xlab(xlabs[n]) +
			ylab(ylabs[n]) +
	    #Aesthethic changes to the plot
			theme(
			  text = element_text(size=10*cex.),
			  panel.grid.major = element_blank(), 
			  panel.grid.minor = element_blank(),
				panel.background = element_blank(), 
				axis.line = element_line(colour = "black")) + 
			scale_x_continuous(
			  # breaks=seq(from=18,to=30,by=1),
			  expand=c(0,0))
	  })
	
	
	##Exporting TIFF images
	set<-2250/3
	# svg(
	tiff(
	   file='.\\Images\\Figure 1_2_Smooth.tiff',
	   # units='cm',
	  res=300,
	  pointsize=12,
	  width=set*3,height=set*2)
	
	grid.arrange(ggplotGrob(plotsSmooth[[2]]), ggplotGrob(plotsSmooth[[1]]), ggplotGrob(plotsSmooth[[6]]), 
	             ggplotGrob(plotsSmooth[[5]]), ggplotGrob(plotsSmooth[[4]]), ggplotGrob(plotsSmooth[[3]]), 
	             ncol=3,nrow=2)
	dev.off()
	
	tiff(
	   file='.\\Images\\Figure 1_2_Raw.tiff',
	   # units='cm',
	  res=300,
	  pointsize=12,
	  width=set*3,height=set*2)
	
	grid.arrange(ggplotGrob(plotsRaw[[2]]), ggplotGrob(plotsRaw[[1]]), ggplotGrob(plotsRaw[[6]]), 
	             ggplotGrob(plotsRaw[[5]]), ggplotGrob(plotsRaw[[4]]), ggplotGrob(plotsRaw[[3]]), 
	             ncol=3,nrow=2)
	dev.off()

	