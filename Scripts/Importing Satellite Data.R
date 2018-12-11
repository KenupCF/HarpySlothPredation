###################################
###Caio Kenup######################
###Contact: caio.kenup@gmail.com###
###2018-12-11######################
###################################
require(plyr)
require(dplyr)
require(magrittr)
require(stringr)
require(raster)
require(rgeos)
require(rgdal)
require(lubridate)
require(devtools)
require(RStoolbox)
options(stringsAsFactors=FALSE)

#Import custom functions
source_url("https://raw.githubusercontent.com/KenupCF/Kenup-Repo/master/Quick%20Functions.R")

#Set working directory to the project folder on your machine
setwd()

#Define UTM projection
p4s.utm.17n<-'+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs +ellps=WGS84'

#Import waterbodies shapefile
water<-readOGR(dsn='.\\GeoData',layer='Water_Bodies')

#Import Soberania National Park files
pnsob<-readOGR(dsn='.\\GeoData',layer='Soberania NP')

#Import prey location fixes
# prey<-read.csv(".\\Data\\Data 200 Presas\\Prey Events (latlong, bird, date).csv")
prey<-read.csv(".\\Data\\Full_Data_Environmental_Variables.csv")

prey$Date<-paste(prey$Year,zero_pad(prey$Month,2),zero_pad(prey$Day,2),sep="-")
prey$ID<-1:nrow(prey)

#Subset only to detectoins with valid coordinates
prey.coord<-prey[!is.na(prey$UTM.E),] 

#Convert to SpatialPointsDataFrame
prey.sp<-SpatialPointsDataFrame(
	data=prey.coord,coords=prey.coord[,c('UTM.E','UTM.N')],proj4string=CRS(p4s.utm.17n))

prey.sp@data$ID<-1:nrow(prey.sp@data)

#Create a 1-,5- and 10-km buffer around each location
prey.sp.01buff<-gBuffer(prey.sp,
                         byid=TRUE,id=prey.sp@data$ID,width=1e3)
prey.sp.05buff<-gBuffer(prey.sp,
                         byid=TRUE,id=prey.sp@data$ID,width=5e3)
prey.sp.10buff<-gBuffer(prey.sp,width=10e3)


##Arrange all satellite images

##Because the satellite image files are too large, we are unable to host them on GitHub.
##However, they can easily be downloaded at https://earthexplorer.usgs.gov/

all.bands<-list(
	B6.2=list.files(recursive=TRUE,pattern='_B6_VCID_2.TIF'),
	B1=list.files(recursive=TRUE,pattern='_B1.TIF'),
	B2=list.files(recursive=TRUE,pattern='_B2.TIF'),
	B3=list.files(recursive=TRUE,pattern='_B3.TIF'),
	B4=list.files(recursive=TRUE,pattern='_B4.TIF'),
	B5=list.files(recursive=TRUE,pattern='_B5.TIF')
	)
	

##Create a single data.frame with all satelite image file names
all.bands<-rbind.fill(lapply(all.bands,function(x){
		x<-x[!str_detect(x,'.gz')] #remove zipped files from list
		y<-gsub(x=x,"^.*/","")		#remove path, keeping only file name
		y<-gsub(x=y,'B6_VCID_2','B6.2') #Renamed Band 6 objects
		z<-str_split_fixed(y,'_',14)	#split file name creating a information matrix 
		empty.col<-which(apply(z,2,function(x){all(x=='')}))	#determine which matrix columns are empty
		z<-z[,-empty.col] #remove empty columns from matrix
		empty.row<-which(apply(z,1,function(x){any(x=='')})) #determine which rows have a trailing column
		z[empty.row,ncol(z)]<-z[empty.row,ncol(z)-1] #repeat second-last column to last column on such rows
		colnames(z)[c(4,ncol(z)-1,ncol(z))]<-c('date','gap_mask','band') #rename important columns
		z<-cbind(z,path=x) #bind file path to matrix
		z<-z[,!is.na(colnames(z))] #keep only named columns on matrix
		z<-as.data.frame(z) #convert to data.frame
		z$gap_mask<-z$gap_mask=='GM' #convert 'gap_mask' to a logical column
		z$band<-gsub(x=z$band,'.TIF','') #remove file extension from 'band' column
		z$path<-as.character(z$path)
		return(z)
	}))

#Arrange them by date and band
all.bands%<>%arrange(date,gap_mask,band)
	
#Count how many bands are available on each date
all.bands.summ<-all.bands%>%
	dplyr::group_by(date)%>%
	dplyr::summarise(n=sum(!gap_mask),n_gm=sum(gap_mask))

#Detect dates with incomplete band information	
incomplete<-all.bands.summ%>%
	filter(n!=6 | !n_gm %in%c(0,6) )%>%
	pull(date)

#Only keep dates with complete band sets
all.bands%<>%filter(!date%in%incomplete)


#Split bands by date and wether or not they have a gap mask	
all.bands.split<-split(all.bands,f=list(all.bands$date,all.bands$gap_mask))

#Detect which splices are from gap masls
gm<-str_detect(names(all.bands.split),'TRUE')

#Separated gap-masks from non-gap masks
all.bands.split.gm<-all.bands.split[gm]
all.bands.split<-all.bands.split[!gm]

#Renaming splices
names(all.bands.split)<-gsub(x=names(all.bands.split),'.FALSE','')
names(all.bands.split.gm)<-gsub(x=names(all.bands.split),'.TRUE','')

#Empty list to be filled in for loop (for satellite images)
sat.list<-list()
#Empty list to be filled in for loop (for gap mask (GM) images)
gm.list<-list()

###Importing and stacking all rasters 


###SUBSETTING FOR TESTING
all.bands.split<-all.bands.split[sample(1:67,size = 12)]

#For each satellite date
for(r in 1:length(all.bands.split)){
  
  #GM and non-GM files
	temp<-all.bands.split[[r]]
	temp.gm<-all.bands.split.gm[[names(all.bands.split)[r]]]
	
	#Import and rename non-GM bands
	sat.list[[r]]<-stack(lapply(temp$path,function(x){x}))
	names(sat.list[[r]])<-gsub(x=names(sat.list[[r]]),"^.*_B","B")
	names(sat.list[[r]])<-gsub(x=names(sat.list[[r]]),"B6_VCID_2","B6.2")
	
	#if there are any GM bands
	if(nrow(temp.gm)>0){
	  
	 	#Import and rename GM bands
		gm.list[[r]]<-stack(lapply(temp.gm$path,function(x){x}))
		names(gm.list[[r]])<-gsub(x=names(gm.list[[r]]),"^.*_B","Bmk")
		names(gm.list[[r]])<-gsub(x=names(gm.list[[r]]),"6_VCID_2","6.2")
		
		#Add GM bands to non-GM stack
		sat.list[[r]]<-stack(sat.list[[r]],gm.list[[r]])
	}
	
	#For the first date
	if(r==1){
		#Rasterize and crop water bodies
	  water.r<-rasterize(x=water,y=sat.list[[1]],field=1)
		water.r<-crop(water.r,extent(prey.sp.10buff))
		names(water.r)<-'Water'
	}
	
	#Crop stack to the extent of the buffer of prey locations (discarding unncessary pixels)
	sat.list[[r]]<-crop(sat.list[[r]],extent(prey.sp.10buff))
	
	#Add water-mask band to stack
	sat.list[[r]]<-stack(sat.list[[r]],water.r)
	
	#Create a cloud-mask using cloud detection from RStoolbox
	cloudTemp<-RStoolbox::cloudMask(sat.list[[r]],blue='B1',tir='B6.2',plot=FALSE)
	
	#Add cloud-mask to stack
	sat.list[[r]]<-stack(sat.list[[r]],cloudTemp)
	
	print(r)
}

names(sat.list)<-names(all.bands.split)

datsats<-as.Date(names(sat.list),format="%Y%m%d")

###SUBSETTING FOR TESTING
prey.sp<-prey.sp[sample(x = 1:nrow(prey.sp),size = 6),]

###Calculating NDVI values for each predation location

	sat.values<-list()
	dat.list<-list()
	
	for(loc in 1:nrow(prey.sp)){
	  
	  #extracting 1 and 5 buffers of location
	  buf1k<-prey.sp.01buff[loc,]
		buf5k<-prey.sp.05buff[loc,]
		
	  sat.values.5k<-list()
		sat.values.1k<-list()
		#Get location date
		datLoop<-as.Date(prey.sp.01buff@data[loc,'Date'])
		#Difference (in days) between location date and satelite images
		datdiff<-as.numeric(difftime(datLoop,datsats,units='days'))
		#Which are the 5 closest satellite dates from the location dates
	    d<-which(abs(datdiff) %in% sort(abs(datdiff))[1:5])
	    dat.list[[loc]]<-datsats[d]
		  names(dat.list)[loc]<-as.character(datLoop)
		
		#For each of these 5 dates
  		for(s in d){
  		  #Extract all pixels in a 1km and 5km buffers
  			sat.values.1k[[s]]<-extract(sat.list[[s]],buf1k)[[1]]
  			sat.values.5k[[s]]<-extract(sat.list[[s]],buf5k)[[1]]
  			cols<-colnames(sat.values.1k[[s]])	
  			
  			#If there aren't any gap mask bands, create dummy gap mask bands
    			if(! any(str_detect(cols,'mk')) ){
    				sat.values.5k[[s]]%<>%cbind(Bmk3=1)
    				sat.values.5k[[s]]%<>%cbind(Bmk4=1)
    				sat.values.1k[[s]]%<>%cbind(Bmk3=1)
    				sat.values.1k[[s]]%<>%cbind(Bmk4=1)
    			}
  			
  			#Convert extracted values to data.frame, add date variable
  			sat.values.5k[[s]]%<>%as.data.frame%>%
  				mutate(date_sat=datsats[s])
  			sat.values.1k[[s]]%<>%as.data.frame%>%
  				mutate(date_sat=datsats[s])
  			
  			print(s)
  			}
		
		#Aggregated extracted values   
		sat.values.1k%<>%rbind.fill%>%
				mutate(buffer.size.km=1)
		sat.values.5k%<>%rbind.fill%>%
				mutate(buffer.size.km=5)
    sat.values[[loc]]<-rbind.fill(sat.values.1k,sat.values.5k)	
    
    sat.values[[loc]]%<>%
        #Removed invalid items (Clouds, Water, and Masked values)
				dplyr::filter(!is.na(CMASK),Bmk4!=0,Bmk3!=0)%>%
				dplyr::filter(is.na(Water))%>%
        #Calculate NDVI of each pixel
				dplyr::mutate(
					NDVI=((B4-B3)/(B3+B4)))%>%
        #Remove NAs
				dplyr::filter(!is.na(NDVI))%>%
        #Add Date and Prey information to all pixels
				dplyr::mutate(Date=prey.sp.01buff@data$Date[loc],
				  ID=prey.sp.01buff@data$ID[loc])	
	
		print(paste0('loc = ',loc))
  
	}
	
	##Summarise NDVI values for 1 and 5 km buffers in a single data.frame
	##(for each Harpy Eagle GPS fix)
	sat.values.df<-rbind.fill(sat.values)%>%
	  dplyr::group_by(ID,Date,date_sat)%>%
	  dplyr::summarise(NDVI_1=mean(NDVI[buffer.size.km==1]),
	            NDVI_5=mean(NDVI[buffer.size.km==5]))%>%
	  dplyr::group_by(ID,Date)%>%
	  dplyr::mutate(diffDate=as.numeric(difftime(time1 = Date,time2 = date_sat,units = "days")))
	  
	#Split data.frame by each prey location
	temp<-split(sat.values.df,sat.values.df$ID)
	
	
	#For each ID
	sat.values.df<-lapply(temp,function(df){
	  datesats<-paste(df$date_sat,collapse=", ") #create a vector of all date of satellites used
	  
	  #If prey was located in a day that has an available satellite image, we used the NDVIs for that satelite image
	  if(any(df$diffDate==0)){
	    df<-df[which(df$diffDate==0),]%>%
	      mutate(date_sat=Date,diffDate=NULL)
	  #Otherwise, we linearly interpolate the NDVI values from the 5 closest dates
	  }else{
	    ndvi1k<-approx(x=df$diffDate,y=df$NDVI_1,xout=0)$y
	    ndvi5k<-approx(x=df$diffDate,y=df$NDVI_5,xout=0)$y
	    df<-df[1,]%>%
	      mutate(NDVI_1=ndvi1k,NDVI_5=ndvi5k,date_sat=datesats,diffDate=NULL)
	    }
	  return(df)
	  })%>%
	  rbind.fill() #bind all data.frames together again
	
  #Add these NDVI values to main dataset
	prey<-merge(prey,sat.values.df,by="ID")
	
	










	

