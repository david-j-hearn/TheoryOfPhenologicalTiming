#_________________________________________________________
#
# Author: David J. Hearn
# Citation: 
# Use Policy: Free to use and modify provided citation is included in any resulting works
#
# Description: 
#	Functions to transform data to correct for biases in:
#		latitude
#		sample size 
#	 All functions require a data frame as input
#		The data frame needs the following named columns: 
#			ID, species, latitude, date, dayOfYear
#		There should only be a single species in the data frame
#	All functions return a data frame with transformed data corrected for biases
#
#FUNCTIONS:
#
# getDOY_LatitudeAdjusted
#	corrects day of year for latitudinal sampling biases
#	returns a data frame with the following fields:
#		ID, species, date, latitude, dayOfYear, dayOfYear.LatitudeAdjusted
#
# getMinDOY_LatitudeSSAdjusted
#	produces the minimum day of year for each input year
#	corrects for latitudinal and sampling size biases
#	returns a data frame with the following fields:
#		year, 
#		dayOfYear.minimum,
#		dayOfYear.minimum.LatitudeAdjusted,
#		dayOfYear.minimum.SampleSizeAdjusted,
#		dayOfYear.minimum.LatitudeAdjusted.SampleSizeAdjusted
#
#_________________________________________________________



#_________________________ADJUST DOY FOR LATITUDINAL BIAS__________________

getDOY_LatitudeAdjusted = function(data) {

dayOfYear = data$dayOfYear
latitude = data$latitude
date = data$date
ID = data$ID
species = data$species

Lat_DOY = lm(dayOfYear ~ latitude)
dayOfYear.LatitudeAdjusted = dayOfYear - (Lat_DOY$coefficients[[2]] * latitude + Lat_DOY$coefficients[[1]])

return( data.frame( ID, species, date, latitude, dayOfYear, dayOfYear.LatitudeAdjusted))

}

#__MINIMUM DOY ADJUSTED FOR LATITUDINAL AND SAMPLE SIZE BIASES______________

getMinDOY_LatitudeSSAdjusted = function(data) {

doy = data$dayOfYear
date = data$date
lat = data$latitude

dateFloor = floor(date)
yearList = c(min(dateFloor):max(dateFloor))


#No adjustments

minDOYs = sapply(yearList, returnMinDOY, years=dateFloor, DOYs=doy)
tb = !is.na(minDOYs)
minDOYs= minDOYs[tb]
yearList0 = yearList[tb]


#Adjustments for sample size

mu = mean(doy)
stdev = sd(doy)

minDOYsT0 = sapply(yearList,returnTransformedMinDOY, years=dateFloor, DOYs=doy, mu=mu, sd=stdev)
tb = !is.na(minDOYsT0)
minDOYsT0= minDOYsT0[tb]
yearList0.1 = yearList[tb]

#Adjustments for latitude 

Lat_DOY = lm(doy ~ lat)
doyT1 = doy - (Lat_DOY$coefficients[[2]] * lat + Lat_DOY$coefficients[[1]])

minDOYsT1 = sapply(yearList,returnMinDOY, years=dateFloor, DOYs=doyT1)
tb = !is.na(minDOYsT1)
minDOYsT1= minDOYsT1[tb]
yearList1 = yearList[tb]

#Adjustments for latitude and sample size

mu = mean(doyT1)
stdev = sd(doyT1)

minDOYsT2 = sapply(yearList,returnTransformedMinDOY, years=dateFloor, DOYs=doyT1 ,mu=mu, sd=stdev)
tb = !is.na(minDOYsT2)
minDOYsT2 = minDOYsT2[tb]
yearList2 = yearList[tb]

if(length(yearList1) != length(yearList2)) {
      print("year lists 1 and 2 don't match in length. Quitting.")
      return(NA)
	}
if(length(yearList0) != length(yearList1)) {
      print("year lists 0 and 1 don't match in length. Quitting.")
      return(NA)
	}
if(length(yearList0) != length(yearList0.1)) {
      print("year lists 0 and 0.1 don't match in length. Quitting.")
      return(NA)
	}

year = yearList0
dayOfYear.minimum = minDOYs
dayOfYear.minimum.LatitudeAdjusted = minDOYsT1
dayOfYear.minimum.SampleSizeAdjusted = minDOYsT0
dayOfYear.minimum.LatitudeAdjusted.SampleSizeAdjusted = minDOYsT2

return(data.frame(year, dayOfYear.minimum, dayOfYear.minimum.LatitudeAdjusted, dayOfYear.minimum.SampleSizeAdjusted, dayOfYear.minimum.LatitudeAdjusted.SampleSizeAdjusted))

}

