#_________________________________________________________
#
# Author: David J. Hearn
# Citation: 
# Use Policy: Free to use and modify provided citation is included in any resulting works
#
# Description: 
#	Wrapper function of R functions to identify and remove outliers
#
#FUNCTION:
#
# removeOutliers(data)
#	takes a data frame as input that, minimally, has the following fields: 
#		latitude, dayOfYear
#	provides the original data frame as output 
#		with rows removed that represent outliers
#
#_________________________________________________________



removeOutliers = function(data, DOY=TRUE, Latitude=TRUE)
{
	if(DOY)
	{
		outliers = boxplot(data$dayOfYear, plot=F)$out
		if(length(outliers)> 0) {
			data = data[-which(data$dayOfYear %in% outliers),]
		}
	}
	if(Latitude)
	{
		outliers = boxplot(data$latitude, plot=F)$out
		if(length(outliers) > 0) {
			data = data[-which(data$latitude %in% outliers),]
		}
	}
	return(data)
}

