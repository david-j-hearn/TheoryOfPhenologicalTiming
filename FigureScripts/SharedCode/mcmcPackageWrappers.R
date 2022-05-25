#Wrapper to Geyer's MCMC R package. Not used for manuscript. Use with caution! Results appear to be biased.

source("../SharedCode/inferPhenology_BayesianMCMC.txt.R")
source("../SharedCode/removeOutliers.txt.R")
source("../SharedCode/biasCorrectionScripts.txt.R")

if(!require(fitdistrplus)) 
{ 
	install.packages("fitdistrplus") 
}
if(!require(mcmc)) 
{ 
	install.packages("mcmc") 
}

library(mcmc)
library(fitdistrplus)


#general function to get an R function of the posterior distribution given data and a function of the prior with a list of its hyperparameters
#Prior can be any of the fillowing:
#logDPrior.normalForBetaBeta = function(oParList, dParList, priorParList = list(vmo = 25, mo = 100, vvo = 25, vo = 50, vmd = 25, md = 25, vvd = 25, vd = 50, min=0, max=365))
#logDPrior.uniformExpForBetaBeta = function(oParList, dParList, priorParList=list(mo=50,Mo=365,md=7,Md=365,mvo=1,Mvo=365*365,mvd=1,Mvd=365*365,lamb=1,min=0,max=365))
#logDPrior.uniform = function(oParList, dParList, priorParList= list(min=1,max=1000)) 
posteriorDistribution_FunctionGetter = function(data, priorFunction, priorHyperparameterList, min, max)
{

	count=0

	posteriorDistribution = function(params)
	{
		count <<- count + 1
		mean.onset = params[1]
		variance.onset = params[2]
		mean.duration = params[3]
		variance.duration = params[4]
		oParList = getShapeParameters(mean=mean.onset, variance=variance.onset, m=min, M=max)	
		dParList = getShapeParameters(mean=mean.duration, variance=variance.duration, m=0, M=max-mean.onset)	
		pr = calculateLogPrior(oParList=oParList, dParList=dParList, logDPrior = priorFunction, priorParList = priorHyperparameterList, verbose=FALSE)
		f_X = get.f_X(m=min, M=max, n=100, onset="beta", duration="beta", oParList=oParList, dParList=dParList)
		if(missing(f_X)) { stop("Couldn't generate f_X function. Quitting MCMC") }
		lk = calculateLogLikelihood(data=data, f_X=f_X)
		post = calculateLogPosterior.unnormalized( logPrior=pr, logLikelihood=lk )
		print(paste(c(count,params,pr,lk,post), collapse=" "))
		return(post)
	}
	return(posteriorDistribution)
}

#wrapper for mcmc package metrop function. Initializes parameter values to the maximum likelihood estimates
#the ... are the parameters for the metrop function that include:
#	nbatch 
#	blen 
#	nspac 
#	scale 
#	outfun 

runGeyerMCMC = function(dataFile=NA, data=NA, fileAddendum = "../Data/", priorFunction, priorHyperparameterList, min, max, initial=NA, ...)
{


	#get phenology data
	if(!is.na(dataFile))
	{
		message(paste("Retrieving data from ", fileAddendum, dataFile))
		data = read.table(paste(fileAddendum,dataFile,sep=""), header=T, sep='\t')
		data = removeOutliers(data)
		data = getDOY_LatitudeAdjusted(data)

		DOY = data$dayOfYear.LatitudeAdjusted
		data = DOY + mean(data$dayOfYear)
	}
	else if(length(data)<=1)
	{
		stop("No data file name and no data provided as input in runGeyerMCMC function. Quitting.")
	}

	#get initial parameters using maximum likelihood
	if(length(initial)<=1)
	{
		message("Running ML inference to obtain initial parameter values")
		MLParams = getShapeParams.ML.estimate(data, min=min, max=max)
		oParList = MLParams$oParList
		dParList = MLParams$dParList
		oMean = min + (max-min)*betaMean(oParList[[1]], oParList[[2]])
		oVar = ((max-min)^2) * betaVariance(oParList[[1]], oParList[[2]])
		dMean = (max-oMean)*betaMean(dParList[[1]], dParList[[2]])
		dVar = ((max-oMean)^2) * betaVariance(dParList[[1]], dParList[[2]])
		initial = c(oMean,oVar,dMean,dVar)
	}

	#get posterior distribution function
	message("Obtaining posterior distribution function")
	obj = posteriorDistribution_FunctionGetter(data=data, priorFunction = priorFunction, priorHyperparameterList = priorHyperparameterList, min=min, max=max)

	#run the MCMC
	message("running MCMC")
	mh.out = metrop(obj = obj, initial=initial, ...)
	return(list(mh = mh.out, meanX = mean(data), varX = var(data)))
}

