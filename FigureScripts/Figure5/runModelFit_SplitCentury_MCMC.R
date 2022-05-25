
source("../SharedCode/theoryOfPhenology.txt.R")
source("../SharedCode/inferPhenology_BayesianMCMC.txt.R")
source("../SharedCode/simulatePhenology_edit.txt.R")
#source("../Weibull_PeaseEstimation.txt.r")
source("../SharedCode/removeOutliers.txt.R")
source("../SharedCode/biasCorrectionScripts.txt.R")

library(fitdistrplus)

run_BayesianReplicate = function(file, fileAddendum="../Data/", yearSplit=1950, min=0, max=365, onset="beta", duration="beta", n=365, firstHalf=T, nIter=100000)
{

	data = read.table(paste(fileAddendum,file,sep=""), header=T, sep='\t')
	data = removeOutliers(data)
	data = getDOY_LatitudeAdjusted(data)

	DOY = data$dayOfYear.LatitudeAdjusted
	DOY = DOY + mean(data$dayOfYear)

	if(firstHalf)
	{
		DOY = DOY[data$date<yearSplit]
	} else
	{
		DOY = DOY[data$date>=yearSplit]
	}

	#based on the observations, estimate the parameters for the f_X PDF
	#mlParams = getShapeParams.ML.estimate(data=DOY, min=min, max=max)
	#mlParams = list(oParList=list(shape1=59.7591442117949, shape2=136.075425886224), dParList=list( shape1=1.00344445554363, shape2=32.80174257825))

	mcmc.output.true = mcmc(
				data = DOY,
				onset = onset,
				duration = duration,
				#oParList= mlParams$oParList,
				#dParList= mlParams$dParList,
				min = min,
				max = max,
				logDPrior = logDPrior.uniformExpForBetaBeta,
				priorParList=list(mo=50,Mo=365,md=14,Md=365,mvo=1,Mvo=365*365,mvd=1,Mvd=365*365,lamb=3,min=0,max=365),
				#logDPrior = logDPrior.uniform,
				#priorParList= list(min=1, max=2000),
			        usePriorToInit = FALSE, 
			        useMLToInit=TRUE, 
			        #useMLToInit=FALSE, 
				useInputParams = FALSE,
				#useInputParams = TRUE,
				proposalRate = 0.0025,
				proposalAdapt = 0.0,
				printInterval=10,
				useBetaProposal=TRUE,
				nIter = nIter,
				graphResults = FALSE,
				#graphInterval = 10,
				n = n       #number of increments to use in numerical integration routines
	)

	printMCMCDataToFile(paste("../MCMCReplicates/",file,".firstHalf.",firstHalf,".txt",sep=""), mcmcOutput=mcmc.output.true, nIter=nIter, thin=2, type="betabeta",min=min, max=max)	
}






