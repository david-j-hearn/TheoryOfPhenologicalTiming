source("../SharedCode/helperFunctions.txt.R")


#gen	onset.shape1	onset.shape2	duration.shape1	duration.shape2	lk	prior	post	mean.onset	variance.onset	mean.duration	variance.duration			




#Assumes the first generation uses maximum likelihood parameters
create_PDFGraph_PosteriorAverage_beta_beta = function(fileData, fileMCMCFirst, fileMCMCSecond, m=0, M=365, yearCutoff=1950, n=100, outFile=NA, thinInterval=10, burnin=300, textPos=0.9, minX=50, maxX=200, lineColor="yellow", durationFill=rgb(0,0,0,0.1), dx=0.1, minD=0, maxD=60)
{
	print(paste("Working on ", fileData))
	d = read.table(fileData, header=T, sep='\t')
	d = removeOutliers(d)
	d = getDOY_LatitudeAdjusted(d)
	DOY = d$dayOfYear.LatitudeAdjusted + mean(d$dayOfYear)
	DOY1 = DOY[d$date<yearCutoff]
	DOY2 = DOY[d$date>=yearCutoff]

	d1 = read.table(fileMCMCFirst, header=T, sep='\t')
	d2 = read.table(fileMCMCSecond, header=T, sep='\t')

	ML.Onset = list(shape1=d1$onset.shape1[1], shape2=d1$onset.shape2[1])
	ML.Duration = list(shape1=d1$duration.shape1[1], shape2= d2$duration.shape2[1])

	print(ML.Onset)
	print(ML.Duration)

	xv = seq(from=m,to=M,by=dx)
	len = length(xv)
	fX.ML = get.f_X(m=m, M=M, n=n, onset="beta", duration="beta", oParList=ML.Onset, dParList=ML.Duration)
	dx1 = (M-m) / n

	#get the average curve - early
	d1 = thin(MCMCData = d1, burnin=burnin, thinInterval=thinInterval)
	print(paste("There are ", nrow(d1), "posterior samples to process from the earlier years."))
	yX.early = rep(0, length(xv))
	yOnset.early = rep(0, length(xv))
	yDuration.early = rep(0, length(xv))
	yCessation.early = rep(0, (n+1))
	meanO.early = 0
	meanD.early = 0
	meanC.early = 0
	for(i in 1:nrow(d1))
	{
		print(paste("Early", i, fileData))
		oS1 = d1$onset.shape1[i]
		oS2 = d1$onset.shape2[i]
		dS1 = d1$duration.shape1[i]
		dS2 = d1$duration.shape2[i]
		oParList = list(shape1 = oS1, shape2 = oS2)
		dParList = list(shape1 = dS1, shape2 = dS2)

		fs = get.f_X.f_C(m=m, M=M, n=n, onset="beta", duration="beta", oParList=oParList, dParList=dParList)
		fX = fs$f_X
		fCs.early = fs$fCs
		yX.early = yX.early + fX(xv)
		yOnset.early = yOnset.early + dbeta((xv-m)/(M-m), shape1=oS1, shape2=oS2)/(M-m)
		
		meanO.early.t =  betaMean(shape1=oS1, shape2=oS2)*(M-m)+m
		meanO.early = meanO.early + meanO.early.t
		yCessation.early = yCessation.early + fCs.early$y
		meanC.early = meanC.early+(fCs.early$x %*% (fCs.early$y*dx1)) - fCs.early$x[1]*fCs.early$y[1]*dx1/2 - fCs.early$x[length(fCs.early$x)]*fCs.early$y[length(fCs.early$x)]*dx1/2

		f_Dy = f_D(xv, m=m, M=M, onset="beta", duration="beta", oParList=oParList, dParList=dParList, n=n)
		yDuration.early = yDuration.early + f_Dy
		meanD = (M-meanO.early.t)*betaMean(shape1=dS1, shape2=dS2)
		meanD.early = meanD.early + meanD
	}
	yX.early = yX.early / nrow(d1)
	yOnset.early = yOnset.early / nrow(d1)
	yDuration.early = yDuration.early / nrow(d1)
	yCessation.early = yCessation.early / nrow(d1) 
	meanO.early = meanO.early / nrow(d1)
	meanD.early = meanD.early / nrow(d1)
	meanC.early = meanC.early / nrow(d1)

	#get the average curve - late
	d2 = thin(MCMCData = d2, burnin=burnin, thinInterval=thinInterval)
	print(paste("There are ", nrow(d2), "posterior samples to process from the more recent years."))
	yX.late = rep(0, length(xv))
	yOnset.late = rep(0, length(xv))
	yDuration.late = rep(0, length(xv))
	yCessation.late = rep(0, (n+1))
	meanO.late = 0
	meanD.late = 0
	meanC.late = 0
	for(i in 1:nrow(d2))
	{

		print(paste("Late", i, fileData))
		oS1 = d2$onset.shape1[i]
		oS2 = d2$onset.shape2[i]
		dS1 = d2$duration.shape1[i]
		dS2 = d2$duration.shape2[i]
		oParList = list(shape1 = oS1, shape2 = oS2)
		dParList = list(shape1 = dS1, shape2 = dS2)

		fs = get.f_X.f_C(m=m, M=M, n=n, onset="beta", duration="beta", oParList=oParList, dParList=dParList)
		fX = fs$f_X
		fCs.late = fs$fCs
		yX.late = yX.late + fX(xv)
		yOnset.late = yOnset.late + dbeta((xv-m)/(M-m), shape1=oS1, shape2=oS2)/(M-m)
		meanO.late.t = betaMean(shape1=oS1, shape2=oS2)*(M-m)+m
		meanO.late = meanO.late + meanO.late.t
		yCessation.late = yCessation.late + fCs.late$y
		meanC.late = meanC.late + (fCs.late$x %*% (fCs.late$y*dx1)) - fCs.late$x[1]*fCs.late$y[1]*dx1/2 - fCs.late$x[length(fCs.late$x)]*fCs.late$y[length(fCs.late$x)]*dx1/2

		f_Dy = f_D(xv, m=m, M=M, onset="beta", duration="beta", oParList=oParList, dParList=dParList, n=n)
		yDuration.late = yDuration.late + f_Dy
		#meanD = calcDMean(m=m, M=M, oS1=oS1, oS2=oS2, dS1=dS1, dS2=dS2,n=n)
		meanD = (M-meanO.late.t)*betaMean(shape1=dS1, shape2=dS2)
		meanD.late = meanD.late + meanD

	}
	yX.late = yX.late / nrow(d2)
	yOnset.late = yOnset.late / nrow(d2)
	yDuration.late = yDuration.late / nrow(d2)
	yCessation.late = yCessation.late / nrow(d2) 
	meanO.late = meanO.late / nrow(d2)
	meanD.late = meanD.late / nrow(d2)
	meanC.late = meanC.late / nrow(d2)

	#duration calculations
	d.x.lim = c(minD,maxD)

	if(!is.na(outFile)) pdf(outFile)

	layout(matrix(c(1,1,1,1,2,2,3,3,4,5,6,7), nrow=3, ncol=4, byrow=T))

	#ML all data
	y=fX.ML(xv)
	hd = hist(DOY, xlim=c(minX,maxX), probability=T, plot=F, breaks=seq(minX,maxX,by=5),col=rgb(1,0,1,0.1)) 
	hist(DOY, xlim=c(minX,maxX), ylim=c(0,max(y,hd$density)), probability=T, main=NULL, xlab="Observed DOY", breaks=seq(minX,maxX,by=5), col=rgb(1,0,1,0.1)) 
	#points(xv, y, type='l',col=rgb(1,0,1,1))
	abline(v=mean(DOY), col=lineColor)
	text(mean(DOY), textPos*max(fX.ML(xv)), paste(round(mean(DOY),1)), pos=4)

	#MCMC average early
	hd = hist(DOY1, xlim=c(minX,maxX), probability=T, plot=F, breaks=seq(minX,maxX,by=5)) 
	hist(DOY1, xlim=c(minX,maxX), ylim=c(0, max(yX.early,hd$density)), probability=T, main=NULL, xlab=paste("Observed DOY ( <",yearCutoff, ")"), ylab="Density", breaks=seq(minX,maxX,by=5),col=rgb(1,0,1,0.1)) 
	points(xv, yX.early, type='l',col=rgb(1,0,1,1))
	abline(v=mean(DOY1), col=lineColor)
	text(mean(DOY1), textPos*max(yX.early), paste(round(mean(DOY1),1)), pos=4)

	#MCMC average late
	hd = hist(DOY2, xlim=c(minX,maxX), probability=T, breaks=seq(minX,maxX,by=5), plot=F) 
	hist(DOY2, xlim=c(minX,maxX), ylim = c(0, max(yX.late, hd$density)), probability=T, main=NULL, xlab=paste("Observed DOY ( >=", yearCutoff, ")"), ylab="Density", breaks=seq(minX,maxX,by=5),col=rgb(1,0,1,0.1)) 
	points(xv, yX.late, type='l',col=rgb(1,0,1,1))
	abline(v=mean(DOY2), col=lineColor)
	text(mean(DOY2), textPos*max(yX.late), paste(round(mean(DOY2),1)), pos=4)

	#Onset and Cessation early
	plot(xv, yOnset.early, type="l", ylab="Density", xlab="Onset DOY (early years)", xlim=c(minX,maxX),col=rgb(1,0,0,1))
	polygon(xv, yOnset.early, col=rgb(1,0,0,0.1), border=NA)
	points(fCs.early$x,yCessation.early, type='l',col=rgb(0,0,1,1))
	polygon(fCs.early$x, yCessation.early, col=rgb(0,0,1,0.1), border=NA)
	abline(v=meanO.early,col=lineColor)
	text(meanO.early, textPos*max(yOnset.early), paste(round(meanO.early,1)), pos=2)
	abline(v=meanC.early,col=lineColor)
	text(meanC.early, textPos*max(yOnset.early), paste(round(meanC.early,1)), pos=4)

	#Duration early
	plot(xv, yDuration.early, type="l", xlim=d.x.lim, ylab="Density", xlab="Duration (days, early years)")
	polygon(xv, yDuration.early, col=durationFill, border=NA)
	abline(v=meanD.early,col=lineColor)
	text(meanD.early, textPos*max(yDuration.early), paste(round(meanD.early,1)), pos=4)

	#Onset and Cessation late
	plot(xv, yOnset.late, type="l", ylab="Density", xlab="Onset DOY (recent years)", xlim=c(minX,maxX),col=rgb(1,0,0,1))
	polygon(xv, yOnset.late, col=rgb(1,0,0,0.1), border=NA)
	points(fCs.late$x,yCessation.late, type='l',col=rgb(0,0,1,1))
	polygon(fCs.late$x, yCessation.late, col=rgb(0,0,1,0.1), border=NA)
	abline(v=meanO.late,col=lineColor)
	text(meanO.late, textPos*max(yOnset.late), paste(round(meanO.late,1)), pos=2)
	abline(v=meanC.late,col=lineColor)
	text(meanC.late, textPos*max(yOnset.late), paste(round(meanC.late,1)), pos=4)

	#Duration late
	plot(xv, yDuration.late, type="l", xlim=d.x.lim, ylab="Density", xlab="Duration (days, recent years)")
	polygon(xv, yDuration.late, col=durationFill, border=NA)
	abline(v=meanD.late,col=lineColor)
	text(meanD.late, textPos*max(yDuration.late), paste(round(meanD.late,1)), pos=4)

	if(!is.na(outFile)) dev.off()
}


create_PDFGraph_MAP_beta_beta = function(fileData, fileMCMCFirst, fileMCMCSecond, m=0, M=365, yearCutoff=1950, n=2000, outFile=NA, minX=50, maxX=200, lineColor="yellow", textPos=0.9)
{
	d = read.table(fileData, header=T, sep='\t')
	d = removeOutliers(d)
	d = getDOY_LatitudeAdjusted(d)
	DOY = d$dayOfYear.LatitudeAdjusted + mean(d$dayOfYear)
	DOY1 = DOY[d$date<yearCutoff]
	DOY2 = DOY[d$date>=yearCutoff]

	d1 = read.table(fileMCMCFirst, header=T, sep='\t')
	d2 = read.table(fileMCMCSecond, header=T, sep='\t')

	ML.Onset = list(shape1=d1$onset.shape1[1], shape2=d1$onset.shape2[1])
	ML.Duration = list(shape1=d1$duration.shape1[1], shape2= d2$duration.shape2[1])

	dx = 0.1
	xv = seq(from=m,to=M,by=dx)
	len = length(xv)
	fX = get.f_X(m=m, M=M, n=n, onset="beta", duration="beta", oParList=ML.Onset, dParList=ML.Duration)
	dx1 = (M-m) / n


	MAP.early = d1[which.max(d1$LnPosterior),]
	MAP.late = d2[which.max(d2$LnPosterior),]
	print(MAP.early)
	print(MAP.early)

	MAP.Onset.early = list(shape1=MAP.early$onset.shape1[1], shape2=MAP.early$onset.shape2[1])
	MAP.Duration.early = list(shape1=MAP.early$duration.shape1[1], shape2= MAP.early$duration.shape2[1])
	print(MAP.Onset.early)
	print(MAP.Duration.early)
	#fX.early = get.f_X(m=m, M=M, n=n, onset="beta", duration="beta", oParList=MAP.Onset.early, dParList=MAP.Duration.early)
	#y = fX.early(xv)
	#EfX.late = (xv %*% (y*dx)) - xv[1]*y[1]*dx/2 - xv[len]*y[len]*dx/2

	oS1 = MAP.Onset.early$shape1
	oS2 = MAP.Onset.early$shape2
	dS1 = MAP.Duration.early$shape1
	dS2 = MAP.Duration.early$shape2
	oParList = list(shape1 = oS1, shape2 = oS2)
	dParList = list(shape1 = dS1, shape2 = dS2)

	fs = get.f_X.f_C(m=m, M=M, n=n, onset="beta", duration="beta", oParList=oParList, dParList=dParList)
	fX.early = fs$f_X
	fCs.early = fs$fCs
	yX.early = fX.early(xv)
	yOnset.early = dbeta((xv-m)/(M-m), shape1=oS1, shape2=oS2)/(M-m)
	meanO.early = betaMean(shape1=oS1, shape2=oS2)*(M-m)+m
	yCessation.early = fCs.early$y
	meanC.early = (fCs.early$x %*% (fCs.early$y*dx1)) - fCs.early$x[1]*fCs.early$y[1]*dx1/2 - fCs.early$x[length(fCs.early$x)]*fCs.early$y[length(fCs.early$x)]*dx1/2

	#f_Dy = f_D(xv, m=m, M=M, onset="beta", duration="beta", oParList=oParList, dParList=dParList, n=n)
	#yDuration.early = yDuration.early + f_Dy
	#meanD = calcDMean(m=m, M=M, oS1=oS1, oS2=oS2, dS1=dS1, dS2=dS2,n=n)
	#meanD.early = meanD.early + meanD

	MAP.Onset.late = list(shape1=MAP.late$onset.shape1[1], shape2=MAP.late$onset.shape2[1])
	MAP.Duration.late = list(shape1=MAP.late$duration.shape1[1], shape2= MAP.late$duration.shape2[1])
	#fX.late = get.f_X(m=m, M=M, n=n, onset="beta", duration="beta", oParList=MAP.Onset.late, dParList=MAP.Duration.late)
	#y = fX.late(xv)
	#EfX.late = (xv %*% (y*dx)) - xv[1]*y[1]*dx/2 - xv[len]*y[len]*dx/2

	oS1 = MAP.Onset.late$shape1
	oS2 = MAP.Onset.late$shape2
	dS1 = MAP.Duration.late$shape1
	dS2 = MAP.Duration.late$shape2
	oParList = list(shape1 = oS1, shape2 = oS2)
	dParList = list(shape1 = dS1, shape2 = dS2)

	fs = get.f_X.f_C(m=m, M=M, n=n, onset="beta", duration="beta", oParList=oParList, dParList=dParList)
	fX.late = fs$f_X
	fCs.late = fs$fCs
	yX.late = fX.late(xv)
	yOnset.late = dbeta((xv-m)/(M-m), shape1=oS1, shape2=oS2)/(M-m)
	meanO.late = betaMean(shape1=oS1, shape2=oS2)*(M-m)+m
	yCessation.late = fCs.late$y
	meanC.late = (fCs.late$x %*% (fCs.late$y*dx1)) - fCs.late$x[1]*fCs.late$y[1]*dx1/2 - fCs.late$x[length(fCs.late$x)]*fCs.late$y[length(fCs.late$x)]*dx1/2

	#f_Dy = f_D(xv, m=m, M=M, onset="beta", duration="beta", oParList=oParList, dParList=dParList, n=n)
	#yDuration.late = yDuration.late + f_Dy
	#meanD = calcDMean(m=m, M=M, oS1=oS1, oS2=oS2, dS1=dS1, dS2=dS2,n=n)
	#meanD.late = meanD.late + meanD


	if(!is.na(outFile)) pdf(outFile)
	layout(matrix(c(1,1,1,1,2,2,3,3,4,4,5,5), nrow=3, ncol=4, byrow=T))

	#ML all data
	y=fX(xv)
	hd = hist(DOY, xlim=c(minX,maxX), probability=T, plot=F, breaks=seq(minX,maxX,by=5),col=rgb(1,0,1,0.1)) 
	hist(DOY, xlim=c(minX,maxX), ylim=c(0,max(y,hd$density)), probability=T, main=NULL, xlab="Observed DOY", breaks=seq(minX,maxX,by=5), col=rgb(1,0,1,0.1)) 
	points(xv, y, type='l',col=rgb(1,0,1,1))
	abline(v=mean(DOY), col=lineColor)
	text(mean(DOY), textPos*max(fX(xv)), paste(round(mean(DOY),1)), pos=4)

	#MCMC average early
	hd = hist(DOY1, xlim=c(minX,maxX), probability=T, plot=F, breaks=seq(minX,maxX,by=5)) 
	hist(DOY1, xlim=c(minX,maxX), ylim=c(0, max(yX.early,hd$density)), probability=T, main=NULL, xlab=paste("Observed DOY ( <",yearCutoff, ")"), ylab="Density", breaks=seq(minX,maxX,by=5),col=rgb(1,0,1,0.1)) 
	points(xv, yX.early, type='l',col=rgb(1,0,1,1))
	abline(v=mean(DOY1), col=lineColor)
	text(mean(DOY1), textPos*max(yX.early), paste(round(mean(DOY1),1)), pos=4)

	#MCMC average late
	hd = hist(DOY2, xlim=c(minX,maxX), probability=T, breaks=seq(minX,maxX,by=5), plot=F) 
	hist(DOY2, xlim=c(minX,maxX), ylim = c(0, max(yX.late, hd$density)), probability=T, main=NULL, xlab=paste("Observed DOY ( >=", yearCutoff, ")"), ylab="Density", breaks=seq(minX,maxX,by=5),col=rgb(1,0,1,0.1)) 
	points(xv, yX.late, type='l',col=rgb(1,0,1,1))
	abline(v=mean(DOY2), col=lineColor)
	text(mean(DOY2), textPos*max(yX.late), paste(round(mean(DOY2),1)), pos=4)

	#Onset and Cessation early
	plot(xv, yOnset.early, type="l", ylab="Density", xlab="Onset DOY (early years)", xlim=c(minX,maxX),col=rgb(1,0,0,1))
	polygon(xv, yOnset.early, col=rgb(1,0,0,0.1), border=NA)
	points(fCs.early$x,yCessation.early, type='l',col=rgb(0,0,1,1))
	polygon(fCs.early$x, yCessation.early, col=rgb(0,0,1,0.1), border=NA)
	abline(v=meanO.early,col=lineColor)
	text(meanO.early, textPos*max(yOnset.early), paste(round(meanO.early,1)), pos=2)
	abline(v=meanC.early,col=lineColor)
	text(meanC.early, textPos*max(yOnset.early), paste(round(meanC.early,1)), pos=4)

	#Duration early
	#plot(xv, yDuration.early, type="l", xlim=d.x.lim, ylab="Density", xlab="Duration (days, early years)")
	#polygon(xv, yDuration.early, col=durationFill, border=NA)
	#abline(v=meanD.early,col=lineColor)
	#text(meanD.early, textPos*max(yDuration.early), paste(round(meanD.early,1)), pos=4)

	#Onset and Cessation late
	plot(xv, yOnset.late, type="l", ylab="Density", xlab="Onset DOY (recent years)", xlim=c(minX,maxX),col=rgb(1,0,0,1))
	polygon(xv, yOnset.late, col=rgb(1,0,0,0.1), border=NA)
	points(fCs.late$x,yCessation.late, type='l',col=rgb(0,0,1,1))
	polygon(fCs.late$x, yCessation.late, col=rgb(0,0,1,0.1), border=NA)
	abline(v=meanO.late,col=lineColor)
	text(meanO.late, textPos*max(yOnset.late), paste(round(meanO.late,1)), pos=2)
	abline(v=meanC.late,col=lineColor)
	text(meanC.late, textPos*max(yOnset.late), paste(round(meanC.late,1)), pos=4)

	#Duration late
	#plot(xv, yDuration.late, type="l", xlim=d.x.lim, ylab="Density", xlab="Duration (days, recent years)")
	#polygon(xv, yDuration.late, col=durationFill, border=NA)
	#abline(v=meanD.late,col=lineColor)
	#text(meanD.late, textPos*max(yDuration.late), paste(round(meanD.late,1)), pos=4)

	#layout(matrix(c(1,1,1,1,2,2,3,3,5,5,6,6), nrow=3, ncol=4, byrow=T))

	##ML all data
	#hist(DOY, xlim=c(m,M), ylim=c(0,0.05), probability=T, main=NULL, xlab="Observed DOY") 
	#points(xv, fX(xv), type='l')
	#abline(v=mean(DOY), col="red")
	#text(mean(DOY), 0.045, paste(round(mean(DOY),1)), pos=4)
	#
	##MAP early
	#hist(DOY1, xlim=c(m,M), ylim=c(0,0.05), probability=T, main=NULL, xlab=paste("Observed DOY ( <",yearCutoff, ")"), ylab="Density") 
	#points(xv, fX.early(xv), type='l')
	#abline(v=mean(DOY1), col="red")
	#text(mean(DOY1), 0.045, paste(round(mean(DOY1),1)), pos=4)
	#
	##MAP late
	#hist(DOY2, xlim=c(m,M), ylim=c(0,0.05), probability=T, main=NULL, xlab=paste("Observed DOY ( >=", yearCutoff, ")"), ylab="Density") 
	#points(xv, fX.late(xv), type='l')
	#abline(v=mean(DOY2), col="red")
	#text(mean(DOY2), 0.045, paste(round(mean(DOY2),1)), pos=4)
	#
	##Onset early
	#y=dbeta((xv-m)/(M-m), MAP.Onset.early$shape1, MAP.Onset.early$shape2)/(M-m)
	#plot(xv, y, type="l", ylab="Density", xlab="Onset DOY (early years)")
	#mean = betaMean(MAP.Onset.early$shape1, MAP.Onset.early$shape2)*(M-m)+m
	#abline(v=mean,col="red")
	#text(mean, 0.9*max(y), paste(round(mean,1)), pos=4)
	#
	##duration calculations
	##mean.d.late = betaMean(MAP.Duration.late$shape1, MAP.Duration.late$shape2)*(mean-m)+m
	##sd.d.late = sqrt(betaVariance(MAP.Duration.late$shape1, MAP.Duration.late$shape2))*(mean-m)
	##mean.d.early = betaMean(MAP.Duration.early$shape1, MAP.Duration.early$shape2)*(mean-m)+m
	##sd.d.early = sqrt(betaVariance(MAP.Duration.early$shape1, MAP.Duration.early$shape2))*(mean-m)
	##
	##d.x.lim = c(0,max(mean.d.late,mean.d.early)+5*max(sd.d.late,sd.d.early))
	#
	##Duration early
	##!!these calculations are wrong!!
	##y = dbeta((xv-m)/(M-m), MAP.Duration.early$shape1, MAP.Duration.early$shape2)/(M-m)
	##plot(xv, y, type="l", xlim=d.x.lim, ylab="Density", xlab="Duration (days, early years)")
	##abline(v=mean.d.early,col="red")
	##text(mean.d.early, 0.9*max(y), paste(round(mean.d.early,1)), pos=4)
	#
	##Onset late
	#y=dbeta((xv-m)/(M-m), MAP.Onset.late$shape1, MAP.Onset.late$shape2)/(M-m)
	#plot(xv, y, type="l", ylab="Density", xlab="Onset DOY (recent years)")
	#mean = betaMean(MAP.Onset.late$shape1, MAP.Onset.late$shape2)*(M-m)+m
	#abline(v=mean,col="red")
	#text(mean, 0.9*max(y), paste(round(mean,1)), pos=4)
	#
	##Duration late
	##y = dbeta((xv-m)/(M-m), MAP.Duration.late$shape1, MAP.Duration.late$shape2)/(M-m)
	##plot(xv, y, type="l", xlim=d.x.lim, ylab="Density", xlab="Duration (days, recent years)")
	##abline(v=mean.d.late,col="red")
	##text(mean.d.late, 0.9*max(y), paste(round(mean.d.late,1)), pos=4)
	#
	if(!is.na(outFile)) dev.off()
}


create_PDFGraph_PosteriorAverage_beta_beta(fileData="../Data/PrimulaMeadia.txt", fileMCMCFirst="../MCMCReplicates/PrimulaMeadia.txt.FirstHalf.SO.txt.Cleaned.ForTracer.txt", fileMCMCSecond="../MCMCReplicates/PrimulaMeadia.txt.SecondHalf.SO.txt.Cleaned.ForTracer.txt", outFile="PrimulaMeadia.PosteriorAverage.PDFDecomposition.pdf")

