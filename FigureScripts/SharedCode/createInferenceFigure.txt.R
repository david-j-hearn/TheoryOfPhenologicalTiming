source("../SharedCode/simulatePhenology_edit.txt.R")
source("../SharedCode/theoryOfPhenology.txt.R")
source("../SharedCode/inferPhenology_BayesianMCMC.txt.R")

#must be beta onset, beta duration
createInferenceFigure = function(data, mcmcSample, burnin, trueValues, trueOParList, trueDParList, min=0, max=365, N=100000, n=2000, SS=20, plotPDF=FALSE)
{
	onset="beta"
	duration="beta"

	mcmcItrs = length(mcmcSample$oParams)

	sampleable = mcmcItrs-burnin

	indices = sample(x=1:sampleable, size=SS, replace=FALSE)+burnin

	m=min
	M=max

	xl = min(data) - 50
	xh = max(data) + 50
	yl = 0
	yh = 0.1

	hist(data,probability=T,col=rgb(1,0,1,0.1), breaks=seq(xl,xh,5), ylim=c(yl,yh), main=NULL, xlab="Day of Year", ylab="Frequency") 
	legend("topleft", c(expression('O'['k=1']), "O", "X", "C", expression('C'['k=N']), "Theory", "True"), 
	       bty = "n", 
	       angle = c(0, 0, 0, 0, 0, 0, 0), density = c(100, 100, 100, 100, 100, 100, 100),
	       fill=c(rgb(1,204/255,0,1),rgb(1,0,0,1),rgb(1,0,1,1),rgb(0,0,1,1), rgb(0,1,1,1), "yellow", "black"))

	#Create Bayesian model averages (models weighted by posterior probabilities)
	for(i in indices) 
	{
		oParList = mcmcSample$oParams[[i]]
		dParList = mcmcSample$dParams[[i]]

		#make curves of inferred distributions
		fCs = f_C.vector(m=m, M=M, n=n, onset=onset, duration=duration, oParList=oParList, dParList=dParList)
		f_C = f_C.function(x=fCs$x, y=fCs$y)
		CDF = F_C.vector(x=fCs$x, y=fCs$y)
		F_C = F_C.function(x=fCs$x, y=CDF)
		Pts = P_t.vector(m=m, M=M, n=n, onset=onset, oParList=oParList, F_C.func = F_C)
		fXs = f_X.vector(x=fCs$x, Pts)

		#f_x
		if(plotPDF) { points(fCs$x, fXs, type="l",col="purple") }
		abline(v=sum(fCs$x*fXs/sum(fXs)),col="purple")
		#print(paste("Estimate of X:", sum(fCs$x*fXs/sum(fXs)))) 


		#f_o_scaled
		xs = seq(0,1,0.001)
		ys = dbeta(xs, shape1=oParList[[1]], shape2=oParList[[2]]) / (M - m)
		xs = m + (M-m)*xs
		if(plotPDF) { points(xs, ys , type="l", col="red",lwd=0.1) } #col=rgb(1,0,0,1))
		abline(v=sum(xs*ys/sum(ys)),col="red")
		#polygon(xs, ys, col=rgb(1,0,0,0.5) , border=NA)
		#print(paste("Estimate of O:", sum(xs*ys/sum(ys)))) #f_c_scaled

		#f_c_scaled
		if(plotPDF) { points(fCs$x, fCs$y, type="l", col="blue",lwd=0.1)  }#col=rgb(0,0,1,1))
		abline(v=sum(fCs$x*fCs$y/sum(fCs$y)),col="blue")
		#polygon(fCs$x, fCs$y, col=rgb(0,0,1,0.5) , border=NA)
		#print(paste("Estimate of C:", sum(fCs$x*fCs$y/sum(fCs$y)))) 

		#fxk=1 
		k=1
		fXk1_N = N*dbeta((fCs$x-m)/(M-m),oParList[[1]],oParList[[2]])*choose(N-1,k-1)*(pbeta((fCs$x-m)/(M-m),oParList[[1]],oParList[[2]])^(k-1))*((1-pbeta((fCs$x-m)/(M-m),oParList[[1]],oParList[[2]]))^(N-k))
		abline(v=sum(fCs$x*fXk1_N/sum(fXk1_N)), col="orange")
		#print(paste("Estimate of Xk1:", sum(fCs$x*fXk1_N/sum(fXk1_N)))) 
		fXk1_N = yl + (yh-yl)*fXk1_N/max(fXk1_N)
		if(plotPDF) { points(fCs$x, fXk1_N, type="l", col="orange",lwd=0.1)  }#col=rgb(1,204/255,0,1))
		#polygon(fCs$x, fXk1_N, col=rgb(1,204/255,0,0.5) , border=NA)


		#fxk=N
		k=N
		fXkN_N = f_X_k(fCs$x, N=N, k=k, f_X.func=f_C, F_X.func=F_C)
		abline(v=sum(fCs$x*fXkN_N/sum(fXkN_N)), col="cyan")
		#print(paste("Estimate of XkN:", sum(fCs$x*fXkN_N/sum(fXkN_N)))) 
		fXkN_N = yl + (yh-yl)*fXkN_N/max(fXkN_N)
		if(plotPDF) { points(fCs$x, fXkN_N, type="l", col="cyan",lwd=0.1) } #col=rgb(0,1,1,1))
		#polygon(fCs$x, fXkN_N, col=rgb(0,1,1,0.5) , border=NA)
	}


	#create the plots of the true values

	#population
	abline(v=trueValues$k1,col="yellow",lwd=2)
	abline(v=trueValues$o,col="yellow",lwd=2)
	abline(v=trueValues$x,col="yellow",lwd=2)
	abline(v=trueValues$c,col="yellow",lwd=2)
	abline(v=trueValues$kn,col="yellow",lwd=2)

	#theory
	fCs = f_C.vector(m=m, M=M, n=n, onset=onset, duration=duration, oParList=trueOParList, dParList=trueDParList)
	f_C = f_C.function(x=fCs$x, y=fCs$y)
	CDF = F_C.vector(x=fCs$x, y=fCs$y)
	F_C = F_C.function(x=fCs$x, y=CDF)
	Pts = P_t.vector(m=m, M=M, n=n, onset=onset, oParList=trueOParList, F_C.func = F_C)
	fXs = f_X.vector(x=fCs$x, Pts)

	#f_x
	points(fCs$x, fXs, type="l",col="black")
	mx = sum(fCs$x*fXs/sum(fXs))
	abline(v=mx,col="black", lwd=1)
	polygon(fCs$x, fXs, col=rgb(1,0,1,0.1) , border=NA) 

	#f_o_scaled
	xs = seq(0,1,0.001)
	ys = dbeta(xs, shape1=trueOParList[[1]], shape2=trueOParList[[2]]) / (M - m)
	xs = m + (M-m)*xs
	points(xs, ys , type="l", col="black")
	mo = sum(xs*ys/sum(ys))
	abline(v=mo,col="black", lwd=1)
	polygon(xs, ys, col=rgb(1,0,0,0.1) , border=NA) 

	#f_c_scaled
	points(fCs$x, fCs$y, type="l", col="black")
	mc = sum(fCs$x*fCs$y/sum(fCs$y))
	abline(v=mc,col="black", lwd=1)
	polygon(fCs$x, fCs$y, col=rgb(0,0,1,0.1) , border=NA) 

	#fxk=1 
	k=1
	fXk1_N = N*dbeta((fCs$x-m)/(M-m),trueOParList[[1]],trueOParList[[2]])*choose(N-1,k-1)*(pbeta((fCs$x-m)/(M-m),trueOParList[[1]],trueOParList[[2]])^(k-1))*((1-pbeta((fCs$x-m)/(M-m),trueOParList[[1]],trueOParList[[2]]))^(N-k))
	mxk1= sum(fCs$x*fXk1_N/sum(fXk1_N))
	abline(v=mxk1, col="black", lwd=1)
	fXk1_N = yl + (yh-yl)*fXk1_N/max(fXk1_N)
	points(fCs$x, fXk1_N, type="l", col="black")
	polygon(fCs$x, fXk1_N, col=rgb(1,204/255,0.1,0.1) , border=NA) 

	#fxk=N
	k=N
	fXkN_N = f_X_k(fCs$x, N=N, k=k, f_X.func=f_C, F_X.func=F_C)
	mxkN= sum(fCs$x*fXkN_N/sum(fXkN_N))
	abline(v=mxkN, col="black", lwd=1)
	fXkN_N = yl + (yh-yl)*fXkN_N/max(fXkN_N)
	points(fCs$x, fXkN_N, type="l", col="black")
	polygon(fCs$x, fXkN_N, col=rgb(0.1,1,1,0.1) , border=NA) 

}

