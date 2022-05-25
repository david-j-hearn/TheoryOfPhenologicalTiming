source("../SharedCode/theoryOfPhenology.txt.R")
source("../SharedCode/simulatePhenology_edit.txt.R")
source("../SharedCode/removeOutliers.txt.R")
source("../SharedCode/inferPhenology_BayesianMCMC.txt.R")

m=0 		#minimum time in days
M=365 		#maximum time in days
n=2000		#for numerical precision
N=100000	#simulated population size
N_Pt=1000	#population size for graph of P_N_t
k=1		
t = 110		#day of year used for graph of P_N_t
subSS = 10000	#sample size from simulated data
onset="beta"
duration="beta"

#used in figure - similar to Mertensia with no corrections for bias - used simply as an illustration of densities
oParList=list(shape1=90.27611,shape2=244.8296)	
dParList=list(shape1=42.84,shape2=233.6238)	

#get the densities based on the true parameters
fCs = f_C.vector(m=m, M=M, n=n, onset=onset, duration=duration, oParList=oParList, dParList=dParList)
f_C = f_C.function(x=fCs$x, y=fCs$y)
CDF = F_C.vector(x=fCs$x, y=fCs$y)
F_C = F_C.function(x=fCs$x, y=CDF)
Pts = P_t.vector(m=m, M=M, n=n, onset=onset, oParList=oParList, F_C.func = F_C)
P_t =   P_t.function(x=fCs$x, y=Pts)
fXs = f_X.vector(x=fCs$x, Pts)
f_X = f_X.function(x=fCs$x, y=fXs)
FXs = F_X.vector( x=fCs$x, y=fXs) 
F_X = F_X.function(x=fCs$x, y=FXs)

Nts = f_N_t(t=t, N=N_Pt, P_t.func = P_t)

#simulate the data based on true parameters
observations = simulatePhenology(N=N, subSS = subSS, min=m, max=M, onset = onset, duration=duration, oParList=oParList, dParList=dParList)

#read in actual data (note: there were no corrections for bias -> results in slightly different parameters than rigorous analyses that generated Table S1)
data = read.table("../Data/MertensiaVirginica.txt", header=T, sep='\t')
data = removeOutliers(data)
data = data$dayOfYear

#set colors for figure
cols = c(rgb(1,204/255,0,1), rgb(1,0,0,1), rgb(1,0,1,1), rgb(1,0,1,1), rgb(0,0,1,1), rgb(0,1,1,1), rgb(1,0,1,1), rgb(0,1,0,1), rgb(155/255,103/255,60/255,1))
cols.t = c(rgb(1,204/255,0,0.1), rgb(1,0,0,0.1), rgb(1,0,1,0.1), rgb(1,0,1,0.1), rgb(0,0,1,0.1), rgb(0,1,1,0.1), rgb(1,0,1,0.1), rgb(0,1,0,0.1),  rgb(155/255,103/255,60/255,0.1))

pdf("Figure4_MertensiaVirginica_TOPTDensities.pdf")

layout(matrix(c(1,1,1,2,2, 1,1,1,3,3, 1,1,1,4,4),3,5,byrow=T))

par(mar=c(5,5,0,.5))

message("hist of raw data, probability is true")
pm=50
pM=200
hist(data, probability=T, xlim=c(50,200), ylim=c(0,0.1), breaks=seq(m,M,5), border=cols[4], col=cols.t[4], main=NULL, xlab = "Day of year")

message("f_x")
points(fCs$x, fXs, type="l",col=cols[3])
polygon(fCs$x, fXs, col=cols.t[3], border=NA)

message("f_O_scaled")
xs = seq(0,1,0.001)
points(m + (M-m)* xs, dbeta(xs, shape1=oParList[[1]], shape2=oParList[[2]]) / (M - m), type="l", col=cols[2])
polygon(m + (M-m)* xs, dbeta(xs, shape1=oParList[[1]], shape2=oParList[[2]]) / (M - m), col=cols.t[2] , border=NA)

message("f_C_scaled")
points(fCs$x, fCs$y, type="l", col=cols[5])
polygon(fCs$x, fCs$y, col=cols.t[5] , border=NA)

message("fOk=1")
yl=0
yh = 0.075
        k=1
        fXk1_N = N*dbeta((fCs$x-m)/(M-m),oParList[[1]],oParList[[2]])*choose(N-1,k-1)*(pbeta((fCs$x-m)/(M-m),oParList[[1]],oParList[[2]])^(k-1))*((1-pbeta((fCs$x-m)/(M-m),oParList[[1]],oParList[[2]]))^(N-k))
        fXk1_N = yl + (yh-yl)*fXk1_N/max(fXk1_N)
        points(fCs$x, fXk1_N, type="l", col=cols[1])
        polygon(fCs$x, fXk1_N, col=cols.t[1] , border=NA)

message("fCk=N")
        k=N
        fXkN_N = f_X_k(fCs$x, N=N, k=k, f_X.func=f_C, F_X.func=F_C)
        fXkN_N = yl + (yh-yl)*fXkN_N/max(fXkN_N)
        points(fCs$x, fXkN_N, type="l", col=cols[6])
        polygon(fCs$x, fXkN_N, col=cols.t[6] , border=NA)

segments(250,0,365,0,col="white",lwd=3)
segments(0,0,250,0,col="black",lwd=1)

legend(pm,0.1, legend=c(expression('PDF O'['k=1']*': N=100,000'), "PDF O", "PDF X", "Histogram observed data", "PDF C", expression('PDF C'['k=N']*': N=100,000'), "Histogram simulated data", expression('PMF P'['N'['t']]*': N=1,000 t=110'), "PDF R, N=100,000"), col=cols, lwd=1.5)

message("inset plot 1: simulated data")
xm=50
xM=200
plot(fCs$x,fXs, type='l', col=cols[7], xlab="Day of year of observation", ylab="Frequencies of simulated values", xlim=c(xm,xM))
hist(observations$observations, probability=T, add=T,border=rgb(0,0,0,0.5), col=cols.t[7], breaks=seq(xm,xM, 5))

message("inset plot 2: P_N_t")
x=850:950
plot(x,Nts[x],type="l", ylab=paste("Probability of individuals\nin phenophase at day", t), xlab=paste("Individuals in phenophase"), col=cols[8])
polygon(x, Nts[x], col = cols.t[8], border=NA)

#make the plot of the duration of phenophase at the population level
message("Making ranges PDF")
fXk1_N = fXk1_N/sum(fXk1_N)
fXkN_N = fXkN_N/sum(fXkN_N)

O1.func =  function_Factory(x=fCs$x, y=fXk1_N)
CN.func =  function_Factory(x=fCs$x, y=fXkN_N)

message("getting pdf of ranges. This may take a while")
y.ranges = populationPhenologyDuration.vector(func.O1 = O1.func, func.CN = CN.func, min=m, max=M, N=N, n=1000)

#simulate 1000 populations of size n
message(paste("Simulating 1000 populations of size. This may take a while" , N,collate=" "))
reps=1000
ranges = simulateRanges(reps = reps, N=N, min=m, max=M, onset = onset, duration=duration, oParList=oParList, dParList=dParList)
message("Histogram of the ranges")
h=hist(ranges$ranges, plot=FALSE)
mh = max(h$density)
#scaled to fit histogram data
y.ranges$density = (mh+0.0075*mh)*y.ranges$density/max(y.ranges$density)
plot(y.ranges$tRs, y.ranges$density, type="l", col=cols[9], xlab="Temporal range of phenophase in population", ylab="Frequencies of simulated values", xlim=c(100,150))
hist(ranges$ranges, probability=T, add=T,border=rgb(0,0,0,0.5), col=cols.t[9])

dev.off()

