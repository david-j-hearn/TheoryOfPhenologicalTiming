source("../SharedCode/simulatePhenology_edit.txt.R")
source("../SharedCode/theoryOfPhenology.txt.R")
source("../SharedCode/inferPhenology_BayesianMCMC.txt.R")
source("../SharedCode/createInferenceFigure.txt.R")

if(!require(fitdistrplus)) { install.packages("fitdistrplus") }
library(fitdistrplus)

#min, max range of days
m=0
M=365

#simulated population size and sample size from simulated population, assuming beta distributed onset and durations, given parameters
N=1000000
subSS = 500 
onset="beta"
duration="beta"

#parameters similar to Mertensia
oParListTrue = list(shape1=88.0676336684406, shape2= 235.433702711231 ) 
dParListTrue = list(shape1=77.7122061756863, shape2= 518.474949834894 )


message(paste("RMVD: True params for simulation "))
print(oParListTrue)
print(dParListTrue)

print(paste("OMean:", betaMean(oParListTrue[[1]], oParListTrue[[2]])))
print(paste("Ovar:", betaVariance(oParListTrue[[1]], oParListTrue[[2]])))
print(paste("dMean:", betaMean(dParListTrue[[1]], dParListTrue[[2]])))
print(paste("dvar:", betaVariance(dParListTrue[[1]], dParListTrue[[2]])))

message("RMVD: Simulating phenology data")
simData = simulatePhenology(N=N, subSS = subSS, min=m, max=M, onset = onset, duration=duration, oParList=oParListTrue, dParList=dParListTrue)
data = simData$observations
trueParams = list(k1=min(simData$onsets),o=mean(simData$onsets),x=mean(data),c=mean(simData$onsets+simData$durations),kn=max(simData$onsets+simData$durations),ovar=var(simData$onsets),cvar=var(simData$onsets+simData$durations),xvar=var(data))

#MCMC iterations
nIter = 10000

#graphing parameters
graphInterval = 100

message("RMVD: Running MCMC to infer parameters of simulated data")
mcmc.output.est = mcmc(
		       data = data,
		       onset = onset,
		       duration = duration,
		       min = m,
		       max = M,
		       logDPrior = logDPrior.uniform, 		#prior function
		       priorParList = list(min=1, max=2000),	#prior hyperparameters
		       useMLToInit=TRUE,
		       proposalRate = 0.025,
		       proposalAdapt = 0.0,
		       useBetaProposal = TRUE,	#parameterizes in terms of beta variance and mean rather than shape parameters; tends to get better estimates
		       nIter = nIter,
		       graphResults = FALSE, #set to true to visualize progress of MCMC if you are running from an R instance with graphics 
		       graphInterval = graphInterval,
		       graphPopulationSize = N,
		       trueOParList = oParListTrue,
		       trueDParList = dParListTrue,
		       trueParams = trueParams
)

printMCMCDataToFile(file="MCMC.output.Fig3.txt", mcmcOutput=mcmc.output.est, nIter=nIter, thin=10, type="betabeta",min=m, max=M, trueOParams=oParListTrue, trueDParams=dParListTrue)


message("Creating inference graph")

#increments for numerical routines
n=2000	#for numerical precision
burnin=nIter/2
SS=10

pdf("Figure3_BayesianInference.pdf")
createInferenceFigure(
		      data=data, 
		      mcmcSample=mcmc.output.est, 
		      burnin=burnin, 
		      trueValues = trueParams, 
		      trueOParList = oParListTrue, 
		      trueDParList = dParListTrue, 
		      min=m, 
		      max=M, 
		      N=N, 
		      n=n, 
		      SS=SS)
dev.off()



