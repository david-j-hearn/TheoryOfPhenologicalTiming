source("../SharedCode/theoryOfPhenology.txt.R")
source("../SharedCode/inferPhenology_BayesianMCMC.txt.R")
source("../SharedCode/removeOutliers.txt.R")
source("../SharedCode/biasCorrectionScripts.txt.R")

thin = function(MCMCData, burnin=1500, thinInterval=25)
{
	MCMCData = MCMCData[burnin:nrow(MCMCData), ]
	MCMCData = MCMCData[seq(1, nrow(MCMCData), thinInterval), ]
	return(MCMCData)
}

#marginalized D
calcDMean = function(m=0, M=365, oS1, oS2, dS1, dS2,n=2000)
{
dx = (M-m)/n
t = seq(m, M, by = dx)
tot = 0
len = length(t)
for(i in 1:len)
	{
	for(j in 1:len)
		{
		if(t[i] + t[j] <= M && t[i]<M)
			{
			mult = 1
			if(i==1 || i==len) {
			       if(j==1 || j==len) { 
					mult = 0.25 
				} else {
					mult=0.5
			       }
			} else if(j==1 || j==len) { 
				mult = 0.5 
				}
			pd = dbeta(t[j]/(M-t[i]), shape1=dS1, shape2=dS2)/(M-t[i])
			po = dbeta((t[i]-m)/(M-m), shape1=oS1, shape2=oS2)/(M-m)
			tot = tot + mult * (t[j] * pd * po) * dx * dx
			}
		}
	}
return(tot)
}

#marginalized D
calcDVar = function(m=0, M=365, oS1, oS2, dS1, dS2,n=2000, mean)
{
dx = (M-m)/n
t = seq(m, M, by = dx)
tot = 0
len = length(t)
for(i in 1:len)
	{
	for(j in 1:len)
		{
		if(t[i] + t[j] <= M && t[i]<M)
			{
			mult = 1
			if(i==1 || i==len) {
			       if(j==1 || j==len) { 
					mult = 0.25 
				} else {
					mult=0.5
			       }
			} else if(j==1 || j==len) { 
				mult = 0.5 
				}
			pd = dbeta(t[j]/(M-t[i]), shape1=dS1, shape2=dS2)/(M-t[i])
			po = dbeta((t[i]-m)/(M-m), shape1=oS1, shape2=oS2)/(M-m)
			tot = tot + mult * (t[j] * t[j] * pd * po) * dx * dx
			}
		}
	}
return( tot - mean*mean )
}

#defunct
get_CI_lim = function(alpha = 0.05, upper=T, data)
{
	#print( data )
	alpha = alpha/2
	if(upper) alpha = 1-alpha
	data = sort(data)
	len = length(data)
	ind =  len*alpha 
	if(ind<1) ind = ind+1
	if(ind>=len) ind = len-1
	indl = floor(ind)
	if(ind == indl) return(data[ind])
	indh = indl+1
	lim = (data[indl] + data[indh])/2
	return( lim )
}

