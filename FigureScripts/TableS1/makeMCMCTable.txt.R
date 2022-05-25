#The below script generates Table S1 without the bells and whistles
#Be sure to have run the scripts to make Fig. 5 first!

#select the files from which to extract results of MCMC runs; results of MCMC runs are stored in ../MCMCReplicates
files = c("AnemoneQuinquefolia.txt", "DicentraCucullaria.txt", "PrimulaMeadia.txt",
	  "CamassiaScilloides.txt", "EnemionBiternatum.txt", "SanguinariaCanadensis.txt",
	  "CardamineConcatenata.txt","ErythroniumAmericanum.txt","ThalictrumThalictroides.txt",
	  "ClaytoniaVirginica.txt", "MertensiaVirginica.txt", "CollinsiaVerna.txt", "PodophyllumPeltatum.txt");

if(!require(HDInterval)) { install.packages("HDInterval") } 
library(HDInterval)

lines = "Species\tObserved.mean.early\tObserved.var.early\tOnset.mean.early\tHPD.Onset.mean.early\tOnset.var.early\tHPD.Onset.var.early\tDuration.mean.early\tHPD.Duration.mean.early\tDuration.var.early\tHPD.Duration.var.early\tObserved.mean.late\tObserved.var.late\tOnset.mean.late\tHPD.Onset.mean.late\tOnset.var.late\tHPD.Onset.var.late\tDuration.mean.late\tHPD.Duration.mean.late\tDuration.var.late\tHPD.Duration.var.late"

burnin = 5000
interval = 100
yearSplit=1950

cnt = 2
for(file in files)
{
	message(file)

	data1 = read.table(paste("../Data/",file,sep=""), header=T, sep='\t')
	data1 = removeOutliers(data1)
	data1 = getDOY_LatitudeAdjusted(data1)	#remove latitudinal biases

	DOY = data1$dayOfYear.LatitudeAdjusted
	DOY = DOY + mean(data1$dayOfYear) #recenter data on the mean value

	data = read.table(paste("../MCMCReplicates/", file, ".FirstHalf.SO.txt.Cleaned.ForTracer.txt",sep=""), header=T, sep='\t')
	message("pre-1950")

	len = nrow(data)
	inds = seq(burnin,len,interval)

	hpd.om.f = round(hdi(data$mean.onset[inds]),1)
	hpd.ov.f = round(hdi(data$variance.onset[inds]),1)
	hpd.dm.f = round(hdi(data$mean.duration[inds]),1)
	hpd.dv.f = round(hdi(data$variance.duration[inds]),1)
	om.f = round(mean(data$mean.onset[inds]),1)
	ov.f = round(mean(data$variance.onset[inds]),1)
	dm.f = round(mean(data$mean.duration[inds]),1)
	dv.f = round(mean(data$variance.duration[inds]),1)

	DOY = DOY[data1$date<yearSplit]
	meanX.f = round(mean(DOY),1)
	varX.f = round(var(DOY),1)

	data = read.table(paste("../MCMCReplicates/", file, ".SecondHalf.SO.txt.Cleaned.ForTracer.txt",sep=""), header=T, sep='\t')
	message("post-1950")

	len = nrow(data)
	inds = seq(burnin,len,interval)

	hpd.om.s = round(hdi(data$mean.onset[inds]),1)
	hpd.ov.s = round(hdi(data$variance.onset[inds]),1)
	hpd.dm.s = round(hdi(data$mean.duration[inds]),1)
	hpd.dv.s = round(hdi(data$variance.duration[inds]),1)
	om.s = round(mean(data$mean.onset[inds]),1)
	ov.s = round(mean(data$variance.onset[inds]),1)
	dm.s = round(mean(data$mean.duration[inds]),1)
	dv.s = round(mean(data$variance.duration[inds]),1)

	DOY = data1$dayOfYear.LatitudeAdjusted
	DOY = DOY + mean(data1$dayOfYear)
	DOY = DOY[data1$date>=yearSplit]
	meanX.s = round(mean(DOY),1)
	varX.s = round(var(DOY),1)

	lines[cnt] = paste(file,"\t",meanX.f,"\t",varX.f,"\t",om.f,"\t(",hpd.om.f[[1]],",",hpd.om.f[[2]],")\t", ov.f,"\t(",hpd.ov.f[[1]],",",hpd.ov.f[[2]],")\t",dm.f,"\t(",hpd.dm.f[[1]],",",hpd.dm.f[[2]],")\t",dv.f,"\t(",hpd.dv.f[[1]],",",hpd.dv.f[[2]],")\t",meanX.s,"\t",varX.s,"\t",om.s,"\t(",hpd.om.s[[1]],",",hpd.om.s[[2]],")\t", ov.s,"\t(",hpd.ov.s[[1]],",",hpd.ov.s[[2]],")\t",dm.s,"\t(",hpd.dm.s[[1]],",",hpd.dm.s[[2]],")\t",dv.s,"\t(",hpd.dv.s[[1]],",",hpd.dv.s[[2]],")",sep="")
	cnt = cnt+1
}

fileConn = file("MCMC_results.txt")
writeLines(lines, fileConn)
close(fileConn)

