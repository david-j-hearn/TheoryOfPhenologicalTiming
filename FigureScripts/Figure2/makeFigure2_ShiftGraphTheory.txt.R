source("../SharedCode/helperFunctions.txt.R")

if(!require(shape)) { install.packages("shape") }
library(shape)

m = 0
M = 365
n = 1000
xv = seq(0,1,(1-0)/n)
xv.s = seq(m,M,(M-m)/n)

arr.length=0.25
arr.adj=1

getYs = function(oM,oV,dM,dV,m=0,M=365,n=1000)
{
xv = seq(0,1,1/n)
xv.s = seq(m,M,(M-m)/n)
oPars = getShapeParameters(oM, oV, m=m, M=M)
dPars = getShapeParameters(dM, dV, m=0, M=M-oM)


fOv = dbeta(xv, shape1=oPars$shape1, shape2=oPars$shape2)/(M-m)
fDv = f_D(xv.s, m=m, M=M, onset="beta", duration="beta", oParList=oPars, dParList=dPars, n=n)
dM.s = calcDMean(m=m, M=M, oS1=oPars$shape1, oS2=oPars$shape2, dS1=dPars$shape1, dS2=dPars$shape2,n=n)
fX = get.f_X(m=m, M=M, n=n, onset="beta", duration="beta", oParList=oPars, dParList=dPars)
fXv = fX(xv.s)
return(list(xv.s=xv.s, fOv=fOv, fDv=fDv, dM.s=dM.s, fXv=fXv))
}

createPlots = function(text, e, l, oM.e, oM.l, arr.length=0.25, arr.adj=1, minO=60, maxO=150, minD=0, maxD=90, minX=60, maxX=180, text.col="black")
{
message(text)
xv = e$xv.s
len = length(xv)
dx = xv[2]-xv[1]


#onset
fOv.e = e$fOv
fOv.l = l$fOv
maxY = max(fOv.e, fOv.l)
xlim = c(minO,maxO)
plot(xv, fOv.e, col="red", xaxt='n',yaxt='n',ylab='',xlab='', type='l', xlim=xlim, ylim=c(0,maxY))
polygon(xv, fOv.e, col=rgb(1,0,0,0.1), border=NA)
points(xv, fOv.l, col="red", xaxt='n',yaxt='n',ylab='',xlab='', type='l')
polygon(xv, fOv.l, col=rgb(1,0,0,0.1), border=NA)
abline(v=oM.e, col="red")
abline(v=oM.l, col="red")
y = max(fOv.e,fOv.l)/2
Arrows(y0=y, x0=oM.e, y1=y, x1=oM.l, col="black", arr.length=arr.length, arr.adj=arr.adj)
text(minO,0.8*maxY, text, pos=4, col=text.col)

#duration
fDv.e = e$fDv
fDv.l = l$fDv
dM.s.e = e$dM.s
dM.s.l = l$dM.s
maxY = max(fDv.e, fDv.l)
xlim = c(minD,maxD)
plot(xv, fDv.e, col="black", xaxt='n',yaxt='n',ylab='',xlab='', type='l', xlim=xlim, ylim=c(0,maxY))
polygon(xv, fDv.e, col=rgb(0,0,0,0.1), border=NA)
points(xv, fDv.l, col="black", xaxt='n',yaxt='n',ylab='',xlab='', type='l')
polygon(xv, fDv.l, col=rgb(0,0,0,0.1), border=NA)
abline(v=dM.s.e, col="black")
abline(v=dM.s.l, col="black")
y = max(fDv.e,fDv.l)/2
Arrows(y0=y, x0=dM.s.e, y1=y, x1=dM.s.l, col="black", arr.length=arr.length, arr.adj=arr.adj)

#observed
fXv.e = e$fXv
fXv.l = l$fXv
maxY = max(fXv.e, fXv.l)
xlim = c(minX,maxX)
xM.s.e = xv %*% (fXv.e*dx) - xv[1]*fXv.e[1]*dx/2 - xv[len]*fXv.e[len]*dx/2
xM.s.l = xv %*% (fXv.l*dx) - xv[1]*fXv.l[1]*dx/2 - xv[len]*fXv.l[len]*dx/2
plot(xv, fXv.e, col="purple", xaxt='n',yaxt='n',ylab='',xlab='', type='l', xlim=xlim, ylim=c(0,maxY))
polygon(xv, fXv.e, col=rgb(1,0,1,0.1), border=NA)
points(xv, fXv.l, col="purple", xaxt='n',yaxt='n',ylab='',xlab='', type='l')
polygon(xv, fXv.l, col=rgb(1,0,1,0.1), border=NA)
abline(v=xM.s.e, col="purple")
abline(v=xM.s.l, col="purple")
y = max(fXv.e,fXv.l)/2
Arrows(y0=y, x0=xM.s.e, y1=y, x1=xM.s.l, col="black", arr.length=arr.length, arr.adj=arr.adj)
}

#______________________________________
#_______CREATE THE PLOT________________

pdf("Figure2_ShiftTheoryGraph.pdf")

layout( matrix(c( 1,2,3, 4,5,6, 7,8,9, 10,11,12, 13,14,15, 16,17,18), nrow=6, byrow=T))

par(mar = c(0,0,0,0))

#case 1 + + + 
oM1.e=90
oM1.l=120
c1.e = getYs(oM=oM1.e,oV=100,dM=20,dV=100,m=m,M=M,n=n)
c1.l = getYs(oM=oM1.l,oV=100,dM=50,dV=100,m=m,M=M,n=n)
createPlots(text="+ + +", e=c1.e, l=c1.l, oM.e=oM1.e, oM.l=oM1.l, arr.length=arr.length, arr.adj=arr.adj)

#case 2 + - +
oM2.e=90
oM2.l=120
c2.e = getYs(oM=oM2.e,oV=100,dM=50,dV=100,m=m,M=M,n=n)
c2.l = getYs(oM=oM2.l,oV=100,dM=20,dV=100,m=m,M=M,n=n)
createPlots(text="+ - +", e=c2.e, l=c2.l, oM.e=oM2.e, oM.l=oM2.l, arr.length=arr.length, arr.adj=arr.adj)

#case 3 + - - 
oM3.e=90
oM3.l=100
c3.e = getYs(oM=oM3.e,oV=100,dM=80,dV=100,m=m,M=M,n=n)
c3.l = getYs(oM=oM3.l,oV=100,dM=20,dV=100,m=m,M=M,n=n)
createPlots(text="+ - -\nParadox", e=c3.e, l=c3.l, oM.e=oM3.e, oM.l=oM3.l, arr.length=arr.length, arr.adj=arr.adj, text.col="purple")

#case 4 - + +
oM4.e=100
oM4.l=90
c4.e = getYs(oM=oM4.e,oV=100,dM=20,dV=100,m=m,M=M,n=n)
c4.l = getYs(oM=oM4.l,oV=100,dM=80,dV=100,m=m,M=M,n=n)
createPlots(text="- + +\nParadox", e=c4.e, l=c4.l, oM.e=oM4.e, oM.l=oM4.l, arr.length=arr.length, arr.adj=arr.adj, text.col="purple")

#case 5 - + -
oM5.e=120
oM5.l=90
c5.e = getYs(oM=oM5.e,oV=100,dM=20,dV=100,m=m,M=M,n=n)
c5.l = getYs(oM=oM5.l,oV=100,dM=50,dV=100,m=m,M=M,n=n)
createPlots(text="- + -", e=c5.e, l=c5.l, oM.e=oM5.e, oM.l=oM5.l, arr.length=arr.length, arr.adj=arr.adj)

#case 6 - - - 
oM6.e=120
oM6.l=90
c6.e = getYs(oM=oM6.e,oV=100,dM=50,dV=100,m=m,M=M,n=n)
c6.l = getYs(oM=oM6.l,oV=100,dM=20,dV=100,m=m,M=M,n=n)
createPlots(text="- - -", e=c6.e, l=c6.l, oM.e=oM6.e, oM.l=oM6.l, arr.length=arr.length, arr.adj=arr.adj)

dev.off()

#_____________
#case 7 0 0 0 
#case 8 0 + +
#case 9 0 - - 
#case10 - 0 - 
#case11 + 0 +
#_____________
#case12 + - 0
#case13 - + 0


