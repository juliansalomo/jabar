library(rgdal)
library(raster)
library(openxlsx)
library(spdep)

data.jabar = read.xlsx("Jabar Data (gabung).xlsx")
jabar2 <- readOGR(dsn = "petaJabar2" , layer = "Jabar2")

w<-poly2nb(jabar2)
ww<-nb2listw(w)
moran(jabar2$p.miskin16, ww, n=length(ww$neighbours), S0=Szero(ww))

moran(data.jabar$p.miskin16, ww, n=length(ww$neighbours), S0=Szero(ww))
moran.test(data.jabar$p.miskin16, ww,randomisation=T, alternative="greater")
moran.plot(data.jabar$p.miskin16, ww, labels=data.jabar$KABKOT)


k=16
colfunc <- colorRampPalette(c("green", "yellow","red"))
color <- colfunc(k)
library(sp)
jabar2$miskin2<- data.jabar$p.miskin16
spplot(jabar2, "miskin2", col.regions=color)

#eksplorasi data
library(arules)
k=4
p.miskin.baru<-discretize(data.jabar$p.miskin16,method="interval",categories=k)
jabar2$miskin<-p.miskin.baru
palette(gray(seq(0.5,0.8,len = k)))
legend("topright", levels(jabar2$miskin), pch=15,col=1:k)

#regresi klasik
reg.klasik<-lm(p.miskin16~EYS2016,data.jabar)
err.regklasik<-residuals(reg.klasik)

#uji asumsi regresi klasik
library(nortest)
library(car)
library(DescTools)
library(lmtest)
ad.test(err.regklasik) #menggunakan package "nortest"
hist(err.regklasik)
qqnorm(err.regklasik,datax=T)
qqline(rnorm(length(err.regklasik),mean(err.regklasik),sd(err.regklasik)),datax=T, col="red")
durbinWatsonTest(err.regklasik)
RunsTest(err.regklasik)
bptest(reg.klasik)

moran(err.regklasik, ww, n=length(ww$neighbours), S0=Szero(ww))
moran.test(err.regklasik, ww,randomisation=T, alternative="greater")

#LMtest
LM<-lm.LMtests(reg.klasik, nb2listw(w, style="W"),test=c("LMerr", "LMlag","RLMerr","RLMlag","SARMA"))

#pemodelan SAR
sar<-lagsarlm(p.miskin16~EYS2016,data=data.jabar,nb2listw(w))
err.sar<-residuals(sar)
ad.test(err.sar) #uji kenormalan,hrs menggunakan package "nortest"
bptest.sarlm(sar)    #uji kehomogenan residual
RunsTest(err.sar)

#pemodelan SEM
sem<-errorsarlm(p.miskin16~EYS2016,data=data.jabar,nb2listw(w))
err.sem<-residuals(sem)
ad.test(err.sem) #uji kenormalan,hrs menggunakan package "nortest"
bptest.sarlm(sem)    #uji kehomogenan residual
RunsTest(err.sem)

#pemodelan SEM
gsm<-sacsarlm(p.miskin16~EYS2016,data=data.jabar,nb2listw(w))
err.gsm<-residuals(gsm)
ad.test(err.gsm) #uji kenormalan,hrs menggunakan package "nortest"
bptest.sarlm(gsm)    #uji kehomogenan residual
RunsTest(err.gsm)

## ---- echo=FALSE, include=FALSE-----------------------------------------------
library(raster)
library(rgeos)
library(spdep)
library(rspatial)
library(raster)


kota <- aggregate(jabar2, "KABKOT")
library(latticeExtra)
library(RColorBrewer)

grps <- 10           #mengelompokkan menjadi 10 interval
brks <- quantile(jabar2$miskin2, 0:(grps-1)/(grps-1), na.rm=TRUE)

p <- spplot(jabar2, "miskin2", at=brks, col.regions=rev(brewer.pal(grps, "RdBu")), col="transparent" )
p + layer(sp.polygons(kota))




## -----------------------------------------------------------------------------
hh$fBadP <- pmax(hh$nBadPlumbi, hh$nBadKitche) / hh$nhousingUn
hh$fWhite <- hh$White / hh$Population
hh$age <- 2000 - hh$yearBuilt

f1 <- houseValue ~ age +  nBedrooms 
m1 <- lm(f1, data=hh)
summary(m1)


## -----------------------------------------------------------------------------
y <- matrix(hh$houseValue)
X <- cbind(1, hh$age, hh$nBedrooms)


## -----------------------------------------------------------------------------
ols <- solve(t(X) %*% X) %*% t(X) %*% y
rownames(ols) <- c('intercept', 'age', 'nBedroom')
ols


## ---- spreg6------------------------------------------------------------------
hh$residuals <- residuals(m1)

brks <- quantile(hh$residuals, 0:(grps-1)/(grps-1), na.rm=TRUE)

spplot(hh, "residuals", at=brks, col.regions=rev(brewer.pal(grps, "RdBu")), col="black")


## ---- spreg8------------------------------------------------------------------
library(spdep)
nb <- poly2nb(hh)
nb[[21]] <- sort(as.integer(c(nb[[21]], 38)))
nb[[38]] <- sort(as.integer(c(21, nb[[38]])))
nb

par(mai=c(0,0,0,0))
plot(hh)
plot(nb, coordinates(hh), col='red', lwd=2, add=TRUE)


## ---- spreg10-----------------------------------------------------------------
resnb <- sapply(nb, function(x) mean(hh$residuals[x]))
cor(hh$residuals, resnb)
plot(hh$residuals, resnb, xlab='Residuals', ylab='Mean adjacent residuals')
lw <- nb2listw(nb)


## -----------------------------------------------------------------------------
moran.mc(hh$residuals, lw, 999)


## ----spreg, message=FALSE-----------------------------------------------------
library(spatialreg )


## ----spregplot1---------------------------------------------------------------
m1s = lagsarlm(f1, data=hh, lw, tol.solve=1.0e-30)

summary(m1s)

hh$residuals <- residuals(m1s)
moran.mc(hh$residuals, lw, 999)

brks <- quantile(hh$residuals, 0:(grps-1)/(grps-1), na.rm=TRUE)
p <- spplot(hh, "residuals", at=brks, col.regions=rev(brewer.pal(grps, "RdBu")), col="transparent")
print( p + layer(sp.polygons(hh)) )


## ----spregplotx---------------------------------------------------------------
m1e <- errorsarlm(f1, data=hh, lw, tol.solve=1.0e-30)
summary(m1e)

hh$residuals <- residuals(m1e)
moran.mc(hh$residuals, lw, 999)

brks <- quantile(hh$residuals, 0:(grps-1)/(grps-1), na.rm=TRUE)
p <- spplot(hh, "residuals", at=brks, col.regions=rev(brewer.pal(grps, "RdBu")),
            col="transparent")
print( p + layer(sp.polygons(hh)) )



## ----spregplot3---------------------------------------------------------------
brks <- quantile(hh$residuals, 0:(grps-1)/(grps-1), na.rm=TRUE)

p <- spplot(hh, "residuals", at=brks, col.regions=rev(brewer.pal(grps, "RdBu")),
            col="transparent")

print( p + layer(sp.polygons(hh)) )


