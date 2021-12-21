##########################################################################
##
##	rate()
## 	Calculate coverage rate of one-step rolling forecasts form  
##		an uGARCHforecast class of object
##  Eg: fore = ugarchforecast(fit)
## 		use method: rate(fore)
##
##########################################################################
rate = function(object){
    UseMethod("rate", object)
}
setMethod("rate", 
    signature(object = "uGARCHforecast"),
    function(object){
        dist = object@model$modeldesc$distribution
        distpar = object@model$pars[c("lambda", "skew", "shape"),1];
        cri.u = qdist(dist, c(.95,.975), lambda = distpar["lambda"], skew = distpar["skew"], shape = distpar["shape"])
		cri.l = qdist(dist, c(.05,.025), lambda = distpar["lambda"], skew = distpar["skew"], shape = distpar["shape"])
        T0 = object@model$modeldata$T
        n.roll = object@forecast$n.roll
        ind.new = (T0+1):(T0+1+n.roll)
        y.new = object@model$modeldata$data[ind.new]
        y.hat = fitted(object)[1,]; 
        sig = sigma(object)[1,];
        out.90 = c(mean(y.new < (y.hat + cri.l[1]*sig)), mean(y.new > (y.hat + cri.u[1]*sig)))
        out.90 = c(1-sum(out.90), out.90)
        out.95 = c(mean(y.new < (y.hat + cri.l[2]*sig)), mean(y.new > (y.hat + cri.u[2]*sig)))
        out.95 = c(1-sum(out.95), out.95)
        out = rbind(out.95, out.90)
        dimnames(out)[[1]] = c("95% PI", "90% PI")
        dimnames(out)[[2]] = c("coverage", "below PI", "beyond PI")
        cat(paste("\n      One-Step Rolling Forecast   \n"))
        cat(paste("----------------------------------\n", sep = ""))
        print(round(out,4))
        invisible(out)
        }
)

##########################################################################
##
##	rate0()
## 	Calculate coverage rate of one-step rolling forecasts form  
##		an arfimaforecast class of object
##  Eg: fore = arfimaforecast(fit)
## 		use method: rate(fore)
##
##########################################################################
rate0 = function(object){
    UseMethod("rate0", object)
}
setMethod("rate0", 
    signature(object = "ARFIMAforecast"),
    function(object){
        dist = object@model$modeldesc$distribution
        distpar = object@model$pars[c("lambda", "skew", "shape"),1];
        cri.u = qdist(dist, c(.95,.975), lambda = distpar["lambda"], skew = distpar["skew"], shape = distpar["shape"])
		cri.l = qdist(dist, c(.05,.025), lambda = distpar["lambda"], skew = distpar["skew"], shape = distpar["shape"])
        T0 = object@model$modeldata$T
        n.roll = object@forecast$n.roll
        ind.new = (T0+1):(T0+1+n.roll)
        y.new = object@model$modeldata$data[ind.new]
        y.hat = fitted(object)[1,]; 
        sig = object@model$pars["sigma",1];
        out.90 = c(mean(y.new < (y.hat + cri.l[1]*sig)), mean(y.new > (y.hat + cri.u[1]*sig)))
        out.90 = c(1-sum(out.90), out.90)
        out.95 = c(mean(y.new < (y.hat + cri.l[2]*sig)), mean(y.new > (y.hat + cri.u[2]*sig)))
        out.95 = c(1-sum(out.95), out.95)
        out = rbind(out.95, out.90)
        dimnames(out)[[1]] = c("95% PI", "90% PI")
        dimnames(out)[[2]] = c("coverage", "below PI", "beyond PI")
        cat(paste("\n      One-Step Rolling Forecast   \n"))
        cat(paste("----------------------------------\n", sep = ""))
        print(round(out,4))
        invisible(out)
        }
)

#################################################################################
##
##   R package rugarch by Alexios Ghalanos Copyright (C) 2008-2015.
##   This file is based on R package rugarch by Alexios Ghalanos.
##
##   Define the function showShort for print shorter summary report
##   for ARFIMAfit object
##   
#################################################################################

showShort0 = function(object){
    UseMethod("showShort0", object)
}
setMethod("showShort0",
		signature(object = "ARFIMAfit"),
		function(object){
			model = object@model
			modelinc = model$modelinc
			cat(paste("\n*----------------------------------*", sep = ""))
			cat(paste("\n*          ARFIMA Model Fit        *", sep = ""))
			cat(paste("\n*----------------------------------*", sep = ""))
			cat("\nMean Model\t: ARFIMA(", modelinc[2],",",ifelse(modelinc[4]>0, "d", 0),",",modelinc[3],")\n", sep = "")
			cat("Distribution\t:", model$modeldesc$distribution,"\n")
			if(object@fit$convergence == 0){
				cat("\nOptimal Parameters")
				cat(paste("\n------------------------------------\n",sep=""))
				print(round(object@fit$matcoef,6), digits = 5)
				#cat("\nRobust Standard Errors:\n")
				#print(round(object@fit$robust.matcoef,6), digits = 5)
				#if(!is.null(object@fit$hessian.message)){
				#	cat(paste("\n", object@fit$hessian.message))
				#}
				cat("\nLogLikelihood :", object@fit$LLH, "\n")
				stdresid = object@fit$residuals/coef(object)["sigma"]
				itestm = t(infocriteria(object))
				cat("\nInformation Criteria")
				cat(paste("\n------------------------------------\n",sep=""))
				print(itestm,digits=5)
				cat("\nWeighted Ljung-Box Test on Standardized Residuals")
				cat(paste("\n------------------------------------\n",sep=""))
				tmp1 = .weightedBoxTest(stdresid, p = 1, df = sum(modelinc[2:3]))
				print(tmp1, digits = 4)
				cat("\nH0 : No serial correlation\n")
				cat("\nWeighted Ljung-Box Test on Standardized Squared Residuals")
				cat(paste("\n------------------------------------\n",sep=""))
				tmp2 = .weightedBoxTest(stdresid, p = 2, df = sum(modelinc[8:9]))
				print(tmp2, digits = 4)
				cat("\n\nARCH LM Tests")
				cat(paste("\n------------------------------------\n",sep=""))
				L2 = .archlmtest(stdresid, lags = 2)
				L5 = .archlmtest(stdresid, lags = 5)
				L10 = .archlmtest(stdresid, lags = 10)
				alm = matrix(0,ncol = 3,nrow = 3)
				alm[1,1:3] = c(L2$statistic, L2$parameter, L2$p.value)
				alm[2,1:3] = c(L5$statistic, L5$parameter, L5$p.value)
				alm[3,1:3] = c(L10$statistic, L10$parameter, L10$p.value)
				colnames(alm) = c("Statistic", "DoF", "P-Value")
				rownames(alm) = c("ARCH Lag[2]", "ARCH Lag[5]", "ARCH Lag[10]")
				print(alm,digits = 4)
			} else{
				cat("\nConvergence Problem:")
				cat("\nSolver Message:", object@fit$message,"\n\n")
				
			}
			return(invisible(object))
})

#################################################################################
##
##   R package rugarch by Alexios Ghalanos Copyright (C) 2008-2015.
##   This file is based on R package rugarch by Alexios Ghalanos.
##
##   Define the function showShort for print shorter summary report
##   for uGARCHfit object
##   
#################################################################################
#----------------------------------------------------------------------------------
# univariate show method / seperate for fit,sim and forecast
#----------------------------------------------------------------------------------
# fit show
library(WeightedPortTest)
showShort = function(object){
	UseMethod("showShort", object)
}
setMethod("showShort",
		signature(object = "uGARCHfit"),
		function(object){
			vmodel = object@model$modeldesc$vmodel
			model = object@model
			modelinc = object@model$modelinc
			cat(paste("\n*---------------------------------*", sep = ""))
			cat(paste("\n*          GARCH Model Fit        *", sep = ""))
			cat(paste("\n*---------------------------------*", sep = ""))
			#cat("\n\nConditional Variance Dynamics \t")
			#cat(paste("\n-----------------------------------", sep = ""))
			cat(paste("\nGARCH Model\t: ", vmodel,"(", modelinc[8], ",", modelinc[9], ")\n", sep=""))
			if(vmodel == "fGARCH"){
				cat(paste("fGARCH Sub-Model\t: ", model$modeldesc$vsubmodel, "\n", sep = ""))
			}
			cat("Mean Model\t: ARFIMA(", modelinc[2],",",ifelse(modelinc[4]>0, "d", 0),",",modelinc[3],")\n", sep = "")
			cat("Distribution\t:", model$modeldesc$distribution,"\n")
			if(object@fit$convergence == 0){
				cat("\nOptimal Parameters")
				cat(paste("\n------------------------------------\n",sep=""))
				print(round(object@fit$matcoef,6), digits = 5)
				if(!is.null(object@fit$hessian.message)){
					cat(paste("\n", object@fit$hessian.message))
				}
				cat("\nLogLikelihood :", object@fit$LLH, "\n")
				stdresid = object@fit$residuals/object@fit$sigma
				itestm = t(infocriteria(object))
				cat("\nInformation Criteria")
				cat(paste("\n------------------------------------\n",sep=""))
				print(itestm,digits=5)
				cat("\nWeighted Ljung-Box Test on Standardized Residuals")
				cat(paste("\n------------------------------------\n",sep=""))
				tmp1 = .weightedBoxTest(stdresid, p = 1, df = sum(modelinc[2:3]))
				print(tmp1, digits = 4)
				cat(paste("d.o.f=", sum(modelinc[2:3]), sep=""))
				cat("\nH0 : No serial correlation\n")
				cat("\nWeighted Ljung-Box Test on Standardized Squared Residuals")
				cat(paste("\n------------------------------------\n",sep=""))
				tmp2 = .weightedBoxTest(stdresid, p = 2, df = sum(modelinc[8:9]))
				print(tmp2, digits = 4)
				cat(paste("d.o.f=", sum(modelinc[8:9]), sep=""))
				cat("\n\nWeighted ARCH LM Tests")
				cat(paste("\n------------------------------------\n",sep=""))
				gdf = sum(modelinc[8:9])
				L2 = .weightedarchlmtest(residuals(object), sigma(object), lags  = gdf+1, fitdf=gdf)
				L5 = .weightedarchlmtest(residuals(object), sigma(object), lags  = gdf+3, fitdf=gdf)
				L10 = .weightedarchlmtest(residuals(object), sigma(object), lags = gdf+5, fitdf=gdf)
				alm = matrix(0, ncol = 4, nrow = 3)
				alm[1,1:4] = as.numeric(c(L2$statistic, L2$parameter, L2$p.value))
				alm[2,1:4] = as.numeric(c(L5$statistic, L5$parameter, L5$p.value))
				alm[3,1:4] = as.numeric(c(L10$statistic, L10$parameter, L10$p.value))
				colnames(alm) = c("Statistic", "Shape", "Scale", "P-Value")
				rownames(alm) = c(paste("ARCH Lag[",gdf+1,"]",sep=""), paste("ARCH Lag[",gdf+3,"]",sep=""), paste("ARCH Lag[",gdf+5,"]",sep=""))
				print(alm,digits = 4)

				}			

			invisible(object)
		})
        
#################################################################################
##
##   R package rugarch by Alexios Ghalanos Copyright (C) 2008-2015.
##   This file is part of the R package rugarch.
##
##   The R package rugarch is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package rugarch is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################

# Q-Statistics on Standardized Residuals
.box.test = function(stdresid, p=1, df = 0)
{
	if(any(!is.finite(stdresid))) stdresid[!is.finite(stdresid)]=0
	# p=1 normal case, p=2 squared std. residuals
	# Q-Statistics on Standardized Residuals
	#H0 : No serial correlation ==> Accept H0 when prob. is High [Q < Chisq(lag)]
	box10 = Box.test(stdresid^p, lag = 1, type = "Ljung-Box", fitdf = 0)
	box15 = Box.test(stdresid^p, lag = df+1, type = "Ljung-Box", fitdf = df)
	box20 = Box.test(stdresid^p, lag = df+5, type = "Ljung-Box", fitdf = df)
	LBSR<-matrix(NA,ncol=2,nrow=3)
	LBSR[1:3,1] = c(box10$statistic[[1]],box15$statistic[[1]],box20$statistic[[1]])
	LBSR[1:3,2] = c(box10$p.value[[1]],box15$p.value[[1]],box20$p.value[[1]])
	rownames(LBSR) = c(paste("Lag[1]",sep=""), paste("Lag[p+q+1][",df+1,"]",sep=""), paste("Lag[p+q+5][",df+5,"]",sep=""))
	colnames(LBSR) = c("statistic","p-value")
	return(LBSR)
}

.weightedBoxTest = function(stdresid, p=1, df = 0)
{
	if(any(!is.finite(stdresid))) stdresid[!is.finite(stdresid)]=0
	# p=1 normal case, p=2 squared std. residuals
	# Q-Statistics on Standardized Residuals
	#H0 : No serial correlation ==> Accept H0 when prob. is High [Q < Chisq(lag)]
	box10 = Weighted.Box.test(stdresid, lag = 1, type = "Ljung-Box", fitdf = 0, if(p==2) sqrd.res = TRUE else sqrd.res = FALSE)
	box15 = Weighted.Box.test(stdresid, lag = max(2, 2*df+df-1), type = "Ljung-Box", fitdf = df, if(p==2) sqrd.res = TRUE else sqrd.res = FALSE)
	box20 = Weighted.Box.test(stdresid, lag = max(5, 4*df+df-1), type = "Ljung-Box", fitdf = df, if(p==2) sqrd.res = TRUE else sqrd.res = FALSE)
	LBSR<-matrix(NA,ncol=2,nrow=3)
	LBSR[1:3,1] = c(box10$statistic[[1]],box15$statistic[[1]],box20$statistic[[1]])
	LBSR[1:3,2] = c(box10$p.value[[1]],box15$p.value[[1]],box20$p.value[[1]])
	rownames(LBSR) = c(paste("Lag[1]",sep=""), paste("Lag[2*(p+q)+(p+q)-1][",max(2, 2*df+df-1),"]",sep=""), paste("Lag[4*(p+q)+(p+q)-1][",max(5, 4*df+df-1),"]",sep=""))
	colnames(LBSR) = c("statistic","p-value")
	return(LBSR)
}


.archlmtest = function (x, lags, demean = FALSE)
{
	if(any(!is.finite(x))) x[!is.finite(x)] = 0
	x = as.vector(x)
	if(demean) x = scale(x, center = TRUE, scale = FALSE)
	lags = lags + 1
	mat = embed(x^2, lags)
	arch.lm = summary(lm(mat[, 1] ~ mat[, -1]))
	STATISTIC = arch.lm$r.squared * length(resid(arch.lm))
	names(STATISTIC) = "Chi-squared"
	PARAMETER = lags - 1
	names(PARAMETER) = "df"
	PVAL = 1 - pchisq(STATISTIC, df = PARAMETER)
	METHOD = "ARCH LM-test"
	result = list(statistic = STATISTIC, parameter = PARAMETER,
			p.value = PVAL, method = METHOD)
	class(result) = "htest"
	return(result)
}


.weightedarchlmtest = function (x, sigma, lags, fitdf = 2, demean = FALSE)
{
	if(any(!is.finite(x))) x[!is.finite(x)] = 0
	x = as.vector(x)
	if(demean) x = scale(x, center = TRUE, scale = FALSE)
	result = Weighted.LM.test(x, sigma^2, lag = lags, type = c("correlation", "partial")[1], fitdf = fitdf, weighted=TRUE) 
	return(result)
}

