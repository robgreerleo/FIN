#############################################################
## The first 3 function are requiared for Homework 8
#############################################################
##				coef.plot()
##__________________________________________________________
##
## 	This function is to plot R.Squared and betas from 
##	fitting FF 3 factor model 
##	the usage is coef.plot(R.Sq,betas)
##	R.Sq: length N vector; 
##	beta: 3 by N matrix; coef(fit)[-1]; fit is lm object
##	Change cols only if you want to change colors
##__________________________________________________________
##
coef.plot = function(R.Sq, betas, cols =c("#3690C0", "#EF3B2C" , "#4DAF4A"),  labs = syb, factors = c("Mkt","SMB", "HML")){
	if(dim(betas)[1] != 3) stop("Dimension of betas should be 3 by N", call. = FALSE)
    layout(matrix(c(1,2,3,4), 1, 4, byrow = TRUE), widths=c(1.25,1, 1, 1));
    par(mar = c(5, 4, 4, 0) + 0.1)
    barplot(R.Sq, space = 0.3, las = 2, horiz = T, col = cols[2], border = cols[2])
    mtext("R Squared")
    par(mar = c(5, 1, 4, 0.4) + 0.1)
    barplot(betas[1,], space =  0.3, las = 2, horiz = T, col = cols[1], border = cols[1], axisnames = F); 
    abline(v = 1, col = "grey");  
    mtext(paste("beta-",factors[1]))
    barplot(betas[2,], space = 0.3, las = 2, horiz = T, col = cols[1], border = cols[1], axisnames = F);
    mtext(paste("beta-",factors[2]))
    barplot(betas[3,], space =  0.3, las = 2, horiz = T, col = cols[1], border = cols[1], axisnames = F)
    mtext(paste("beta-",factors[3]))
    par(mar = c(5, 4, 4, 2) + 0.1)
}

#############################################################
##
##				resid.summary()
##__________________________________________________________
##
## 	This function plot the heatmap and test  
##	correlations of residuals
##  usage: resid.summary(res); res is an n by N matrix 
##__________________________________________________________
##

resid.summary<-function(res, labs = syb){
	dimnames(res)[[2]] = labs
	N = dim(res)[2]
	sig = cor.test.all(res)
	sum1 = paste("\nSignificant pairs at 1% level: ", sig[[1]], "of ", choose(N,2), "pairs" )
	sum2 = paste("\nSignificant pairs at 5% level: ", sig[[2]], "of ", choose(N,2), "pairs" )
	sum = cat(sum1, sum2)
	sum.plot = cor.plot(cor(res))
	invisible(list(sum = sum, plot = sum.plot ))
}

#############################################################
##
##				RSq.plot()
##__________________________________________________________
##
## 	This function plot R Squared values from 3 models
##  usage: RSq.plot(RSq.all); 
##  RSq.all is a 3 by N matrix
##	The rows are ordered as models
##__________________________________________________________
##

RSq.plot = function(RSq.all,  models = c("FF", "FA", "PC"), labs = syb, cols =c("#3690C0", "#EF3B2C" , "#4DAF4A") ){
	RSq = RSq.all
	if(dim(RSq)[1] != length(models)) RSq = t(RSq)
	dimnames(RSq)[[2]] = labs; dimnames(RSq)[[1]] = models;
	par(las = 1, pty = "m", mar = c(3.6,4.1,0.6,2.1), cex.axis = 0.8, mgp = c(2,0.75,0))
	barplot(RSq, horiz = T, beside = T, border = cols, col = cols, xlab = "R Squared")
	legend("topright", legend = rev(models), fill = rev(cols), bty = "n")
	par(las = 0, mar = c(5.1, 4.1, 4.1, 2.1), cex.axis = 1, mgp = c(3,1,0))
}
#############################################################
## 	functions for resid.summary
#############################################################
## 
##	You will not need to use these 2 functions
##

## test correlation for all pairs
cor.test.all = function(res){
	N = dim(res)[2]
	sig.01 = 0; sig.05 = 0;
	for(i in 1:(N-1)){
		for(j in (i+1):N){
			test.pval = cor.test(res[,i], res[,j])$p.value
			sig.01 = sig.01 + (test.pval < 0.01)
			sig.05 = sig.05 + (test.pval < 0.05)
		}
	}
	list(sig.01 = sig.01, sig.05 = sig.05)
}


#############################################################
## 	functions for resid.summary
#############################################################
## heatmap for correlatio matrix

cor.plot = function(cor){
  cor0 = cor*upper.tri(cor)
  c.range = round(range(cor0, na.rm = T),2)
	c.max = max(abs(c.range));
	c.ind = (c.max+c.range)*100+1
	cols = cm.colors(n = c.max*200+2, alpha = 1)
	heatmap(cor0, col = cols[c.ind[1]:c.ind[2]],  Rowv = NA, Colv = NA, revC = T,  symm = T, margin = c(5,5))
	legend.txt = as.character(c(c.range[1], 0, c.range[2]))
	legend.col = cols[c(c.ind[1], c.max*100+1, c.ind[2])]
	legend(mean(par('usr')[1:2]), par('usr')[4], bty='n', 
	       legend=legend.txt, fill=legend.col, border=T,   cex=.8)
	}


