####################
##### SET-UP ######
###################

#Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Downloading packages
list.of.packages <- c("readr", "zoo", "tseries", "stargazer", "fUnitRoots", "dplyr",
                      "aTSA", "xtable", "forecast", "ellipse", "graphics", "knitr", 
                      "rmarkdown", "markdown")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

invisible(lapply(list.of.packages, library, character.only = TRUE))

#Reset environment
rm(list=ls())


###############################################################
############## IMPORTING & TRANSFORMING THE DATA ##############
###############################################################

#Importing the data
serie <- read.csv("valeurs_mensuelles.csv", sep=";")

#Creating the serie of dates
dates <- as.yearmon(seq(from=1990+0/12,to=2023+1/12,by=1/12)) 

#Inverting the time values
serie <- serie[dim(serie)[1]:1,]
serie$Time <- seq(1,nrow(serie),1)


#Transforming our serie into a zoo object
prod <- zoo(serie$Values, order.by=dates)

#We drop the last two observations for predictions
prod_train <- prod[1:(length(prod)-2)]
plot(prod_train, xlab = "Time", ylab = "Values")

#We create the differenciated serie
dprod <- diff(prod_train,1)
plot(cbind(prod_train,dprod))


###################################################
############## PART 1 : STATIONARITY ##############
###################################################

#We regress the values of our serie on time to check if there is a trend
reg <- lm(Values ~ Time, data=serie)
stargazer(reg, type="text")

#Time isn't significative, there is no trend, so we do the adf test with only a constant: "c"

adf <- adfTest(prod_train, lag=0, type="c") 
adf 

#The unit root hypothesis is not rejected
#However, this test isn't valid because we didn't consider residuals autocorrelation


#We define two functions to determine the number of lags we have to consider to get rid of autocorrelation

#Qtest
Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}

#adfTest
adfTest_valid <- function(series,kmax,type){ #ADF tests until no more autocorrelated residuals
  k <- 0
  noautocorr <- 0
  while (noautocorr==0){
    cat(paste0("ADF with ",k, " lags: residuals OK? "))
    adf <- adfTest(series,lags=k,type=type)
    pvals <- Qtests(adf@test$lm$residuals,kmax,fitdf=length(adf@test$lm$coefficients))[,2]
    if (sum(pvals<0.05,na.rm=T) == 0) {
      noautocorr <- 1; cat("OK \n")}
    else cat("No \n")
    k <- k + 1
  }
  return(adf)
}

Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients)) 

#All residuals of the adf test with 0 lags are autocorrelated
#We need to fix a number of lags to consider for the adf test
#To do so we use the function adfTest_valid defined above

adfTest_valid(prod_train, 24, "c") 

#The test indicates that we have to consider 5 lags to obtain a valid adfTest
#The adf test with 5 lags reject the non stationarity assumption but with a p-value close to 0.05
#We will consider other tests to confirm the (or not) the stationarity of the series

pp.test(x=as.vector(prod_train), output=TRUE) #Phillips-Perron test
kpss.test(x=as.vector(prod_train)) #KPSS

#We will thus differenciate our serie to obtain a better rejected adf test

#Creating lags and first difference
serie$lag1 <- lag(serie$Values)
serie$dif <- serie$Values - serie$lag1

dprod <- diff(prod_train,1)

#Again, we run a regression to test for the presence of a time trend
reg_diff <- lm(dif ~ Time, data=serie)
stargazer(reg_diff, type="text")

#There is no time trend and significant constant so we run an adf test of type "nc"
dadf <- adfTest(dprod, lag=0, type="nc") 
dadf

#We reject the non stationarity but we didn't consider lags so the test isn't valid
#Again, all the residuals are correlated if we use an adf test with 0 lags

Qtests(dadf@test$lm$residuals, 24, fitdf = length(dadf@test$lm$coefficients)) 

#We run the same test as before to determine the number of lags we have to consider

adfTest_valid(dprod,24,"nc")

#Thus, we consider 4 lags and the test is strongly rejected (pvalue < 0.01)

pp.test(x=as.vector(dprod), output=TRUE) #Phillips-Perron test
kpss.test(x=as.vector(dprod)) #KPSS

#Phillips-Perron and KPSS tests confirm the stationarity of the differenciated series


###################################################
############## PART 2 : ARMA MODELS ###############
###################################################

#We plot the acf and pacf to determine qmax and pmax
par(mfrow=c(1,2))
acf(dprod);pacf(dprod)
dev.off()

#From the plot we see that
qmax <- 2
pmax <- 5

#We create a loop to calculate the AIC/BIC for each possibility in the grid (qmax,pmax)=(2,5)
mat <- matrix(NA,nrow=pmax+1,ncol=qmax+1) #empty matrix to fill
rownames(mat) <- paste0("p=",0:pmax) #renames lines
colnames(mat) <- paste0("q=",0:qmax) #renames columns
AICs <- mat #AIC matrix not filled non remplie
BICs <- mat #BIC matrix not filled non remplie
pqs <- expand.grid(0:pmax,0:qmax) #all possible combinations of p and q
for (row in 1:dim(pqs)[1]){ #loop for each (p,q)
  p <- pqs[row,1] #gets p
  q <- pqs[row,2] #gets q
  estim <- try(arima(dprod,c(p,0,q),include.mean = F)) #tries ARIMA estimation
  AICs[p+1,q+1] <- if (class(estim)=="try-error") NA else estim$aic #assigns the AIC
  BICs[p+1,q+1] <- if (class(estim)=="try-error") NA else BIC(estim) #assigns the BIC
}

#Print the AICs
AICs
AICs==min(AICs)

xtable(AICs) #We get the latex table

#The ARIMA(5,1,0) minimizes the AIC. We thus keep it.
arima510 <- arima(prod_train,c(5,1,0),include.mean=F)

#Print the AICs
BICs
BICs==min(BICs)

xtable(BICs) #We get the latex table

#The ARIMA(1,1,1) minimizes the BIC. We then keep it.
arima111 <- arima(prod_train,c(1,1,1),include.mean=F)

#Tables of the ARIMA
stargazer(arima510, type="text")
stargazer(arima111, type="text")

#However, we will carry out a further study to determine the best model. 
#We determine the significance of the coefficients of the ARs and MAs, and analyse the autocorrelation of the residuals.

#To do so, we define three functions:

p_value <- #This function returns the p-values of the estimated ARMA
  function(estim) {
    coef <- estim$coef 
    se <- sqrt(diag(estim$var.coef)) 
    t <- coef / se 
    pval <- (1 - pnorm(abs(t)))*2 
    return(rbind(coef, se, pval))
  }

model_choice <- #This function estimates an arima and checks the fit and validity of the model with p-value
  function(p, q, data = dxm, k = 24) {
    estim <-
      try(arima(prod_train, c(p, 1, q), optim.control = list(maxit = 20000)))
    if (class(estim) == "try-error")
      return(c(
        "p" = p,
        "q" = q,
        "arsignif" = NA,
        "masignif" = NA,
        "resnocorr" = NA,
        "ok" = NA
      ))
    arsignif <- if (p == 0) 
      NA
    else
      p_value(estim)[3, p] <= 0.05 
    masignif <- if (q == 0) 
      NA
    else
      p_value(estim)[3, p + q] <= 0.05 
    resnocorr <-
      sum(Qtests(estim$residuals, 30, length(estim$coef) - 1)[, 2] <= 0.05, na.rm =
            T) == 0 
    checks <- c(arsignif, masignif, resnocorr) 
    ok <-
      as.numeric(sum(checks, na.rm = T) == (3 - sum(is.na(checks)))) 
    return(
      c(
        "p" = p,
        "q" = q,
        "arsignif" = arsignif,
        "masignif" = masignif,
        "resnocorr" = resnocorr,
        "ok" = ok
      )
    )
  }


arma_model_choice <- #This function runs the previous one with all p<pmax & q<qmax
  function(pmax, qmax) {
    pqs <- expand.grid(0:pmax, 0:qmax) 
    t(apply(matrix(1:dim(pqs)[1]), 1, function(row) { 
      p <- pqs[row, 1]
      q <- pqs[row, 2]
      cat(paste0("Computing ARMA(", p, ",", q, ") \n"))
      model_choice(p, q) 
    }))
  }

arma_models <- arma_model_choice(pmax, qmax) 
arma_models

#We only keep models with the column "ok" equals to 1 (significant coefficients & non-autocorrelated residuals) 

selection <- arma_models[arma_models[, "ok"] == 1 & !is.na(arma_models[, "ok"]),] 
selection 

#Thus, we will study between the four following models: ARMA(5,0), ARMA(1,1), ARMA(4,1), and ARMA(0,2)

pqs <- apply(selection, 1, function(row) list("p" = as.numeric(row[1]), "q" = as.numeric(row[2]))) 
names(pqs) <- paste0("arma(", selection[, 1], ",", selection[, 2], ")") 
models <- lapply(pqs, function(pq) arima(dprod, c(pq[["p"]], 0, pq[["q"]]))) 
vapply(models, FUN.VALUE = numeric(2), function(m) c("AIC" = AIC(m), "BIC" = BIC(m)))

#Again we choose according to the AIC / BIC 
#We know that the ARMA(5,0) minimizes the AIC
#And the ARMA(1,1) minimizes the BIC
#Thus, we keep these two models

#Test of validity of the ARMA:

qtest_arima510 <- Qtests(arima510$residuals,24,fitdf= length(arima510$coef)-1)
xtable(qtest_arima510)
qtest_arima111 <- Qtests(arima111$residuals,24,fitdf= length(arima111$coef)-1)
xtable(qtest_arima111)

#We never reject the absence of residual autocorrelation. The ARIMA(5,1,0) is valid.

#Choice between the two models:
#We create a function that calculate the Adjusted R2
adj_r2 <- function(model){
  ss_res <- sum(model$residuals^2) #sum of the squared residuals
  p <- model$arma[1] #gets the AR order
  q <- model$arma[2] #gets the MA order
  ss_tot <- sum(dprod[-c(1:max(p,q))]^2) #sum of the observations from the sample squared.
  n <- model$nobs-max(p,q) #sample size
  adj_r2 <- 1-(ss_res/(n-p-q-1))/(ss_tot/(n-1)) #r2 ajust´e
  return(adj_r2)
}

adj_r2(arima510)
adj_r2(arima111)
#The ARIMA(5,1,0) has the largest adjusted R2, it gives the best prediction in the sample. 
#We keep it as the best model.

#We take a look at the residuals of our arima, thanks to the function checkresiduals 
checkresiduals(arima510)

#We plot the inverse of our roots
model <- Arima(prod,c(5,1,0))
plot(model)
#All the roots of our model are > 1


#######################################
###### PLOTTING OUR ESTIMATION ########
#######################################

serie$lag1 <- lag(serie$Values)
serie$lag2 <- lag(serie$lag1)
serie$lag3 <- lag(serie$lag2)
serie$lag4 <- lag(serie$lag3)
serie$lag5 <- lag(serie$lag4)
serie$lag6 <- lag(serie$lag5)

phi1 <- arima510$coef[1]
phi2 <- arima510$coef[2]
phi3 <- arima510$coef[3]
phi4 <- arima510$coef[4]
phi5 <- arima510$coef[5]

serie$pred <- (1+phi1)*serie$lag1 + (phi2-phi1)*serie$lag2 +
  (phi3-phi2)*serie$lag3 + (phi4-phi3)*serie$lag4 +
  (phi5-phi4)*serie$lag5 - phi5*serie$lag6

predict <- zoo(serie$pred,order.by=dates)
predict <- predict[1:(length(predict)-2)]

plot(prod_train, col = "black", xlab = "Time", ylab = "Values")
lines(predict, col = "red")
legend("topleft", legend = c(expression(X[t]), "ARIMA(5,1,0)"), col = c("black", "red"), lty = 1)


#We take a look at the QQplot of our residuals
residuals <- prod - predict
residuals <- na.omit(residuals)

qqnorm(residuals)
qqline(residuals)


##########################################
########## PART 3 : PREDICTION ###########
##########################################

#We supposed that our residuals are gaussian
#We want to estimate their standard error
sigma <- sd(residuals)
sigma

#The standard error of our residuals is sigma=11

#We seek to predict the two next values of our serie thanks to the ARIMA model
model_pred <- predict(arima510, n.ahead=2) #prediction
pred <- zoo(model_pred$pred, order.by=as.yearmon(c(2023+0/12,2023+1/12))) #transforming our predictions into a zoo object
link <- rbind(prod_train[length(prod_train)],pred[1]) 

#We create a function that plots the prediction, starting the serie at a given date
plot_pred <- function(start){
  plot(prod, col = "black", ylab = "Values", xlab="Time",
       xlim = c(start, 2023+3/12))
  U <- model_pred$pred + 1.96*model_pred$se
  L <- model_pred$pred - 1.96*model_pred$se
  Upper1 <- model_pred$pred[1] + 1.96*sigma
  Lower1 <- model_pred$pred[1] - 1.96*sigma
  Upper2 <- model_pred$pred[2] + 1.96*sigma*sqrt((1+(1+phi1)^2))
  Lower2 <- model_pred$pred[2] - 1.96*sigma*sqrt((1+(1+phi1)^2))
  Upper <- c(Upper1,Upper2)
  Upper <- zoo(Upper, order.by=as.yearmon(c(2023+0/12,2023+1/12)))
  Lower <- c(Lower1,Lower2)
  xx <- c(time(Upper), rev(time(Upper)))
  yy <- c(Lower,rev(Upper))
  polygon(xx,yy,border="8", col=gray(0.6,alpha=0.2))
  lines(pred, type = "p", col = "red")
  lines(pred, type = "l", col = "red")
  lines(link, type = "l", col = "red")
  legend("topleft", legend=c(expression(X[t]), "Prediction"), col=c("black", "red"), lty=1)
}

#We plot our prediction, starting the serie in January 2020
plot_pred(2020)


######################################
############## ELLIPSE ###############
######################################

Sigma <- matrix(data=c(sigma^2, sigma^2*(1+phi1), sigma^2*(1+phi1), sigma^2*(1+(1+phi1)^2)), nrow=2, byrow=TRUE)

legend_labels <- c("Border of the confidence region at level 95%", "Predicted values", "True values")
legend_colors <- c("black", "black", "red")


ell <- ellipse(x=Sigma, centre = c(model_pred$pred[1], model_pred$pred[2]), t = sqrt(qchisq(0.95, 2)))
plot(ell, type='l', xlab = expression(X[T+1]), ylab = expression(X[T+2]))
points(x=model_pred$pred[1], y=model_pred$pred[2], pch=4)
points(x=as.vector(prod)[397], y=as.vector(prod)[398], pch=4, col="red")
legend("topleft", legend=legend_labels, col=legend_colors, lty=c(1,0,0), pch=c(NA,4,4))

