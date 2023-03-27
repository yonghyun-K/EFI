tmpdata = sapply(as.data.frame(table(Y[,x], useNA = "always")), as.numeric)
s <- prelim.cat(tmpdata[,1:(ncol(tmpdata) - 1)], tmpdata[,ncol(tmpdata)]) # preliminary manipulations
thetahat <- em.cat(s, showits=F)


formula(paste("~", paste(names(Y)[x], collapse = "*")))

library(cvam)
x = c(1,2)
get.fitted(cvam(formula(paste("~", paste(names(Y)[x], collapse = "*"))), data = Y))$fit
get.fitted(cvam(formula(paste("~", paste(names(Y)[x], collapse = "*"))), data = Y))
cvamEstimate(~V1, cvam(formula(paste("~", paste(names(Y)[x], collapse = "*"))), data = Y))
cvamEstimate(~V3, cvam(formula(paste("~", paste(names(Y)[x], collapse = "*"))), data = Y))

get.fitted(cvam(~ V1 * V2 + V3, data = Y))$fit

array(get.fitted(cvam(formula(paste("~", paste(names(Y)[x], collapse = "*"))), data = Y))$fit, dim = supplen[x])



install.packages("imputeMulti")
library(imputeMulti)



data(tract2221)




summary(tract2221)

install.packages("HardyWeinberg")
library(HardyWeinberg)
data(Markers)


install.packages("cvam")
library(cvam)
data(crime)
get.fitted(cvam(~ V1 * V2, freq = n, data = crime))

get.fitted(cvam(~ V1 * V2, freq = n, data = crime))

get.fitted(cvam(~ V1 * V2, freq = n, data = crime))

formula("~ V1 * V2")

get.fitted(cvam(~ V1 + V2, freq = n, data = crime))


loglin(HairEyeColor, list(c(1, 2), c(1, 3), c(2, 3)))

fit <- cvam(~ Hair + Eye, freq = Freq, data = as.data.frame.table(HairEyeColor))
fit$freq
get.fitted(cvam(~ Hair * Eye, freq = Freq, data = as.data.frame.table(HairEyeColor), saturated = T))

get.fitted(cvam(~ Hair * Eye * Sex, freq = Freq, data = as.data.frame.table(HairEyeColor), saturated = T))

get.fitted(cvam(~ .[[1]] + .[[2]], freq = Freq, data = as.data.frame.table(HairEyeColor)))


as.data.frame.table(HairEyeColor)

as.data.frame.table(table(Y, useNA = "ifany"))

get.fitted(cvam(~ V1 + V2, freq = NULL, data = crime))

for(i in 1:1000){
  get.fitted(cvam(~ V1 + V2 + V1 * V2, freq = n, data = crime))
}


get.coef(cvam(~ V1 + V2 , freq = n, data = crime))

get.coef(cvam(~ V1 + V2 , freq = n, data = crime))


s <- prelim.cat(tmpdata[,1:(ncol(tmpdata) - 1)], tmpdata[,ncol(tmpdata)]) # preliminary manipulations
thetahat <- em.cat(s, showits=F)

View(cvam.cvam)

crimes
s <- prelim.cat(crimes[,1:2],crimes[,3]) # preliminary manipulations
thetahat <- em.cat(s)
thetahat
