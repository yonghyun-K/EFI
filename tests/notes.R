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



# install.packages("imputeMulti")
library(imputeMulti)



data(tract2221)
p = ncol(tract2221)
cand.edges = as.list(data.frame(combn(p, 2)))
dp = doublep(tract2221, cand.edges)
EFI = efi(tract2221, dp)

# Y = cbind(tract2221, 1)

Y = as.data.frame.table(table(tract2221, useNA = "ifany"))

plot(dp)

nrow(Y[!(Y$Freq == 0),])
nrow(tract2221)

prod(sapply(tract2221, nlevels))

summary(tract2221)


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

##########03272023

lapply(Y[-ncol(Y)], function(x) xtabs(~ x, Y) / n)

get.fitted(cvam(formula(paste("~", paste(namesY[x], collapse = "*"))), data = Y, freq = Freq),
           mfTrue = F)

get.modelMatrix(cvam(formula(paste("~", paste(namesY[x], collapse = "*"))), data = Y, freq = Freq),
                msgIfNone  = F)

Y

get.mfTrue(cvam(formula(paste("~", paste(namesY[x], collapse = "*"))), data = Y, freq = Freq))

cvamImpute(cvam(formula(paste("~", paste(namesY[x], collapse = "*"))), data = Y, freq = Freq), Y, freq = Freq)

summary(cvam(formula(paste("~", paste(namesY[x], collapse = "*"))), data = Y, freq = Freq))

tmpvec2 = sapply(cand.edges, function(x){
  get.loglik(cvam(formula(paste("~", paste(namesY[x], collapse = "*"))), data = Y, freq = Freq))
})

get.estimates(cvam(formula(paste("~", paste(namesY[x], collapse = "*"))), data = Y, freq = Freq))

cvamres = cvam(~ M * D, data = Y, freq = Freq)
lapply(cvamres, print)
cvamres
get.fitted(cvamres)[1:2]

Y

get.fitted(cvamres)[1:2] %in% c(1,2)
apply(Y, 1, function(y) which(y))

mapply(function, ...)

filter(cvamres, .[1:2] == Y[1,])

as.matrix(cvamPredict(~ M + P, cvam(~ M * P, data = Y, freq = Freq), data = Y))

cvamLik(~ M + P, cvam(~ M * P, data = Y, freq = Freq), data = Y)

cvamLik(~ M + P, cvam(~ M * P, data = Y, freq = Freq), data = Y)$likVal * Y$Freq

cvamLik(~ D, cvam(~ D, data = Y, freq = Freq), data = Y)$likVal
marginalProb$D

cvamLik(~ M, cvam(~ M, data = Y, freq = Freq), data = Y)$likVal
marginalProb$M

cvamLik(~ M, cvam(~ M, data = Y, freq = Freq), data = Y)$likVal

marginalProbmat =sapply(Y[-ncol(Y)], function(x) cvamLik(~ x, cvam(~ x, data = Y, freq = Freq), data = Y)$likVal)
marginalProbmat = as.data.frame(marginalProbmat)

mat2 = sapply(cand.edges, function(x){
  apply(select(marginalProbmat, -x), 1, prod) *
    cvamLik(formula(paste("~", paste(namesY[x], collapse = "+"))), cvam(formula(paste("~", paste(namesY[x], collapse = "*"))), data = Y, freq = Freq), data = Y)$likVal
})

###################
names(Y)

Z = cbind(Y, 1:nrow(Y))[3970, ]
tmpest = cvamEstimate(~ nativity , cvam(~ nativity, data = Y))
tmpest %>% filter(nativity == "born_other_state") %>% select(prob)

Z = cbind(Y, 1:nrow(Y))[3969, ]
tmpest = cvamEstimate(~ geog_mobility + ind_income | pov_status, cvam(~ geog_mobility * ind_income * pov_status, data = Y))
tmpest %>% filter(nativity == "born_other_state") %>% select(prob)


cvamEstimate(~ age | nativity, cvam(~ age * nativity, data = Y))$prob



cvamPredict(~ age , cvam(~ age * nativity, data = Y), data = Y[1,])

head(mat2)

x = c(1,2)
select(marginalProbmat, -x)
apply(select(marginalProbmat, -x), 1, prod)
cvamLik(formula(paste("~", paste(namesY[x], collapse = "+"))), cvam(formula(paste("~", paste(namesY[x], collapse = "*"))), data = Y, freq = Freq), data = Y)$likVal


as.matrix(cvamPredict(~ M + D, cvamres, data = Y, freq = Freq)) %*% get.fitted(cvamres, mfTrue = F)

as.matrix(cvamPredict(~ M + P, cvam(~ M * P, data = Y, freq = Freq), data = Y)) %*% get.fitted(cvamres, mfTrue = F)

as.matrix(cvamPredict(~ G + A, cvam(~ G * A, data = Y, freq = Freq), data = Y))
names(Y)
Y[2,]

cvamImpute(cvamres, Y, freq = Freq)

get.loglik(cvamres)

tmpres = glm(Freq ~ M  * D, data = Y, na.action = na.omit)
tmpres$fitted.values


glm(rep(1, nrow(Y)) ~ 0 + M  * D, data = Y, na.action = na.omit)
cvam(~ M * D, data = Y, freq = Freq)

