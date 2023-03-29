library(foreach)
library(doParallel)
library(CVXR)
library(cat)
# library(mice)
library(igraph)
library(plyr)
library(dplyr)
library(cvam)

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



data(tract2221, package = "imputeMulti")
Y = tract2221
p = ncol(tract2221)
cand.edges = as.list(data.frame(combn(p, 2)))
dp = doublep(tract2221, cand.edges)
plot(dp)
EFI = efi(tract2221, dp)
names(tract2221)
lapply(tract2221, levels)
summary(tract2221)
estimate(EFI, "(marital_status == \"never_mar\") & (edu_attain == \"lt_hs\")")
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

cvamPredict(~ edu_attain , cvam(~ age * edu_attain, data = Y), data = Y[1,])


cvamPredict(~ age + edu_attain , cvam(~ age * edu_attain, data = Y), data = Y[2,])

cvamImpute(cvam(~ age * edu_attain, data = Y), data= Y)

cvamImpute(cvam(~ age * edu_attain, data = Y), data= Y, synthetic = T)

cvamEstimate(~ edu_attain | age, cvam(~ age * edu_attain, data = Y))

cvamEstimate(~ age, cvam(~ age * edu_attain, data = Y))

cvam(~, data = Y)

names(Y)

head(Y)
tail(Y)

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




weightvec = dp$weightvec
weightmat = dp$weightmat
cand.edges = dp$cand.edges
supp = dp$supp
n = dp$n
p = dp$p
marginalProb = dp$marginalProb

Y_FI0 = adply(cbind(Y, 1:nrow(Y)), 1, function(Z){
  y = Z[-c(p+1, p+2)]
  Freq = Z[p+1]
  id = Z[p+2]
  mis_idx = which(is.na(y))
  if(length(mis_idx) == 0){
    return(cbind(y, Freq, id = id))
  }else{
    obs_idx0 = (1:p)[-mis_idx]
    obs_idx = which(colSums(weightmat[mis_idx,, drop = F]) != 0)
    obs_idx = obs_idx[!(obs_idx %in% mis_idx)]

    grid_mis = expand.grid(supp[mis_idx])
    mat_FI = adply(grid_mis, 1, function(z){
      # mis_val = unlist(z)
      # obs_val0 = unlist(y[obs_idx0])
      # obs_val = unlist(y[obs_idx])
      mis_val = z
      obs_val0 = y[obs_idx0]
      obs_val = y[obs_idx]
      y_tmp = y
      y_tmp[mis_idx] <- mis_val; y_tmp[obs_idx0] <- obs_val0;
      # # print(y_tmp)
      #
      #
      # sel_cand = sapply(cand.edges, function(x) all(x %in% c(mis_idx, obs_idx)))
      #
      # # print(grid_tmp[,c(mis_idx, obs_idx),drop = F])
      # # print(z)
      # # print(cbind(mis_val, obs_val))
      # # print(apply(grid_tmp[,c(mis_idx, obs_idx),drop = F], 1, function(x) all(x == cbind(mis_val, obs_val))) )
      # # condP = colSums(mat[apply(t(grid_tmp[,c(mis_idx, obs_idx),drop = F]) == c(mis_val, obs_val), 2, all), sel_cand ,drop = F]) / colSums(mat[apply(t(grid_tmp[,obs_idx,drop = F]) == obs_val, 2, all), sel_cand, drop = F]) # P(Y2 = 1, Y4 = 2 | Y1 = 1, Y2 = 2, Y3 = 2)
      # condP = colSums(mat[apply(grid_tmp[,c(mis_idx, obs_idx),drop = F], 1, function(x) all(x == cbind(mis_val, obs_val))), sel_cand ,drop = F]) / colSums(mat[apply(grid_tmp[,obs_idx,drop = F], 1, function(x) all(x == obs_val)), sel_cand, drop = F]) # P(Y2 = 1, Y4 = 2 | Y1 = 1, Y2 = 2, Y3 = 2)
      #
      # # print(paste("P(", paste(paste("Y", mis_idx, sep = "") , mis_val, sep = "=", collapse = ", "), "|", paste(paste("Y", obs_idx, sep = "") , obs_val, sep = "=", collapse = ", "), ")"));
      # sumtmp = sum(weightvec[sel_cand])
      # if(sumtmp == 0){
      #   weightvec_tmp = rep(1, length(weightvec[sel_cand])) / length(weightvec[sel_cand])
      # }else{
      #   weightvec_tmp = weightvec[sel_cand] / sumtmp
      # }
      #
      # if(length(condP) != 0){
      #   fweight = sum(condP * weightvec_tmp)
      # }else{
      #   # fweight = marginalProb[[mis_idx]][mis_val]
      #   fweight = sapply(marginalProb[mis_idx], function(prob_tmp) prob_tmp[mis_val])
      # }
      #
      # # print(fweight)
      # # print(cbind(y_tmp, fweight * Freq, id))
      cbind(y_tmp)
    })
    return(mat_FI)
  }
})

nrow(Y_FI0)

Z = Y[1,,drop = F]

Y[integer(0),]

Y
data.frame()

Ztmp[,is.na(Ztmp)] = expand.grid(lapply(Ztmp[,is.na(Ztmp)], levels))
Ztmp[,]

cbind(expand.grid(lapply(Ztmp[,is.na(Ztmp)], levels)), Ztmp[,!is.na(Ztmp)])

expand.grid(lapply(Ztmp[,is.na(Ztmp)], levels))
Ztmp[,!is.na(Ztmp)]


cbind.data.frame(NULL, Ztmp[,!is.na(Ztmp)])

data.frame(Ztmp[,!is.na(Ztmp)])

apply(tmpY, 1, function(y) as.numeric(y))

lapply(marginalProb, "[[", as.numeric(tmpY[1,]))

mapply(function(x, y)x[y], marginalProb, as.numeric(tmpY[1,]))

mice::mice(Y)

nrow(Y_FI0)
Y_FI0[Y_FI0$`1:nrow(Y)` == 3,]
tmpY[names(Y)]

namestmp = names(cbind(Y, id = 1:nrow(Y)))
Y_FI0 = adply(cbind(Y, id = 1:nrow(Y)), 1, function(Z){
  if(any(is.na(Z))){
    tmpY = cbind(expand.grid(lapply(Z[,is.na(Z),drop = F], levels)), Z[,!is.na(Z)])[namestmp]
    return(tmpY)
  }else{
    Z
  }
})

marginalProbmat = sapply(Y[-ncol(Y)], function(x) cvamLik(~ x, cvam(~ x, data = Y, freq = Freq), data = Y)$likVal)
marginalProbmat = as.data.frame(marginalProbmat)
tmpmat = sapply(cand.edges[weightvec != 0], function(x){
  apply(select(marginalProbmat, -x), 1, prod) *
    cvamLik(formula(paste("~", paste(namesY[x], collapse = "+"))), cvam(formula(paste("~", paste(namesY[x], collapse = "*"))), data = Y, freq = Freq), data = Y)$likVal
})
tmpmat = as.data.frame(tmpmat)

marginalProbmat2 =sapply(names(Y)[-ncol(Y)], function(x){
  form = as.formula(paste("~", x))
  cvamLik(form, cvam(form, data = Y, freq = Freq), data = Y_FI0)$likVal
})
marginalProbmat2 = as.data.frame(marginalProbmat2)

tmpmat2 = sapply(cand.edges[weightvec != 0], function(x){
  apply(select(marginalProbmat2, -x), 1, prod) *
  cvamLik(formula(paste("~", paste(namesY[x], collapse = "+"))), cvam(formula(paste("~", paste(namesY[x], collapse = "*"))), data = Y, freq = Freq), data = Y_FI0)$likVal
})
tmpmat2 = as.data.frame(tmpmat2)

timestmp = unlist(table(Y_FI0$id))
restmp = mapply(function(x, y) return(x / rep(y, times = timestmp)), tmpmat2, tmpmat)

fw = apply(restmp, 1, function(x) sum(x * weightvec))
cbind(Y_FI0, fw)


