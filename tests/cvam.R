# library(cvam)
library(EFI)

data(crime, package = "cvam")
Y = crime
p = 2
cand.edges = as.list(data.frame(combn(p, 2)))
dp = doublep(Y, cand.edges, freq = T)
plot(dp)
EFI = efi(Y, dp, freq = T)
names(Y)
lapply(Y, levels)
summary(Y)

library(dplyr)
estimate(EFI, "(V1 == \"no\") & (V2 == \"no\")")
estimate(EFI, "(V1 == \"yes\") & (V2 == \"yes\")")

sum(EFI$imp[3:4,5])

library(cvam)

cvamEstimate(~ V1 + V2, cvam(~V1 * V2, data = Y, freq = n), data = Y)
cvamEstimate(~ V1 + V2, cvam(~V1 * V2, data = Y, freq = n), data = Y)$SE - estimate(EFI, "(V1 == \"no\") & (V2 == \"no\")")$sqrtB

EFI$imp

Y[1:8,3] <- 1

Y[3, 3] <- 1000

p = 2
Y2 = do.call("rbind", apply(Y, 1, function(x) matrix(rep(x[1:p], each = x[p+1]), nc = p)))
Y2 = as.data.frame(Y2)
for(k in 1:p){
  Y2[[k]] = factor(Y2[[k]])
}
names(Y2) <- colnames(crime)[1:2]
miceres = mice::mice(Y2, printFlag = F, m = 30)

mice::complete(miceres, action = 2)
mean(mice::complete(miceres, action = 4)[,1] == "yes")




cor(mice::complete(miceres, action = 1) == "no", mice::complete(miceres, action = 4) == "no")

var(mice::complete(miceres, action = 2)[,1] == "no" & mice::complete(miceres, action = 2)[,2] == "no")

miceres = mice::mice(Y2, printFlag = F, m = 30)
tmp2 = sapply(1:30, function(x) mean(apply(mice::complete(miceres, action = x), 1, function(y) all(y == c("no", "no")) )))
plot(tmp2)
points(sapply(1:30, function(x) mean(rbinom(4, 1, 0.25))), col = "red")

# var(tmp2)
# mean(tmp2)
Y

sqrt(0.25 * 0.75 / 4)

sqrt(0.5 * 0.5 / 6)

micefit <-mice:::with.mids(data = miceres, exp = lm((V1 == "no") & (V2 == "no") ~ 1))
micefit <-mice:::with.mids(data = miceres, exp = lm((V1 == "no") ~ 1))
summary(mice::pool(micefit))
tmp <- mice::pool(micefit)
sqrt(tmp$pooled$t)
sqrt(tmp$pooled$b)
sqrt(tmp$pooled$ubar)

sqrt(0.25 * 0.75 / 28.75)

sqrt(tmp$pooled$estimate * (1 - tmp$pooled$estimate) / nrow(Y2))

dp = doublep(Y2, cand.edges, freq = F)
EFI = efi(Y2, dp, freq = F)
names(Y)
lapply(Y, levels)
summary(Y)

estimate(EFI, "(V1 == \"no\") & (V2 == \"no\")")
estimate(EFI, "(V1 == \"yes\") & (V2 == \"yes\")")


data(hivtest, package = "cvam")

Y = hivtest
p = 4
cand.edges = as.list(data.frame(combn(p, 2)))
dp = doublep(Y, cand.edges, freq = T)
plot(dp)
EFI = efi(Y, dp, freq = T)
names(Y)
lapply(Y, levels)
summary(Y)

data(microUCBAdmissions, package = "cvam")
Y = microUCBAdmissions
p = ncol(Y)
cand.edges = as.list(data.frame(combn(p, 2)))
dp = doublep(Y, cand.edges)
plot(dp)
EFI = efi(Y, dp)
names(Y)
lapply(Y, levels)
summary(Y)

data(seatbelt, package = "cvam")
Y = seatbelt[2:8] # Error in Y = seatbelt
p = 6
cand.edges = as.list(data.frame(combn(p, 2)))
dp = doublep(Y, cand.edges, freq = T)
plot(dp)
EFI = efi(Y, dp, freq = T)
names(Y)
lapply(Y, levels)
summary(Y)

data(cig2019, package = "cvam")
Y = cig2019[sapply(cig2019, is.factor)]
Y = Y[-5] # Error in Y = cig2019[sapply(cig2019, is.factor)]
p = ncol(Y)
cand.edges = as.list(data.frame(combn(p, 2)))
dp = doublep(Y, cand.edges)
plot(dp)
EFI = efi(Y, dp)
names(Y)
lapply(Y, levels)
summary(Y)

data(abortion2000, package = "cvam")
Y = abortion2000[sapply(abortion2000, is.factor)]
p = ncol(Y)
cand.edges = as.list(data.frame(combn(p, 2)))
dp = doublep(Y, cand.edges)
plot(dp)
EFI = efi(Y, dp)
names(Y)
lapply(Y, levels)
summary(Y)
estimate(EFI, "(Race == \"White\") & (CenRace == \"White\")")
