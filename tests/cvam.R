# library(cvam)

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
