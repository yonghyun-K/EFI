library(devtools)
install_github("yonghyun-K/EFI", dependencies = T, force = T)

library(EFI)
library(dplyr)

data(older, package = "cat")
# p = 6
# Y = do.call("rbind", apply(older, 1, function(x) matrix(rep(x[1:p], each = x[p+1]), nc = p)))
# Y = as.data.frame(Y)
# for(k in 1:p){
#   Y[[k]] = factor(Y[[k]])
# }
# # p = ncol(Y)
# (n = nrow(Y)); sum(older[,p+1])
#
# cand.edges = as.list(data.frame(combn(p, 2)))
# dp = doublep(Y, cand.edges)
# plot(dp)
# EFI = efi(Y, dp)
# estimate(EFI, "(V1 == 1) & (V2 == 2)")

p = 6
cand.edges = as.list(data.frame(combn(p, 2)))
Y = data.frame(older, stringsAsFactors= T)
for(k in 1:(ncol(Y) - 1)){
  Y[[k]] = factor(Y[[k]])
}
dp = doublep(Y, cand.edges, freq = T)
plot(dp)
EFI = efi(Y, dp, freq = T)
estimate(EFI, "(M == 1) & (P == 2)")

# data(crimes, package = "cat")
#
# p = 2
# Y = as.data.frame(crimes)
# for(k in 1:p){
#   Y[[k]] = factor(Y[[k]])
# }
#
# cand.edges = as.list(data.frame(combn(p, 2)))
# dp = doublep(Y, cand.edges, freq = T)
# EFI = efi(Y, dp)
# estimate(EFI, "(V1 == 1) & (V2 == 2)")

data(belt, package = "cat")
belt
p = 6
Y = as.data.frame(belt)
for(k in 1:p){
  Y[[k]] = factor(Y[[k]])
}

cand.edges = as.list(data.frame(combn(p, 2)))
dp = doublep(Y, cand.edges, freq = T)
plot(dp)
EFI = efi(Y, dp,  freq = T)
estimate(EFI, "(I1 == 1) & (I2 == 2)")
estimate(EFI, "(I2 == 2) & (B2 == 2)")

# summary(glm(as.numeric(Y$I1) == 1 ~ as.numeric(Y$D), family = binomial(link = "identity"), weights = Y$Freq))
# summary(glm(as.numeric(Y$I1) == 1 ~ as.numeric(Y$S), family = binomial(link = "identity"), weights = Y$Freq))

data(HairEyeColor, package = "cat")

p = 3
Y = do.call("rbind", apply(as.data.frame.table(HairEyeColor), 1, function(x) matrix(rep(x[1:p], each = x[p+1]), nc = p)))
Y = as.data.frame(Y)
for(k in 1:p){
  Y[[k]] = factor(Y[[k]])
  # levels(Y[[k]]) <- as.character(1:nlevels(Y[[k]]))
}
(n = nrow(Y)); sum(HairEyeColor)
delta = matrix(rbinom(n * p, 1, 0.9), nr = n, nc = p)
Y[delta == 0] = NA

cand.edges = as.list(data.frame(combn(p, 2)))
dp = doublep(Y, cand.edges, freq = F)
plot(dp)
EFI = efi(Y, dp, freq = F)
# data = EFI
estimate(EFI, "(V1 == \"Black\") & (V2 == \"Brown\")")
estimate(EFI, "(V1 == \"Black\") & (V3 == \"Male\")")
