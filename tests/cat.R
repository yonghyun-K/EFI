library(dplyr)

data(older, package = "cat")
p = 6
Y = do.call("rbind", apply(older, 1, function(x) matrix(rep(x[1:p], each = x[p+1]), nc = p)))
Y = as.data.frame(Y)
for(k in 1:p){
  Y[[k]] = factor(Y[[k]])
}
# p = ncol(Y)

(n = nrow(Y)); sum(older[,p+1])

cand.edges = as.list(data.frame(combn(p, 2)))
dp = doublep(Y, cand.edges)
EFI = efi(Y, dp)
estimate(EFI, "(V1 == 1) & (V2 == 2)")
# with(imp, summary(glm((Var1 == 1) & (Var2 == 2) ~ 1, family = binomial(link = "identity"), weights = w)))


# data(crimes, package = "cat")
# 
# p = 2
# Y = do.call("rbind", apply(crimes, 1, function(x) matrix(rep(x[1:p], each = x[p+1]), nc = p)))
# Y = as.data.frame(Y)
# for(k in 1:p){
#   Y[[k]] = factor(Y[[k]])
# }
# 
# (n = nrow(Y)); sum(crimes[,p+1])
# 
# cand.edges = as.list(data.frame(combn(p, 2)))
# dp = doublep(Y, cand.edges)
# EFI = efi(Y, dp)
# estimate(EFI, "(V1 == 1) & (V2 == 2)")

data(belt, package = "cat")
belt
p = 6
Y = do.call("rbind", apply(belt, 1, function(x) matrix(rep(x[1:p], each = x[p+1]), nc = p)))
Y = as.data.frame(Y)
for(k in 1:p){
  Y[[k]] = factor(Y[[k]])
}

(n = nrow(Y)); sum(belt[,p+1])
Y = Y[sample(1:n, 1000),]

cand.edges = as.list(data.frame(combn(p, 2)))
dp = doublep(Y, cand.edges)
plot(dp)
EFI = efi(Y, dp)
estimate(EFI, "(V1 == 1) & (V2 == 2)")


data(HairEyeColor, package = "cat")

p = 3
Y = do.call("rbind", apply(as.data.frame.table(HairEyeColor), 1, function(x) matrix(rep(x[1:p], each = x[p+1]), nc = p)))
Y = as.data.frame(Y)
for(k in 1:p){
  Y[[k]] = factor(Y[[k]])
  levels(Y[[k]]) <- as.character(1:nlevels(Y[[k]]))
}
(n = nrow(Y)); sum(HairEyeColor)
delta = matrix(rbinom(n * p, 1, 0.9), nr = n, nc = p)
Y[delta == 0] = NA

cand.edges = as.list(data.frame(combn(p, 2)))
dp = doublep(Y, cand.edges)
plot(dp)
EFI = efi(Y, dp)
estimate(EFI, "(V1 == 1) & (V2 == 2)")
