library(devtools)
install_github("yonghyun-K/EFI", dependencies = T, force = T)
library(EFI)

library(cvam)

data(older, package = "cat")
p = 6
edges_list = apply(combn(p, 2), 2, list)
Y = data.frame(older, stringsAsFactors= T)
for(k in 1:(ncol(Y) - 1)){
  Y[[k]] = factor(Y[[k]])
}

dp = doublep(Y, edges_list, freq = T)
plot(dp)
EFI = efi(Y, dp, freq = T)
estimate(EFI, "(M == 1) & (P == 2)")

data(belt, package = "cat")
belt
p = 6
Y = as.data.frame(belt)
for(k in 1:p){
  Y[[k]] = factor(Y[[k]])
}

edges_list = apply(combn(p, 2), 2, list)
dp = doublep(Y, edges_list, freq = T)
plot(dp)
EFI = efi(Y, dp,  freq = T)
estimate(EFI, "(I1 == 1) & (I2 == 2)")
estimate(EFI, "(I2 == 2) & (B2 == 2)")

cvamEstimate(~ I1 + I2, cvam(~ I1 * I2 + I1 * D + B2 * 0, data = Y, freq = Freq,
                             control = cvamControl(iterMaxEM = 1500L))) # prob = 0.6839
cvamEstimate(~ I2 + B2, cvam(~ I1 * I2 + I1 * D + B2 * 0, data = Y, freq = Freq,
                             control = cvamControl(iterMaxEM = 1500L))) # prob = 0.5938
cvamEstimate(~ I2 + B2, cvam(~ I1 * I2 + I1 * D + B2 * B1 + I1 * S + D * S + I2 * D + I1 * B1, data = Y, freq = Freq,
                             control = cvamControl(iterMaxEM = 10000L))) # prob = 0.1698

# summary(glm(B2 == 1 ~ I1 + I2 + D + S + B1, family = binomial("identity"), weights = Freq, data = data.frame(belt)))

p = 6
Y = do.call("rbind", apply(belt, 1, function(x) matrix(rep(x[1:p], each = x[p+1]), nc = p)))
Y = as.data.frame(Y)
for(k in 1:p){
  Y[[k]] = factor(Y[[k]])
}
names(Y) <- colnames(belt)[-(p+1)]
summary(lm(B2 == 1 ~ I1 + I2 + D + S+ B1, data = Y))

(n = nrow(Y)); sum(belt[,p+1])
miceres = mice::mice(Y)
micefit <-mice:::with.mids(data = miceres, exp = lm((I1 == 1) & (I2 == 2) ~ 1))
summary(mice::pool(micefit))
micefit <-mice:::with.mids(data = miceres, exp = lm((I2 == 2) & (B2 == 2) ~ 1))
summary(mice::pool(micefit))

# summary(glm(as.numeric(Y$I1) == 1 ~ as.numeric(Y$D), family = binomial(link = "identity"), weights = Y$Freq))
# summary(glm(as.numeric(Y$I1) == 1 ~ as.numeric(Y$S), family = binomial(link = "identity"), weights = Y$Freq))

data(HairEyeColor)
n = nrow(HairEyeColor)
p = 3
Y = do.call("rbind", apply(as.data.frame.table(HairEyeColor), 1, function(x) matrix(rep(x[1:p], each = x[p+1]), nc = p)))
Y = as.data.frame(Y)
names(Y) <- names(dimnames(HairEyeColor))
mean(Y$Hair == "Black" & Y$Eye == "Brown")
mean(Y$Hair == "Black" & Y$Sex == "Male")
delta = matrix(rbinom(n * p, 1, 0.5), nr = n, nc = p)
for(k in 1:p){
  Y[[k]] = factor(Y[[k]])
  Y[delta[,k] == 0, k] <- NA
}
# (n = nrow(Y)); sum(HairEyeColor)

edges_list = apply(combn(p, 2), 2, list)
dp = doublep(Y, edges_list, freq = F)
plot(dp)
EFI = efi(Y, dp, freq = F)
# data = EFI
estimate(EFI, "(Hair == \"Black\") & (Eye == \"Brown\")")
estimate(EFI, "(Hair == \"Black\") & (Sex == \"Male\")")

Y_true = do.call("rbind", apply(as.data.frame.table(HairEyeColor), 1, function(x) matrix(rep(x[1:p], each = x[p+1]), nc = p)))
Y_true = as.data.frame(Y_true)
for(k in 1:p){
  Y_true[[k]] = factor(Y_true[[k]])
}
names(Y_true) <- names(dimnames(HairEyeColor))
mean(Y_true$Hair == "Black" & Y_true$Eye == "Brown")
mean(Y_true$Hair == "Black" & Y_true$Sex == "Male")

miceres = mice::mice(Y)
micefit <-mice:::with.mids(data = miceres, exp = lm((Hair == "Black") & (Eye == "Brown") ~ 1))
summary(mice::pool(micefit))
micefit <-mice:::with.mids(data = miceres, exp = lm((Hair == "Black") & (Sex == "Male") ~ 1))
summary(mice::pool(micefit))
