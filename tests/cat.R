# library(devtools)
# install_github("yonghyun-K/EFI", dependencies = T, force = T)
# library(EFI)

# older data example ####
data(older, package = "cat")
Y = tidyr::uncount(data.frame(older, stringsAsFactors= T), Freq); n = nrow(Y); p = ncol(Y)
for(k in 1:p) Y[[k]] = factor(Y[[k]])
summary(Y)
edges_list <- construct_tree(p)
dp = doublep(Y, edges_list, freq = F, R = 1)
plot(dp)
EFI = efi(Y, dp, freq = F)
miceres = mice::mice(Y, printFlag = F)
c(EFI = estimate(EFI, "(M == 1) & (P == 2)")$Estimate,
MICE = summary(mice::pool(mice:::with.mids(
  data = miceres, exp = lm((M == 1) & (P == 2) ~ 1))))$estimate)

# belt data example ####
data(belt, package = "cat")
Y = as.data.frame(belt); p = ncol(belt) - 1
for(k in 1:p) Y[[k]] = factor(Y[[k]])
edges_list = construct_tree(p)
dp = doublep(Y, edges_list, freq = T, R = 1)
plot(dp)
EFI = efi(Y, dp,  freq = T)
estimate(EFI, "(I1 == 1) & (I2 == 2)")$Estimate
estimate(EFI, "(I2 == 2) & (B2 == 2)")$Estimate

# cvam::cvamEstimate(~ I1 + I2, cvam::cvam(~ I1 * I2 + I1 * D + B2 * 0, data = Y, freq = Freq,
#                              control = cvam::cvamControl(iterMaxEM = 1500L))) # prob = 0.6839
# cvam::cvamEstimate(~ I2 + B2, cvam::cvam(~ I1 * I2 + I1 * D + B2 * 0, data = Y, freq = Freq,
#                              control = cvam::cvamControl(iterMaxEM = 1500L))) # prob = 0.5938
# cvam::cvamEstimate(~ I2 + B2, cvam::cvam(~ I1 * I2 + I1 * D + B2 * B1 + I1 * S + D * S + I2 * D + I1 * B1, data = Y, freq = Freq,
#                              control = cvam::cvamControl(iterMaxEM = 10000L))) # prob = 0.1698

Y = do.call("rbind", apply(belt, 1, function(x) matrix(rep(x[1:p], each = x[p+1]), nc = p)))
Y = as.data.frame(Y)
for(k in 1:p) Y[[k]] = factor(Y[[k]])
names(Y) <- colnames(belt)[-(p+1)]
# summary(lm(B2 == 1 ~ I1 + I2 + D + S+ B1, data = Y))

(n = nrow(Y)); sum(belt[,p+1])
miceres = mice::mice(Y)
summary(mice::pool(mice:::with.mids(data = miceres, exp = lm((I1 == 1) & (I2 == 2) ~ 1))))$estimate
summary(mice::pool(mice:::with.mids(data = miceres, exp = lm((I2 == 2) & (B2 == 2) ~ 1))))$estimate

# summary(glm(as.numeric(Y$I1) == 1 ~ as.numeric(Y$D), family = binomial(link = "identity"), weights = Y$Freq))
# summary(glm(as.numeric(Y$I1) == 1 ~ as.numeric(Y$S), family = binomial(link = "identity"), weights = Y$Freq))

# HairEyeColor data example ####
data(HairEyeColor)
Y = as.data.frame.table(HairEyeColor, stringsAsFactors = TRUE)
Y = Y_true = tidyr::uncount(Y, Freq)
n = nrow(Y); p = ncol(Y)
delta = matrix(rbinom(n * p, 1, 0.5), nr = n, nc = p); Y[delta == 0] = NA
edges_list = apply(combn(p, 2), 2, list)
dp = doublep(Y, edges_list, freq = F)
plot(dp)
EFI = efi(Y, dp, freq = F)
miceres = mice::mice(Y, printFlag = F)

c(True = mean(Y_true$Hair == "Black" & Y_true$Eye == "Brown"),
EFI = estimate(EFI, "(Hair == \"Black\") & (Eye == \"Brown\")")$Estimate,
MICE = summary(mice::pool(mice:::with.mids(data = miceres, exp = lm((
  Hair == "Black") & (Eye == "Brown") ~ 1))))$estimate)

c(True = mean(Y_true$Hair == "Black" & Y_true$Sex == "Male"),
EFI = estimate(EFI, "(Hair == \"Black\") & (Sex == \"Male\")")$Estimate,
MICE = summary(mice::pool(mice:::with.mids(data = miceres, exp = lm((
  Hair == "Black") & (Sex == "Male") ~ 1))))$estimate)
