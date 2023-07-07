# library(cvam)
# library(EFI)
# library(dplyr)

data(crime, package = "cvam")
Y = tidyr::uncount(data.frame(crime, stringsAsFactors= T), n); n = nrow(Y); p = ncol(Y)
summary(Y)
edges_list <- construct_tree(p)
dp = doublep(Y, edges_list, freq = F, R = 1)
plot(dp)
EFI = efi(Y, dp, freq = F)
miceres = mice::mice(Y, printFlag = F)
c(EFI = estimate(EFI, "(V1 == \"no\") & (V2 == \"no\")")$Estimate,
  MICE = summary(mice::pool(mice:::with.mids(
    data = miceres, exp = lm((V1 == "no") & (V2 == "no") ~ 1))))$estimate)
c(EFI = estimate(EFI, "(V1 == \"yes\") & (V2 == \"yes\")")$Estimate,
  MICE = summary(mice::pool(mice:::with.mids(
    data = miceres, exp = lm((V1 == "yes") & (V2 == "yes") ~ 1))))$estimate)

data(hivtest, package = "cvam")
Y = Y_true = tidyr::uncount(hivtest, COUNT); n = nrow(Y); p = ncol(Y)
summary(Y)
delta = matrix(rbinom(n * p, 1, 0.5), nr = n, nc = p)
Y[delta == 0] <- NA
edges_list <- construct_tree(p)
dp = doublep(Y, edges_list, freq = F, R = 1)
plot(dp)
EFI = efi(Y, dp, freq = F)
miceres = mice::mice(Y, printFlag = F)
c(True = mean(Y_true$A == "neg" & Y_true$B == "neg"),
  EFI = estimate(EFI, "(A == \"neg\") & (C == \"neg\")")$Estimate,
  MICE = summary(mice::pool(mice:::with.mids(
    data = miceres, exp = lm((A == "neg") & (C == "neg") ~ 1))))$estimate)

data(microUCBAdmissions, package = "cvam")
Y = Y_true = microUCBAdmissions; n = nrow(Y); p = ncol(Y)
summary(Y)
delta = matrix(rbinom(n * p, 1, 0.5), nr = n, nc = p)
Y[delta == 0] <- NA
edges_list <- construct_tree(p)
dp = doublep(Y, edges_list, freq = F, R = 1)
plot(dp)
EFI = efi(Y, dp, freq = F)
miceres = mice::mice(Y, printFlag = F)
c(True = mean(Y_true$Admit == "Admitted" & Y_true$Dept == "A"),
  EFI = estimate(EFI, "(Admit == \"Admitted\") & (Dept == \"A\")")$Estimate,
  MICE = summary(mice::pool(mice:::with.mids(
    data = miceres, exp = lm((Admit == "Admitted") & (Dept == "A") ~ 1))))$estimate)


data(seatbelt, package = "cvam")
Y = tidyr::uncount(data.frame(seatbelt, stringsAsFactors= T), n); n = nrow(Y); p = ncol(Y)
summary(Y)
edges_list <- construct_tree(p)
dp = doublep(seatbelt, edges_list, freq = T, R = 1)
plot(dp)
EFI = efi(seatbelt, dp, freq = T)
miceres = mice::mice(Y, printFlag = T)
c(EFI = estimate(EFI, "(injury.p == \"no\") & (injury.f == \"no\")")$Estimate,
  MICE = summary(mice::pool(mice:::with.mids(
    data = miceres, exp = lm((injury.p == "no") & (injury.f == "no") ~ 1))))$estimate)

data(cig2019, package = "cvam")
Y = Y0 = cig2019[sapply(cig2019, is.factor)]; n = nrow(Y); p = ncol(Y)
Y = dplyr::count(Y, !!!Y)
summary(Y)
edges_list <- construct_tree(p)
dp = doublep(Y, edges_list, freq = T, R = 1)
plot(dp)
EFI = efi(Y, dp, freq = T)
miceres = mice::mice(Y0, printFlag = T)
c(EFI = estimate(EFI, "(sex_a == \"Male\") & (smkev_a == \"Yes\")")$Estimate,
  MICE = summary(mice::pool(mice:::with.mids(
    data = miceres, exp = lm((sex_a == "Male") & (smkev_a == "Yes") ~ 1))))$estimate)


# data(abortion2000, package = "cvam")
# Y = abortion2000[sapply(abortion2000, is.factor)]; n = nrow(Y); p = ncol(Y)
# # table(dplyr::count(Y, !!!Y)$n)
# summary(Y)
# edges_list <- construct_tree(p)
# dp = doublep(Y, edges_list, freq = F, R = 1)
# plot(dp)
# EFI = efi(Y, dp, freq = F)
# miceres = mice::mice(Y, printFlag = T)
# c(EFI = estimate(EFI, "(AbNoMore == \"No\") & (AbSingle == \"Yes\")")$Estimate,
#   MICE = summary(mice::pool(mice:::with.mids(
#     data = miceres, exp = lm((AbNoMore == "No") & (AbSingle == "Yes") ~ 1))))$estimate)

