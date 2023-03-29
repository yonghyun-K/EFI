# library(imputeMulti)

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
