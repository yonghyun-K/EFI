
# summary(tract2221)
# data(tract2221, package = "imputeMulti"); p = ncol(tract2221)
# Y = tract2221
# for(k in 1:p)Y[[k]] <- droplevels(Y[[k]])
# summary(Y)
# edges_list <- construct_tree(p, ntree = 10)
# dp = doublep(Y, edges_list, R = 1)
# plot(dp)
# EFI = efi(tract2221, dp)
# estimate(EFI, "(marital_status == \"never_mar\") & (edu_attain == \"lt_hs\")")
