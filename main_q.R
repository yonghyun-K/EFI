# args <- commandArgs(trailingOnly = TRUE)

# args <- c(n, n_B, B, SIMNUM, lambda)
# args <- c(20000, 15000, 5, 12, 100, 30, 2)

library(foreach)
library(doParallel)
library(xtable)
library(doSNOW)
library(mice, warn.conflicts = FALSE)

timenow0 = Sys.time()
timenow = paste(timenow0, ".txt", sep = "")

writeLines(c(""), timenow)
sink(timenow, append=TRUE)

source("par.R")

# p_Y2 = p_Y
# p_Y2 = c((0.8 + p_Y1[1]) / 2, 0.8 * p_Y1[2] + 0.1, 0.9 - 0.5 * p_Y1[3], 0.6)
# print(paste("p_Y2 =", paste(p_Y2, collapse = " ")))

# theta2 = sum(p_x * p_Y2)
# print(paste("theta2 =", theta2))

# theta = c(theta1, theta2)

res = foreach(simnum = 1:SIMNUM,
              .packages = c("mice"),
              .options.snow=opts) %dopar% {
                # for(simnum in 1:SIMNUM){
                
                if(simnum == 1) sink(timenow, append=TRUE)
                print(paste("simnum =", simnum))
                
                # Generate X and Y ####
                
                # is.factor(X[,1]) == FALSE: X can be high-dimensional
                X_num = matrix(rep(1:k_x, n * p)[rmultinom(n * p, 1, p_x) == 1], nr = n, nc = p) 
                colnames(X_num) <- paste("X", 1:p, sep = "")
                
                p_Y_mat = apply(X_num[,1:q, drop = F], 2, function(k) p_Y[k])
                # p_Y_mat = apply((X_num[,1:q, drop = F] + X_num[,(1:q) * 2, drop = F] + X_num[,(1:q) * 3, drop = F]), 2, function(k) p_Y[k])
                colnames(p_Y_mat) <- paste("Y", 1:q, sep = "")
                # p_Y_mat[,2] = rowSums(cbind((0.8 + Y1) / 2, 0.8 * Y1 + 0.1, 0.9 - 0.5 * Y1, 0.6) * model.matrix(~ 0 + as.factor(X[,1])))
                
                Y_num = apply(p_Y_mat, 2, function(k) rbinom(nrow(p_Y_mat), 1, k))
                
                # if(p > 20){
                #   print("Y_num[,1] = ifelse(rowsum(X_num[,1:10]) %% 2 == 0, 0, 1)
                #   Y_num[,2] = ifelse(apply(X_num[,seq(from = 11, to = 19, by = 2)] - X_num[,seq(from = 12, to = 20, by = 2)], 1, prod)
                #                      > 0, 0, 1)")
                #   Y_num[,1] = ifelse(rowSums(X_num[,1:10]) %% 2 == 0, 0, 1)
                #   Y_num[,2] = ifelse(apply(X_num[,seq(from = 11, to = 19, by = 2)] - X_num[,seq(from = 12, to = 20, by = 2)], 1, prod)
                #                      != 0, 0, 1)
                # }

                
                Y = data.frame(Y_num)
                for(k in 1:q){
                  Y[[k]] = factor(Y[[k]])
                }
                # Y = structure(factor(Y_num), dim = dim(Y_num), class = c('matrix', 'factor'))
                # Y = apply(X[,1:q], 2, function(k) rbinom(n, 1, p_Y[k])))
                
                colnames(Y) <- paste("Y", 1:q, sep = "")

                p_delta = sapply(1:q, p_delta_ftn)
                
                print("p_delta_ftn")
                print(p_delta_ftn)
                print(paste("mean(p_delta)", round(mean(p_delta), 5)))
                
                print(summary(p_delta))
                
                # p_delta = sapply(1:q, )

                # p_delta = sapply(1:q,  ) # 0.5
                
                # p_delta[,q] = 1
                # colMeans(p_delta)
                
                delta_num = apply(p_delta, 2, function(k) rbinom(nrow(p_delta), 1, k))
                delta = data.frame(delta_num)
                for(k in 1:q){
                  delta[[k]] = factor(delta[[k]])
                }
                # delta = structure(factor(delta_num), dim = dim(delta_num), class = c('matrix', 'factor'))
                colnames(delta) <- paste("delta", 1:q, sep = "")
                
                Y_ogn = Y
                
                theta_full = colMeans(Y_num)
                
                # is.factor(Y[,1]) == TRUE: Y is low-dimensional and saved as factors.
                Y[delta == 0] = NA
                
                Y_num[delta == 0] = NA
                
                # train_idx = sample(1:n, round(n * 0.9), replace = FALSE)
                X = data.frame(X_num)
                for(k in 1:p){
                  X[[k]] = factor(X[[k]])
                }
                
                train_idx = 1:n
                x = X[train_idx,]
                y = Y[train_idx,, drop = F]
                delta_train = delta[train_idx,, drop = F]
                x_num = X_num[train_idx,]
                y_num = Y_num[train_idx,, drop = F]
                
                # MICE ####
                imp <- mice(cbind(X, Y_num), printFlag = FALSE)
                theta_mice = sapply(1:q, function(k){
                  str = paste("Y", k, "~1", collapse = "", sep = "")
                  summary(pool(with(imp, lm(formula(str)))))$estimate
                })
                
                comp_mice = complete(imp,1:5)
                X_mice = comp_mice[,1:p]
                y_mice = comp_mice[,(p+1):ncol(comp_mice)]

                w_B = NULL
                # bstp_idx_B = NULL
                
                for(b in 1:B){
                  print(paste("b =", b))
                  # select_x = sample(1:p, p_star, replace = F)
                  
                  cand_tmp = ncol(combn(5, p_star))
                   select_x = combn(5, p_star)[,(b+cand_tmp-1) %% cand_tmp + 1]
                  # select_x = c(1, 2)
                  # select_x = c(3, 4)
                  # select_x = c(5, 6)

                  # table(data.frame(cbind(X[,select_x], Y_ogn)), useNA = "ifany")
                  
                  # n_mat_true = table(data.frame(cbind(X_mice[,select_x], y_mice, delta)), useNA = "ifany"); if(b == 1) print("use MICE for initial values")
                  n_mat_true = table(data.frame(cbind(X[,select_x], Y_ogn, delta)), useNA = "ifany"); if(b == 1) print("use full data for initial values")
                  # summary(data.frame(cbind(X[,select_x], Y_ogn, delta)))

                  expand_txt = paste(paste(rep("c(1,0)", q)), collapse = ",")
                  expand_txt = paste("expand.grid(", expand_txt, ")", collapse = "", sep = "")
                  delind = eval(parse(text = expand_txt))
                  
                  p_01s_true = vector(mode = "list", length = nrow(delind)) # P(Y_mis | x1, x2, Y_obs, delta)
                  for(k in 1:nrow(delind)){
                    n_mat_trueTx = paste("n_mat_true[", paste(rep(",", p_star + q), collapse = ""), paste(delind[k,] + 1, collapse = ","), "]")
                    MARGINTx = paste(c(1:p_star, ((p_star + 1):(p_star + q))[delind[k,] == 1]), collapse = ", ")
                    MARGINTx = paste("MARGIN = c(", MARGINTx, ")")
                    applyTx = paste(n_mat_trueTx, MARGINTx, "sum", sep = ", ")
                    applyTx = paste("apply(", applyTx, ")")
                    p_01s_trueTx = paste(n_mat_trueTx, MARGINTx, applyTx,"\"/\"", sep = ", ")
                    p_01s_trueTx = paste("sweep(", p_01s_trueTx, ")")
                    # print(p_01s_trueTx)
                    tmp <- eval(parse(text = p_01s_trueTx))
                    tmp[is.na(tmp)] <- 1/(2^sum(delind[k,] == 0))
                    
                    # tmp[!is.na(tmp)] <- 1/(2^sum(delind[k,] == 0)); if(b == 1 & k == 1) print("use uniform priors for initial") # to be removed
                    
                    p_01s_true[[k]] <- tmp
                  }

                  bstp_idx = sample(1:length(train_idx), n_B, replace = F)
                  x_b = x[bstp_idx,select_x, drop = F]
                  y_b = y[bstp_idx,, drop = F]
                  z_b = cbind(x_num[bstp_idx,select_x], y_num[bstp_idx,, drop = F] + 1)
                  delta_b = delta_train[bstp_idx,, drop = F]
                  
                  x_oob = x_num[!(1:length(train_idx) %in% bstp_idx),select_x]
                  y_oob = y_num[!(1:length(train_idx) %in% bstp_idx),, drop = F]
                  delta_obb = delta_train[!(1:length(train_idx) %in% bstp_idx),, drop = F]
                  
                  # x_oob = x_num[bstp_idx,select_x, drop = F]
                  # y_oob = y_num[bstp_idx,, drop = F]
                  # delta_obb = delta_train[bstp_idx,, drop = F]
                  
                  # EM algorithm to compute \pi_{ijkl} ####
                  n_mat = table(data.frame(cbind(x_b, y_b)), useNA = "always")
                  n_mat = eval(parse(text = paste("n_mat[", paste(rep(-k_x - 1, p_star), collapse = ","), 
                                                  paste(rep(",", q), collapse = ""), "]")))
                  
                  p_01s = p_01s_true
                  
                  diff2 = Inf
                  n_hats = vector(mode = "list", length = nrow(delind)) # sum P(Y_mis | x, Y_obs, delta)
                  cnt = 0
                  while(T){
                    cnt = cnt + 1
                    for(k in 1:nrow(delind)){
                      p_01sTx = paste("p_01s[[", k, "]]")
                      MARGINTx = paste(c(1:p_star, ((p_star + 1):(p_star + q))[delind[k,,drop = F] == 1]), collapse = ", ")
                      MARGINTx = paste("MARGIN = c(", MARGINTx, ")")
                      n_matTx1 = paste(rep(",", p_star), collapse = "")
                      n_matTx2 = paste(ifelse(delind[k,] == 1, -3, 3), collapse = ",")
                      n_matTx = paste("n_mat[", n_matTx1, n_matTx2, "]")
                      n_hatsTx = paste(p_01sTx, MARGINTx, n_matTx,"\"*\"", sep = ", ")
                      n_hatsTx = paste("sweep(", n_hatsTx, ")")
                      # print(n_hatsTx)
                      # print(eval(parse(text = n_hatsTx)))
                      tmp <- eval(parse(text = n_hatsTx))
                      n_hats[[k]] <- tmp
                    }
                    
                    # if(cnt == 1) n_hats_true = n_hats
                    
                    # sum(sapply(n_hats, sum))

                    n_hat = Reduce(`+`, n_hats) # \hat N(X, Y)
                    
                    Py_x = sweep(n_hat, MARGIN = 1:p_star, apply(n_hat, 1:p_star, sum), "/") # P(Y | x)
                    
                    # if(cnt == 1) Py_x_true = sweep(n_hat, MARGIN = 1:p_star, apply(n_hat, 1:p_star, sum), "/") # P(Y | x)
                    
                    Pdel1_xy = vector(mode = "list", length = q) # P(delta_j = 1 | x, y_j)
                    Pdel0_xy = vector(mode = "list", length = q) # P(delta_j = 0 | x, y_j)

                    for(k in 1:q){
                      tmp = c(1:p_star, k + p_star)
                      # print(delind[delind[,k] == 1,])
                      np1_hat = Reduce(`+`, n_hats[delind[,k,drop = F] == 1])
                      np0_hat = Reduce(`+`, n_hats[delind[,k,drop = F] == 0])
                      
                      Pdel1_xy[[k]] <- apply(np1_hat, tmp, sum) / apply(n_hat, tmp, sum)
                      Pdel0_xy[[k]] <- apply(np0_hat, tmp, sum) / apply(n_hat, tmp, sum)
                    }
                    
                    Pdel_xy = list(Pdel0_xy, Pdel1_xy)
                    
                    # if(cnt == 1) Pdel_xy_true = list(Pdel0_xy, Pdel1_xy)
                    
                    # P(y_mis | x, y_obs, delta_mis = 0, delta_obs = 1) \propto 
                    # P(del1 | x, y1) * \cdots * P(delq | x, yq) * p(Y | x)
                    p_01s_new = vector(mode = "list", length = nrow(delind))
                    for(k in 1:nrow(delind)){
                      ydex = delind[k,,drop = F]
                      tmp <- Py_x
                      # print(ydex)
                      for(l in 1:q){
                        MARGIN = c(1:p_star, l + p_star)
                        tmp <- sweep(tmp, MARGIN = MARGIN, Pdel_xy[[ydex[,l] + 1]][[l]], "*")
                      }
                      MARGIN2 = c(1:p_star, ((p_star + 1):(p_star + q))[delind[k,,drop = F] == 1])
                      tmp = sweep(tmp, MARGIN = MARGIN2, apply(tmp, MARGIN2, sum), "/")
                      
                      # If # of (x & y = NA) = 0, we ignore the associated probs.
                      # if(any(is.na(tmp))) print(sum(is.na(tmp)))
                      tmp[is.na(tmp)] <- 1/(2^sum(delind[k,,drop = F] == 0))
                      p_01s_new[[k]] <- tmp
                    }
                    
                    p_mat = n_hat / nrow(x_b)
                    
                    # if(cnt %% 100 == 0){ ####
                    #   lik_seg_all = 0
                    #   for(k in 1:nrow(delind)){
                    #     ydex = delind[k,,drop = F]
                    #     z_b_sub = z_b[apply(delta_b, 1, function(id) all(id == ydex)),]
                    #     if(nrow(z_b_sub) != 0){
                    #       tmpftn = function(x){
                    #         if(!is.na(x)) x
                    #         else c(1,2)
                    #       }
                    #       lik_seg = sum(apply(z_b_sub, 1, function(rows){
                    #         if(k != 1) cands = expand.grid(sapply(rows, tmpftn))
                    #         else cands = t(rows)
                    #         # print(cands)
                    #         log(sum(apply(cands, 1, function(k2) {
                    #           p_mat[t(k2)] * prod(sapply(1:q, function(l){
                    #             texts = paste(k2[c(1:p_star, p_star + l)], collapse = ",")
                    #             texts = paste("Pdel_xy[[ydex[,", l, "] + 1]][[", l, "]][", texts, "]")
                    #             # print(texts)
                    #             eval(parse(text = texts))
                    #           }  ))
                    #           # if(p_mat[t(k2)] == 0){
                    #           #   print(t(k2))
                    #           #   stop()
                    #           # } 
                    #         })))
                    #       }))
                    #       # cands = expand.grid(sapply(z_oob_sub[1,], tmpftn))
                    #       # log(sum(apply(cands, 1, function(k) {
                    #       #   p_mat[t(k)] * prod(sapply(1:q, function(l){
                    #       #     texts = paste(k[c(1:p_star, p_star + l)], collapse = ",")
                    #       #     texts = paste("Pdel_xy[[ydex[,", l, "] + 1]][[", l, "]][", texts, "]")
                    #       #     # print(texts)
                    #       #     eval(parse(text = texts))
                    #       #   }  ))
                    #       # })))
                    #       
                    #     }else{
                    #       lik_seg = 0
                    #     }
                    #     # print(lik_seg)
                    #     lik_seg_all = lik_seg_all + lik_seg
                    #   }
                    #   print(lik_seg_all)
                    # }
                    #####
                    
                    p_tmp = sapply(p_01s, cbind)
                    p_tmp_new = sapply(p_01s_new, cbind)

                    diff = norm(p_tmp_new - p_tmp, "2")
                    # print(diff)
                    if(diff < 10^(-4) | cnt > 1000){
		                # if(cnt > 10){
                    #  if(diff < 10^(-1)){
                    # if(diff2 < 10^(-1)){
                      # if(T){
                      #  if(diff < 10^(-1)){
                      break
                    } 
                    else{
                      p_01s = p_01s_new
                      diff2 = diff
                    }
                  }
                  
                  # print("EM is done!")
                  
                  p_01s = p_01s_new
                  
                  # p_tmp_true = sapply(p_01s_true, cbind) ####
                  # 
                  # hist(p_tmp_true[,2])
                  # hist(p_tmp_new[,2])
                  # plot(p_tmp_true[,2], p_tmp_new[,2])
                  # 
                  # hist(p_tmp_true[,3])
                  # hist(p_tmp_new[,3])
                  # plot(p_tmp_true[,3], p_tmp_new[,3])
                  # 
                  # hist(p_tmp_true[,4])
                  # hist(p_tmp_new[,4])
                  # plot(p_tmp_true[,4], p_tmp_new[,4])
                  # 
                  # plot(cbind(Py_x_true, Py_x))
                  # 
                  # plot(Pdel_xy[[2]][[1]], Pdel_xy_true[[2]][[1]])
                  # plot(Pdel_xy[[2]][[2]], Pdel_xy_true[[2]][[2]])
                  # 
                  # boxplot(sapply(n_hats, cbind) - sapply(n_hats_true, cbind))
                  # plot(sapply(n_hats, cbind)[,2], sapply(n_hats_true, cbind)[,2])
                  # plot(sapply(n_hats, cbind)[,3], sapply(n_hats_true, cbind)[,3])
                  # plot(sapply(n_hats, cbind)[,4], sapply(n_hats_true, cbind)[,4])
                  # 
                  # hist(p_tmp_true[,4] - p_tmp_new[,4], freq  = F)
                  # ####

                  p_mat = n_hat / nrow(x_b)
                  
                  # Compute observed likelihood l_{obs}^{(b)} for each bag b ####
                  z_oob = cbind(x_oob, y_oob + 1)
                  # table(data.frame(delta_obb))
                  
                  lik_seg_all = 0
                  for(k in 1:nrow(delind)){
                    ydex = delind[k,,drop = F]
                    z_oob_sub = z_oob[apply(delta_obb, 1, function(id) all(id == ydex)),]
                    if(nrow(z_oob_sub) != 0){
                      tmpftn = function(x){
                        if(!is.na(x)) x
                        else c(1,2)
                      }
                    lik_seg = sum(apply(z_oob_sub, 1, function(rows){
                      if(k != 1) cands = expand.grid(sapply(rows, tmpftn))
                      else cands = t(rows)
                      # print(cands)
                      log(sum(apply(cands, 1, function(k2) {
                        p_mat[t(k2)] * prod(sapply(1:q, function(l){
                          texts = paste(k2[c(1:p_star, p_star + l)], collapse = ",")
                          texts = paste("Pdel_xy[[ydex[,", l, "] + 1]][[", l, "]][", texts, "]")
                          # print(texts)
                          eval(parse(text = texts))
                        }  ))
                        # if(p_mat[t(k2)] == 0){
                        #   print(t(k2))
                        #   stop()
                        # } 
                      })))
                    }))
                      # cands = expand.grid(sapply(z_oob_sub[1,], tmpftn))
                      # log(sum(apply(cands, 1, function(k) {
                      #   p_mat[t(k)] * prod(sapply(1:q, function(l){
                      #     texts = paste(k[c(1:p_star, p_star + l)], collapse = ",")
                      #     texts = paste("Pdel_xy[[ydex[,", l, "] + 1]][[", l, "]][", texts, "]")
                      #     # print(texts)
                      #     eval(parse(text = texts))
                      #   }  ))
                      # })))
                      
                    }else{
                      lik_seg = 0
                    }
                    # print(lik_seg)
                    lik_seg_all = lik_seg_all + lik_seg
                  }

                  
                  l_obs = sum(lik_seg_all) / nrow(z_oob) # observed likelihood
                  if(is.nan(l_obs)) next() 
                  
                  expl_obs = exp(lambda * l_obs)
                  
                  # print(paste("select_x =", select_x[1], select_x[2]))
                  # print(paste("expl_obs =", expl_obs))
                  
                  # Find the conditional probabilities ####
                  # P(y_mis | y_obs, X_s^{(b)})
                  
                  phats = vector(mode = "list", length = nrow(delind))
                  for(k in 1:nrow(delind)){
                    MARGIN2 = c(1:p_star, ((p_star + 1):(p_star + q))[delind[k,,drop = F] == 1])
                    # print(MARGIN2)
                    phats[[k]] <- sweep(p_mat, MARGIN = MARGIN2, apply(p_mat, MARGIN2, sum), FUN = "/") 
                  }
                  
                  # sapply(phats, function(k) any(is.nan(k)))

                  w_B = c(w_B, expl_obs)
                  
                  # Fractional Imputation step ####
                  z = cbind(x_num[,select_x],y_num + 1)
                  z_imp = NULL
                  
                  if(p_star == 1) colnames(z)[1] <- "Var"
                  
                  lik_seg_all = 0
                  for(k in 1:nrow(delind)){
                    ydex = delind[k,,drop = F]
                    z_sub = z[apply(delta, 1, function(id) all(id == ydex)),]
                    if(nrow(z_sub) != 0){
                      tmpftn = function(x){
                        if(!is.na(x)) x
                        else c(1,2)
                      }
                      
                      seg = apply(z_sub, 1, function(rows){
                        if(k != 1) cands = expand.grid(sapply(rows, tmpftn))
                        else cands = t(rows)
                        # print(cands)
                        # cands = t(z_sub[1,])
                        p_tmp = apply(cands, 1, function(k2) {
                          p_mat[t(k2)] * prod(sapply(1:q, function(l){
                            texts = paste(k2[c(1:p_star, p_star + l)], collapse = ",")
                            texts = paste("Pdel_xy[[ydex[,", l, "] + 1]][[", l, "]][", texts, "]")
                            # print(texts)
                            eval(parse(text = texts))
                          }  ))
                        })
                        if(sum(p_tmp)== 0){
                          lenp = length(p_tmp)
                          p_tmp = rep(1 / lenp, lenp)
                        }
                        cbind(cands, w = p_tmp / sum(p_tmp))
                      })
                      if(k != 1){
                        z_imp = rbind(z_imp, do.call(rbind, seg))
                      }else{
                        row.names(seg) <- c(colnames(z), "w")
                        z_imp = rbind(z_imp, t(seg))
                      }
                    }
                  }
                  
                  # table(data.frame(z))
                  if(b == 1){
                    z_imp_res <- expl_obs * as.matrix(z_imp)
                  }else{
                    z_imp_res = z_imp_res + expl_obs * as.matrix(z_imp)
                  }
                  
                  if(any(is.nan(z_imp_res))) stop("any(is.nan(z_imp_res))")
                }
                
                if(sum(w_B) == 0){
                  warning(sum(w_B) == 0)
                  next
                }
                
                z_imp_res = z_imp_res / sum(w_B)
                
                wcol = p_star + q + 1
                
                # print(paste("sum(z_imp_res[,5]) = ", sum(z_imp_res[,wcol])))
                
                theta_prop = colSums(z_imp_res[,(p_star+1):(p_star+q) ,drop = F] * z_imp_res[,wcol]) / sum(z_imp_res[,wcol]) - 1
                
                theta_full
                
                # print(theta_prop)
                
                theta_cc = colMeans(y_num, na.rm = T)
                
                rbind(theta_full, theta_cc, theta_mice, theta_prop)
                
              }

timenow2 = Sys.time()
print(timenow2 - timenow0)

for(k in 1:q){
  res_fin = sapply(res, function(x) x[,k])
  
  BIAS = rowMeans(res_fin - theta)
  SE = apply(res_fin, 1, function(x) sqrt(var(x) * (length(x)-1)/length(x) ))
  RMSE = apply(res_fin - theta, 1, function(x) sqrt(mean(x^2)))
  
  tmp_tbl <- round(cbind(BIAS = BIAS, SE = SE, RMSE = RMSE), 4)
  print(tmp_tbl)
  # xtable::xtable(tmp_tbl, digits = 4)
}

stopCluster(cl)
