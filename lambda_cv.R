n = 500
p = 2
n_B = 300
B = 10
SIMNUM = 50
res = NULL
K = 10

lambda_vec = seq(from = 0.0, to = 20, len = 100)

set.seed(1)



p_x = c(0.3, 0.25, 0.25, 0.2)
X = matrix(rep(1:4, n * p)[rmultinom(n * p, 1, p_x) == 1], nr = n, nc = p)
colnames(X) <- paste("X", 1:p, sep = "")
# table(X) / n

p_Y1 = c(0.4, 0.3, 0.4, 0.8)
p1 = p_Y1[X[, 1]]
Y1 = rbinom(n, 1, p1)
theta1 = sum(p_x * p_Y1)

p2 = c((0.8 + Y1) / 2, 0.8 * Y1 + 0.1, 0.9 - 0.5 * Y1, 0.6)[X[, 1]]
Y2 = rbinom(n, 1, p2)

p_Y2 = c((0.8 + theta1) / 2, 0.8 * theta1 + 0.1, 0.9 - 0.5 * theta1, 0.6)
theta2 = sum(p_x * p_Y2)

theta = c(theta1, theta2)

delta1 = rbinom(n, 1, 0.8)
delta2 = rbinom(n, 1, 0.7)

Y_ogn = cbind(Y1, Y2)
theta_full = colMeans(Y_ogn)
delta = cbind(delta1, delta2)

Y = Y_ogn
Y[delta == 0] = NA

# lambda = 0.1




bstp_idx_B = NULL
for (b in 1:B) {
  bstp_idx = sample(1:((1 - 1 / K) * n), n_B, replace = TRUE)
  bstp_idx_B = cbind(bstp_idx_B, bstp_idx)
}


for (lambda in lambda_vec) {
  # set.seed(5)
  print(lambda)
  
  res_K = rep(0, K)
  library(caret)
  cuts = cut(sample(1:n, n),breaks=10,labels=FALSE)
  for(k in 1:K){
    # flds = createFolds(1:n, k = K, list = TRUE)
    # test_idx = flds[[k]]
    # train_idx = do.call(c, flds[-k])


    test_idx = which(cuts == k)
    train_idx = which(cuts != k)
    
    names(train_idx) = NULL
    
    x = X[-test_idx, ]
    y = Y[-test_idx, ]
    
    x_test = X[test_idx, ]
    y_test = Y[test_idx, ]
    
    w_B = NULL
    
    w_01 = array(0, dim = c(4, 4, 2, 2))
    w_10 = array(0, dim = c(4, 4, 2, 2))
    w_00 = array(0, dim = c(4, 4, 2, 2))
    
    p_mat_lambda = array(0, dim = c(4, 4, 2, 2))
    
    for (b in 1:B) {
      bstp_idx = bstp_idx_B[,b]
      x_b = x[bstp_idx, ]
      y_b = y[bstp_idx, ]
      
      x_oob = x[!(1:length(train_idx) %in% bstp_idx), ]
      y_oob = y[!(1:length(train_idx) %in% bstp_idx), ]
      
      # unique(cbind(x_b, y_b))
      # unique(cbind(x_b, y_b))[complete.cases(unique(cbind(x_b, y_b))),]
      
      n_mat = table(data.frame(cbind(x_b, y_b)), useNA = "ifany")
      # apply(n_mat, (p+1) : (p+2), sum)
      # apply(n_mat, 1:p, sum)
      
      # p_mat = n_mat[,,-3,-3] / sum(n_mat[,,-3,-3])
      
      # (n_mat[,,-3,-3] / sum(n_mat[,,-3,-3]))[,,2,2]
      # n_mat[,,3,3]
      
      # n_mat[,,3,-3, drop = F] # n^(1101)
      # p_mat[,,1,] * n_mat[,,3,-3] / apply(p_mat, c(1,2,4), sum)
      
      p_mat = array(1, dim = dim(n_mat[, , -3, -3])) / length(n_mat[, , -3, -3])
      
      while (T) {
        nhat_11 = n_mat[, , -3, -3] # n^(1111)
        nhat_01 = sweep(p_mat,
                        MARGIN = c(1, 2, 4),
                        n_mat[, , 3, -3] / apply(p_mat, c(1, 2, 4), sum),
                        FUN = "*") #n_hat^(1101)
        nhat_10 = sweep(p_mat,
                        MARGIN = c(1, 2, 3),
                        n_mat[, , -3, 3] / apply(p_mat, c(1, 2, 3), sum),
                        FUN = "*") #n_hat^(1110)
        nhat_00 = sweep(p_mat,
                        MARGIN = c(1, 2),
                        n_mat[, , 3, 3] / apply(p_mat, c(1, 2), sum),
                        FUN = "*") #n_hat^(1100)
        
        nhat_01[is.nan(nhat_01)] <- 0
        nhat_10[is.nan(nhat_10)] <- 0
        nhat_00[is.nan(nhat_00)] <- 0
        
        p_mat2 = (nhat_11 + nhat_01 + nhat_10 + nhat_00) / n_B
        
        # print(norm(p_mat2 - p_mat, "2"))
        if (norm(p_mat2 - p_mat, "2") < 10 ^ (-3))
          break
        else
          p_mat = p_mat2
      }
      
      # plot(c(table(data.frame(cbind(X[train_idx,][bstp_idx,], Y_ogn[train_idx,][bstp_idx,]))) / n_B),
      #      c(p_mat)); abline(0,1)
      
      # plot(c(n_mat[,,-3,-3] / sum(n_mat[,,-3,-3])), c(p_mat2)); abline(0,1)
      # plot(c(table(data.frame(cbind(X[train_idx,][bstp_idx,], Y_ogn[train_idx,][bstp_idx,]))) / n_B),
      #      c(n_mat[,,-3,-3] / sum(n_mat[,,-3,-3])))
      
      # unique(cbind(x_oob, y_oob))
      
      z_oob = cbind(x_oob, y_oob)
      z_oob_11 = z_oob[!is.na(z_oob[, 3]) & !is.na(z_oob[, 4]),]
      z_oob_01 = z_oob[is.na(z_oob[, 3]) & !is.na(z_oob[, 4]),]
      z_oob_10 = z_oob[!is.na(z_oob[, 3]) & is.na(z_oob[, 4]),]
      z_oob_00 = z_oob[is.na(z_oob[, 3]) & is.na(z_oob[, 4]),]
      
      l_11 = apply(z_oob_11, 1, function(x)
        log(p_mat[t(as.character(x))]))
      l_01 = apply(z_oob_01, 1, function(x) {
        tmp = as.character(x)
        tmp1 = tmp; tmp1[3] <- "0"
        tmp2 = tmp; tmp2[3] <- "1"
        log(p_mat[t(tmp1)] + p_mat[t(tmp2)])
      })
      l_10 = apply(z_oob_10, 1, function(x) {
        tmp = as.character(x)
        tmp1 = tmp; tmp1[4] <- "0"
        tmp2 = tmp; tmp2[4] <- "1"
        log(p_mat[t(tmp1)] + p_mat[t(tmp2)])
      })
      l_00 = apply(z_oob_00, 1, function(x) {
        tmp = as.character(x)
        tmp1 = tmp; tmp1[3:4] <- c("0", "0")
        tmp2 = tmp; tmp2[3:4] <- c("0", "1")
        tmp3 = tmp; tmp3[3:4] <- c("1", "0")
        tmp4 = tmp; tmp4[3:4] <- c("1", "1")
        log(p_mat[t(tmp1)] + p_mat[t(tmp2)] + p_mat[t(tmp3)]+ p_mat[t(tmp4)])
      })
      
      # if (nrow(z_test_01) != 0) {
      #   l_01 = apply(z_test_01, 1, function(x) {
      #     tmp = as.character(x)
      #     tmp1 = tmp
      #     tmp1[3] <- "0"
      #     tmp2 = tmp
      #     tmp2[3] <- "1"
      #     log(p_mat[t(tmp1)] + p_mat[t(tmp2)])
      #   })
      # } else{
      #   l_01 = 0
      # }
      # if (nrow(z_test_10) != 0) {
      #   l_10 = apply(z_test_10, 1, function(x) {
      #     tmp = as.character(x)
      #     tmp1 = tmp
      #     tmp1[4] <- "0"
      #     tmp2 = tmp
      #     tmp2[4] <- "1"
      #     log(p_mat[t(tmp1)] + p_mat[t(tmp2)])
      #   })
      # } else{
      #   l_10 = 0
      # }
      # if (nrow(z_test_00) != 0) {
      #   l_00 = apply(z_test_00, 1, function(x) {
      #     tmp = as.character(x)
      #     tmp1 = tmp
      #     tmp1[3:4] <- c("0", "0")
      #     tmp2 = tmp
      #     tmp2[3:4] <- c("0", "1")
      #     tmp3 = tmp
      #     tmp3[3:4] <- c("1", "0")
      #     tmp4 = tmp
      #     tmp4[3:4] <- c("1", "1")
      #     log(p_mat[t(tmp1)] + p_mat[t(tmp2)] + p_mat[t(tmp3)] + p_mat[t(tmp4)])
      #   })
      # } else{
      #   l_00 = 0
      # }
      
      l_obs = sum(c(l_11, l_01, l_10, l_00)) / nrow(z_oob)
      l_obs
      
      expl_obs = exp(lambda * l_obs)
      
      phat_01 = sweep(p_mat,
                      MARGIN = c(1, 2, 4),
                      apply(p_mat, c(1, 2, 4), sum),
                      FUN = "/") # Y1: mis, Y2: obs
      phat_10 = sweep(p_mat,
                      MARGIN = c(1, 2, 3),
                      apply(p_mat, c(1, 2, 3), sum),
                      FUN = "/") # Y1: obs, Y2: mis
      phat_00 = sweep(p_mat,
                      MARGIN = c(1, 2),
                      apply(p_mat, c(1, 2), sum),
                      FUN = "/") # Y1: mis, Y2: mis
      
      if (any(is.nan(phat_01))) {
        # print("is.nan(phat_01)")
        # print(expl_obs)
        next
      } else if (any(is.nan(phat_10))) {
        # print("is.nan(phat_10)")
        # print(expl_obs)
        next
      } else if (any(is.nan(phat_00))) {
        # print("is.nan(phat_00)")
        # print(expl_obs)
        next
      }
      
      bstp_idx_B = cbind(bstp_idx_B, bstp_idx)
      
      w_B = c(w_B, expl_obs)
      
      w_01 = w_01 + expl_obs * phat_01
      w_10 = w_10 + expl_obs * phat_10
      w_00 = w_00 + expl_obs * phat_00
      
      p_mat_lambda = p_mat_lambda + expl_obs * p_mat
    }
    
    w_01 = w_01 / sum(w_B)
    w_10 = w_10 / sum(w_B)
    w_00 = w_00 / sum(w_B)
    
    p_mat_lambda = p_mat_lambda / sum(w_B)
    
    py = apply(p_mat_lambda, c(3, 4), sum) # P(y = a)
    
    tmp_ftn = function(y, py) {
      if (!is.na(y[1]) & !is.na(y[2])) {
        res = log(py[y[1] + 1, y[2] + 1])
      } else if (is.na(y[1]) & !is.na(y[2])) {
        res = log(sum(py[, y[2] + 1]))
      } else if (!is.na(y[1]) & is.na(y[2])) {
        res = log(sum(py[y[1] + 1, ]))
      } else if (is.na(y[1]) & is.na(y[2])) {
        res = 0
      }
      return(res)
    }
    
    l_lambda = sum(apply(y_test, 1, tmp_ftn, py = py))
    
    res_K[k] = l_lambda
  }
  
  res = c(res, mean(res_K))
}

plot(lambda_vec, res)