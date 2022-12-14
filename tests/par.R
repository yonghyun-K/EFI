# Simulation  setup ####
n = 100
print(paste("n(Sample size) =", n))

n_B = n
print(paste("n_B(Bootstrap sample size) =", n_B))

p = 5
print(paste("p =", p))

p_star = 2
print(paste("p_star =", p_star))

q = 2
print(paste("q =", q))

B = choose(p, p_star)
# B = 5
print(paste("B(The number of bootstraps) =", B))

SIMNUM = round(min(c(detectCores() / 3 * 2, 100)))
# SIMNUM = 10
print(paste("SIMNUM(MC size) =", SIMNUM))

k_x = 2
print(paste("k_x(Number of categories of X) =", k_x))

# lambda = 120
# print(paste("lambda =", lambda))

maxit = 500
print(paste("maxit =", maxit))

lambda_vec_text = "seq(from = 0, to = 300, len = SIMNUM)"
lambda_vec = eval(parse(text = lambda_vec_text))
print(paste("lambda_vec_text =", lambda_vec_text))

cores=SIMNUM
print(paste("cores =", cores))
cl <- makeSOCKcluster(cores)
registerDoSNOW(cl)

progress <- function(n) cat(sprintf("task %d is complete\n", n))
opts <- list(progress=progress)

set.seed(1)

# p_x = runif(k_x, min = 0.25, max = 0.75)
# p_x = p_x / sum(p_x)
# print(paste("p_x =", paste(round(p_x,5), collapse = " ")))
p_x = sapply(1:p, function(x) runif(k_x, min = 0.25, max = 0.75))
p_x = sweep(p_x, MARGIN = 2, apply(p_x, 2, sum), "/")
print("p_x = "); print(p_x)

# p_Y = runif(k_x, min = 0.25, max = 0.75)
# print(paste("p_Y[1:k_x] =", paste(round(p_Y[1:k_x],5), collapse = " ")))
# p_Y = sapply(1:q, function(x) runif(k_x, min = 0.25, max = 0.75))
p_Y = sapply(1:q, function(x) {
  while(T){
    tmp = runif(k_x, min = 0.2, max = 0.8) 
    if(max(tmp) - min(tmp) > 0.4) break
  }
  return(tmp)
})
# p_Y = sapply(1:q, function(x) 1: k_x / (k_x + 1))
print("p_Y = "); print(p_Y)

# theta = sum(p_x * p_Y[1:k_x])
# print(paste("theta =", round(theta, 5)))
theta = colSums(p_x[,1:q] * p_Y)
# theta = colSums(rowMeans(p_x[,1:3]) * p_Y)
print("theta = "); print(theta)

# misstype = "NMAR"
misstype = "MAR"
# misstype = "SCens"
print(paste("misstype =", misstype))

# p_delta_ftn = function(k) rep(0.5, n); print("MCAR") # 0.5

# p_delta_ftn = function(k) X_num[,k] / 5 + 0.2; print("MAR")# 0.5
 # p_delta_ftn = function(k) X_num[,k] / 3 + 1/3; print("MAR") # 0.84

# p_delta_ftn = function(k) 1 / (1 + exp(-(X_num[,k] + X_num[,k+1] - X_num[,k+2]))); print("MAR") # 0.8
# p_delta_ftn = function(k) 1 / (1 + exp(-(-1 + X_num[,k] + X_num[,k+1] - X_num[,k+2]))); print("MAR") # 0.61
# p_delta_ftn = function(k) 1 / (1 + exp(-(-X_num[,k] / 2 + X_num[,k + 3] / 4 + X_num[,k + 4] / 4))); print("MAR") # 0.415
# p_delta_ftn = function(k) 1 / (1 + exp(-(-X_num[,k] / 2 + X_num[,k + 2] / 2))); print("MAR") # 0.415
p_delta_ftn = function(k) 1 / (1 + exp(-(X_num[,k] * 2 - 3))); print("MAR") # 0.415
# p_delta_ftn = function(k) 1 / (1 + exp(-(X_num[,k+1] / 2 - 1 / 4))); print("MAR") # 0.415

# p_delta_ftn = function(k) 1 / (1 + exp(-(1 - X_num[,k] / 4 + 1 / 4 - Y_num[,k] / 4))); print("NMAR") # 0.415
# p_delta_ftn = function(k) 1 / (1 + exp( -(1.5 + X_num[,k] + X_num[,k+1] - X_num[,k+2] - 5 * Y_num[,k]))); print("NMAR") # 0.561
# p_delta_ftn = function(k) 1 / (1 + exp( -(X_num[,k] + X_num[,k+1] - 5 * Y_num[,k]))); print("NMAR")
# p_delta_ftn = function(k) 1 / (1 + exp( -( X_num[,k] / 4 - Y_num[,k] / 4))); print("NMAR")
# p_delta_ftn = function(k) 1 / (1 + exp( -(Y_num[,k] / 4))); print("NMAR")

print("p_delta_ftn")
print(p_delta_ftn)
