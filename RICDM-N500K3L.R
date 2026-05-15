## working dir 
setwd("C:/Users/dell/Desktop/Ising/Res/ResN500K3L")

## used packages
library(IsingSampler)
library(rstan)

## true values
set.seed(12)
N <- 500 # Number of samples
J <- 18
K <- 3 # Number of attributes
CNUM <- 50
Q <- matrix(c(
  1,0,0,
  0,1,0,
  0,0,1,
  1,0,0,
  0,1,0,
  0,0,1,
  1,0,0,
  0,1,0,
  0,0,1,
  1,1,0, 
  1,0,1,
  0,1,1,
  1,1,0,
  1,0,1,
  0,1,1,
  1,1,1,
  1,1,1,
  1,1,1
), nrow=J, ncol=K, byr=TRUE)

## generate all attribute patterns
generate_alpha_patterns <- function(K){
  alpha <- as.matrix(expand.grid(rep(list(0:1), K)))
  colnames(alpha) <- paste0("att", 1:K)
  return(alpha)
}

alpha <- generate_alpha_patterns(K)

pK <- c(0.5,0.5,0.5)
pC0 <- rep(1,2^K)

for(c in 1:2^K){
  for(k in 1:K){
    pC0[c] <- pC0[c]*(pK[k]^alpha[c,k])*((1-pK[k])^(1-alpha[c,k]))
  }
}

h_true <- runif(J,-2,-1)
BT0 <-  runif(J,3,4)
J_true <- matrix(0, J, J)
J_true[7,8] <- J_true[8,7] <- 1.0
J_true[8,9] <- J_true[9,8] <- 1.0
J_true[12,15] <- J_true[15,12] <- 1.0


# Stan model
stanm <- stan_model(model_code = "
data {
  int<lower=1> N;              // sample size
  int<lower=1> J;              // number of items
  int<lower=1> K;              // number of attributes
  int<lower=1> C;              // number of attribute profiles
  int<lower=0,upper=1> Y[N, J];     // item resposes
  int<lower=0,upper=1> Q[J, K];     // Q matrix
  int<lower=0,upper=1> alpha[C, K]; // all attribute profiles
}

parameters {
  vector<lower=0>[J] guess;  
  vector[J] h;   
  vector[(J-2*K)*(J-2*K-1)/2] J_vec; 
  vector<lower=0>[(J-2*K)*(J-2*K-1)/2] lam;
  simplex[C] pC;  
}

transformed parameters {
  matrix[J,J] Jmat; 
  
  {int idx = 1;
    for(i in 1:(J-1)){
      for(j in (i+1):J){
        if(i>(2*K)){
          Jmat[i, j] = J_vec[idx];
          Jmat[j, i] = J_vec[idx];
          idx += 1;
        }else{
          Jmat[i, j] = 0;
          Jmat[j, i] = 0;
        }        
      }
    }
    for(i in 1:J) Jmat[i, i] = 0;
  }
}

model {
  guess ~ normal(0, 2);
  pC ~ dirichlet(rep_vector(1.0, C));   
  h ~ normal(0, 2);
  J_vec ~ double_exponential(0, 1./sqrt(lam));   
  lam ~ gamma(9, 1);
  
  // Marginal likelihood: sum over all latent classes for each subject
  for(n in 1:N){
    vector[C] lp;
    for(c in 1:C){
      lp[c] = log(pC[c]);
      for(j in 1:J){
        lp[c] += bernoulli_lpmf(Y[n, j] | inv_logit(h[j] + guess[j]* prod(alpha[c, ] .^ Q[j,]) + dot_product(to_vector(Y[n,]), Jmat[, j])) );
      }
    }
    target += log_sum_exp(lp);  
  }
}

generated quantities {
  matrix[N, C] prob_class;
  matrix[N, K] prob_attr;
  for(n in 1:N){
    vector[C] log_postc;
    real logmarg = negative_infinity();
    for(c in 1:C){
      log_postc[c] = log(pC[c]);
      for(j in 1:J){
        log_postc[c] += bernoulli_lpmf(Y[n, j] | inv_logit(h[j] + guess[j]* prod(alpha[c, ] .^ Q[j,]) + dot_product(to_vector(Y[n,]), Jmat[, j])) );
      }
      logmarg = log_sum_exp(logmarg, log_postc[c]);
    }
   for(c in 1:C){
     prob_class[n,c] = exp(log_postc[c] - logmarg);
   }
   for(k in 1:K){
     prob_attr[n,k] = 0;
     for(c in 1:C){
       prob_attr[n,k] += prob_class[n,c]*alpha[c,k];
     }
   }
 }
}
", verbose = T )

## generated data
for(CIR in 1:CNUM){
  X <- rmultinom(N, size = 1, prob = pC0)
  Y <- matrix(0, nrow = N, ncol = J)
  for(n in 1:N){
    Y[n, ] <- IsingSampler(1, J_true,( h_true + BT0 * apply(matrix(rep(alpha[which(X[,n]==1),], J), nrow=J, byr=TRUE)^Q, 1, prod) ))
  }
  write(Y, file="Y.txt", ncol=dim(Y)[1], append=T)
  write(t(X), file="X.txt", ncol=dim(t(X))[1], append=T)
  print(CIR)
}

## repetations
for(CIR in 1:CNUM){
  bt <- proc.time()
  Y <- matrix(0, nrow=N, ncol=J)
  Y <- matrix(scan("Y.txt", skip=(CIR-1)*J, nlines=J), nrow=N, ncol=J)

  # model input and starting values
  # data for rstan
  fdata <- list(
    N = nrow(Y),
    J = ncol(Y),
    K = K,
    C = 2^K,
    Y = Y,
    Q = Q,
    alpha = alpha
  )

  fit_mcmc <- sampling(stanm, data = fdata, iter = 3000, chains = 3, cores = 3, control = list(adapt_delta = 0.95))  

  aa <- summary(fit_mcmc, probs = c(0.025, 0.975), pars = c("guess","h","J_vec","pC","prob_class","prob_attr") )

  write(get_num_divergent(fit_mcmc), file="diver.txt", ncol=1, append=TRUE, sep="\t")
  write(aa$summary[,c(7)], file="Rhh.txt", ncol=length(aa$summary[,c(7)]), append=TRUE, sep="\t")
  write(aa$summary[,c(1)], file="amean.txt", ncol=length(aa$summary[,c(1)]), append=TRUE, sep="\t")
  write(aa$summary[,c(3)], file="asd.txt", ncol=length(aa$summary[,c(3)]), append=TRUE, sep="\t")
  write(aa$summary[,c(4)], file="ainl.txt", ncol=length(aa$summary[,c(4)]), append=TRUE, sep="\t")
  write(aa$summary[,c(5)], file="ainr.txt", ncol=length(aa$summary[,c(5)]), append=TRUE, sep="\t")

  print(CIR)
  et<-proc.time()
  print((et-bt)[3])
  write((et-bt)[3], file="time.txt", ncol=1, append=TRUE, sep="\t")
}

