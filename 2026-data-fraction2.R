##clear screen
rm(list=ls())

## used packages
library(rstan)
library(CDM)

set.seed(12)
Yt <- data.fraction2$data
Qt <- data.fraction2$q.matrix3

Y <- cbind(Yt[,4],Yt[,1],Yt[,2],Yt[,9],Yt[,5],Yt[,3],Yt[,c(6,7,8,10,11)])
Q <- rbind(Qt[4,],Qt[1,],Qt[2,],Qt[9,],Qt[5,],Qt[3,],Qt[c(6,7,8,10,11),])


N <- dim(Y)[1] # Number of samples
J <- dim(Y)[2]
K <- dim(Q)[2] # Number of attributes

# generate all attribute patterns
generate_alpha_patterns <- function(K) {
  alpha <- as.matrix(expand.grid(rep(list(0:1), K)))
  colnames(alpha) <- paste0("att", 1:K)
  return(alpha)
}

alpha <- generate_alpha_patterns(K)


## DINA model
fit <- gdina(data = Y, q.matrix = Q, rule = "DINA", linkfct = "logit")

coef(fit)$est
fit$attribute.patt

fit$ic

Ypa <- matrix(0,nrow=J,ncol=J)

Je1 <- matrix(coef(fit)$est,ncol=2,byr=TRUE)
for(i in 1:(J-1)){
  for(j in (i+1):J){
    itep <- c(i,j)
    pp00 <- 0
    pp01 <- 0
    pp10 <- 0
    pp11 <- 0
    for(c in 1:2^K){
      p11 <- exp(Je1[itep[1],1] + Je1[itep[1],2]*prod(alpha[c,]^Q[itep[1],]))/(1+exp(Je1[itep[1],1] + Je1[itep[1],2]*prod(alpha[c,]^Q[itep[1],])))
      p21 <- exp(Je1[itep[2],1] + Je1[itep[2],2]*prod(alpha[c,]^Q[itep[2],]))/(1+exp(Je1[itep[2],1] + Je1[itep[2],2]*prod(alpha[c,]^Q[itep[2],])))
  
      pp00 <- pp00 + fit$attribute.patt[c,1]*(1-p11)*(1-p21)
      pp01 <- pp01 + fit$attribute.patt[c,1]*(1-p11)*p21
      pp10 <- pp10 + fit$attribute.patt[c,1]*p11*(1-p21)
      pp11 <- pp11 + fit$attribute.patt[c,1]*p11*p21
    }

    E00 <- N*pp00
    E01 <- N*pp01
    E10 <- N*pp10
    E11 <- N*pp11

  Ypa[i,j] <- (table(Y[,itep])[1,1]-E00)^2/E00 + (table(Y[,itep])[1,2]-E01)^2/E01 + (table(Y[,itep])[2,1]-E10)^2/E10 + (table(Y[,itep])[2,2]-E11)^2/E11
  }
}


personDINA <- do.call(rbind, lapply(strsplit(fit$pattern$map.est, ""), as.integer))


##RI-CDM
# Stan model
stanm2 <- stan_model(model_code = "
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


fdata <- list(
  N = nrow(Y),
  J = ncol(Y),
  K = K,
  C = 2^K,
  Y = Y,
  Q = Q,
  alpha = alpha
)

fit_mcmc2 <- sampling(stanm2, data = fdata, iter = 3000, chains = 3, cores = 3, control = list(adapt_delta = 0.95))  

aa2 <- summary(fit_mcmc2, probs = c(0.025, 0.975), pars = c("pC","h","guess","J_vec","prob_class","prob_attr") )

get_num_divergent(fit_mcmc2)


person <- matrix(aa2$summary[(2^K+2*J+(J-2*K)*(J-2*K-1)/2+N*(2^K)+1):(2^K+2*J+(J-2*K)*(J-2*K-1)/2+N*(2^K)+N*K),1],ncol=K,byr=TRUE)

person1 <- ifelse(person>0.5, 1, 0)

sum( rowSums(personDINA == person1)==K )/N


YpaIs <- matrix(0,nrow=J,ncol=J)

Ypattern <- as.matrix(expand.grid(rep(list(0:1), J)))
colnames(Ypattern) <- paste0("Y", 1:J)

pCe <- aa2$summary[1:2^K,1]
he <- aa2$summary[(2^K+1):(2^K+J),1]
ge <- aa2$summary[(2^K+J+1):(2^K+2*J),1]

J_vece <- aa2$summary[(2^K+2*J+1):(2^K+2*J+(J-2*K)*(J-2*K-1)/2),1]
Jmate <- matrix(0, nrow=J, ncol=J)
ind <- 1
for(i in (2*K+1):(J-1)){
  for(j in (i+1):J){
    Jmate[i,j] <- J_vece[ind]
    Jmate[j,i] <- J_vece[ind]
    ind <- ind +1
  }
}

for(i in 1:(J-1)){
  for(j in (i+1):J){
    itep <- c(i,j)

    p00 <- 0
    p01 <- 0
    p10 <- 0
    p11 <- 0

    for(c in 1:2^K){
      hec <-  he
      for(jj in 1:J){
        hec[jj] <-  hec[jj] + ge[jj]*prod(alpha[c, ]^Q[jj,])
      }
      energy <- numeric(nrow(Ypattern))
      for(k in 1:nrow(Ypattern)){
        y <- as.numeric(Ypattern[k,])
        e1 <- sum(hec*y)
        e2 <- 0
        for(ii in 1:(J-1)){
          for(jj in (ii+1):J){
            e2 <- e2 + Jmate[ii,jj]*y[ii]*y[jj]
          }
        }
      energy[k] <- e1 + e2
    }
    weight <- exp(energy)
    Zw <- sum(weight)

    p00 <- p00 + pCe[c]*(sum(weight[(Ypattern[,paste0("Y",itep)[1]]==0)&(Ypattern[,paste0("Y",itep)[2]]==0)])/Zw)
    p01 <- p01 + pCe[c]*(sum(weight[(Ypattern[,paste0("Y",itep)[1]]==0)&(Ypattern[,paste0("Y",itep)[2]]==1)])/Zw)
    p10 <- p10 + pCe[c]*(sum(weight[(Ypattern[,paste0("Y",itep)[1]]==1)&(Ypattern[,paste0("Y",itep)[2]]==0)])/Zw)
    p11 <- p11 + pCe[c]*(sum(weight[(Ypattern[,paste0("Y",itep)[1]]==1)&(Ypattern[,paste0("Y",itep)[2]]==1)])/Zw)
   }

   E00 <- N*p00
   E01 <- N*p01
   E10 <- N*p10
   E11 <- N*p11

   YpaIs[i,j] <- (table(Y[,itep])[1,1]-E00)^2/E00 + (table(Y[,itep])[1,2]-E01)^2/E01 + 
                (table(Y[,itep])[2,1]-E10)^2/E10 + (table(Y[,itep])[2,2]-E11)^2/E11

  }
}


lik <- rep(0, N)
for(n in 1:N){
  for(c in 1:2^K){
      hec <-  he
      for(jj in 1:J){
        hec[jj] <-  hec[jj] + ge[jj]*prod(alpha[c, ]^Q[jj,])
      }
      energy <- numeric(nrow(Ypattern))
      for(k in 1:nrow(Ypattern)){
        y <- as.numeric(Ypattern[k,])
        e1 <- sum(hec*y)
        e2 <- 0
        for(ii in 1:(J-1)){
          for(jj in (ii+1):J){
            e2 <- e2 + Jmate[ii,jj]*y[ii]*y[jj]
          }
        }
      energy[k] <- e1 + e2
    }
    weight <- exp(energy)
    Zw <- sum(weight)

  lik[n] <- lik[n] + pCe[c]*(weight[which(rowSums(Ypattern==matrix(rep(Y[n, ],each=2^J),nrow=2^J))==J)]/Zw)
  }
}

##AIC 
sum(-2*log(lik))+2*(2^K-1+2*J+(J-2*K)*(J-2*K-1)/2) 

##BIC 
sum(-2*log(lik))+(2^K-1+2*J+(J-2*K)*(J-2*K-1)/2)*log(N)




