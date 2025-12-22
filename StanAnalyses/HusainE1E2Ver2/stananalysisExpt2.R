library(rstan)
library(parallel)

CP1<-read.table("expt2critdata.txt",header=T)

e2data <- list(mu_prior=c(0,0,0,0),
               subj=sort(as.integer(factor(CP1$subj))),
               item=sort(as.integer(factor(CP1$item))),
               lrt = log(CP1$rt),
               distance = CP1$dist,
               expectation = CP1$exp,
               interaction = CP1$int,
               N = nrow(CP1),
               I = length(unique(CP1$subj)),
               K = length(unique(CP1$item))
)  


expt2_code <-'
data {
    int<lower=0> N;
    real lrt[N];                     //outcome
real distance[N];                     //predictor
real expectation[N];                     //predictor
real interaction[N];                     //predictor
int<lower=1> I;                 //number of subjects
int<lower=1> K;                 //number of items
int<lower=1, upper=I> subj[N];    //subject id
int<lower=1, upper=K> item[N];    //item id
vector[4] mu_prior;             //vector of zeros passed in from R
}
transformed data {
real ZERO;                      // like #define ZERO 0 in C/C++
ZERO <- 0.0;
}
parameters {
vector[4] beta;                 // intercept and slope
vector[4] u[I];                 // random intercept and slopes subj
vector[4] w[K];
real<lower=0> sigma_e;          // residual sd
vector<lower=0>[4] sigma_u;     // subj sd
vector<lower=0>[4] sigma_w;     // item sd
corr_matrix[4] Omega_u;           // correlation matrix for random intercepts and slopes subj
corr_matrix[4] Omega_w;           // correlation matrix for random intercepts and slopes item
}
transformed parameters {
matrix[4,4] D_u;
matrix[4,4] D_w;
D_u <- diag_matrix(sigma_u);
D_w <- diag_matrix(sigma_w);
}
model {
matrix[4,4] L_u;
matrix[4,4] DL_u;
matrix[4,4] L_w;
matrix[4,4] DL_w;
real mu[N]; // mu for likelihood
//priors:
beta ~ normal(0,10);
sigma_e ~ normal(0,10);
sigma_u ~ normal(0,10);
sigma_w ~ normal(0,10);
Omega_u ~ lkj_corr(4.0);
Omega_w ~ lkj_corr(4.0);
L_u <- cholesky_decompose(Omega_u);
L_w <- cholesky_decompose(Omega_w);
for (m in 1:4) {
for (n in 1:m) {
DL_u[m,n] <- L_u[m,n] * sigma_u[m];
}
}
for (m in 1:4){
for (n in (m+1):4){
DL_u[m,n] <- ZERO;
}
}
for (m in 1:4){
for (n in 1:m){
DL_w[m,n] <- L_w[m,n] * sigma_w[m];
}
}
for (m in 1:4){
for (n in (m+1):4){
DL_w[m,n] <- ZERO;
}
}
for (i in 1:I)                  // loop for subj random effects
u[i] ~ multi_normal_cholesky(mu_prior, DL_u);
for (k in 1:K)                  // loop for item random effects
w[k] ~ multi_normal_cholesky(mu_prior, DL_w);    
for (n in 1:N) {
mu[n] <- beta[1] + beta[2]*distance[n] + beta[3]*expectation[n] + beta[4]*interaction[n] 
+ u[subj[n], 1] + u[subj[n], 2]*distance[n] + u[subj[n], 3]*expectation[n] + u[subj[n], 4]*interaction[n]+ w[item[n], 1] + w[item[n], 2]*distance[n] + w[item[n], 3]*expectation[n] + w[item[n], 4]*interaction[n];
}
lrt ~ normal(mu,sigma_e);        // likelihood
}
generated quantities {
cov_matrix[4] Sigma_u;
cov_matrix[4] Sigma_w;
Sigma_u <- D_u * Omega_u * D_u;
Sigma_w <- D_w * Omega_w * D_w;
}
'

sink("expt2resultsrun2.txt")

e2.sm <- stan_model("expt2subjitem.stan", model_name = "e2subjitem")
sflist <- mclapply(1:4, mc.cores = detectCores(),
                   function(i) sampling(e2.sm, data = e2data,chains = 1, chain_id = i, seed = 12345))
e2.sf <- sflist2stanfit(sflist)
print(e2.sf,digits=4)
sink()

