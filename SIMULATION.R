rm = (list = ls())
nrep <- 100
n_1 = n_2 = 200
miu_1 = 0
p = 200        ## could switch to 200, 400, 600, 800.
miu_2 = matrix(0,p,1)
for (s in 1:10){
  miu_2[s]=1
}

## model1
segma_1 = rbind()
omega_1 = cbind()
for (s in 1:p){
  c = cbind()
  c_1 = cbind()
  for(i in 1:p){
    if(i == s){
      c[i] = 1
    }
    else{
      c[i] = 0.5
    }
    c_1[i] = c[i]^-1
  }
  segma_1 = rbind(segma_1, c)
  omega_1 = rbind(omega_1, c_1)
}

## model 2
segma_2 = rbind()
for (s in 1:p){
  c_2 = cbind()
  for(i in 1:p){
    c_2[i] = 0.8 ** abs(i - s)
  }
  segma_2 = rbind(segma_2, c_2)
  }
omega_2 = segma_2^-1

T1<-T2<-rep(NA, nrep)
true<- rep(NA,nrep)


##model1 normal distribution case
miu = (miu_1 + miu_2)/2
se = miu_1 - miu_2
for (s in 1:nrep){
  set.seed(1004145944+s)
  a = sample(1:2, 1, replace=TRUE)
  if(a == 1){
    z = rnorm(1, miu_1, segma_1)
    true[s] = 1
  }
  else if(a ==2){
    z =  rnorm(1, miu_2, segma_1)
    true[s] = 2
  }
  ## LDA
  if(t(z - miu)%*%omega_1%*%se >= 0){
    T1[s] = 1
  }
  else{
    T1[s] = 2
  }
  
  ##LPD
  X = rnorm(n_1, miu_1, segma_1)
  Y = rnorm(n_2, miu_2, segma_1)
  X_BAR = 1/n_1 * (sum(X))
  Y_BAR = 1/n_2 * (sum(Y))
  se_Bar = X_BAR - Y_BAR
  miu_hat = (X_BAR+Y_BAR)/2
  segma_x_hat = 1/n_1 * (sum((X - X_BAR) * t(X - X_BAR)))
  segma_y_hat = 1/n_2 * (sum((Y - Y_BAR) * t(Y - Y_BAR)))
  segma_hat = 1/(n_1+n_2) * (n_1 * segma_x_hat + n_2 * segma_y_hat)
  delta_p =  t(se)%*%omega_1%*%se
  lambda_n = sqrt(delta_p * log(p)/(n_1 = n_2))
  Betahat <- ((- lambda_n + (X_BAR - Y_BAR))/ segma_hat)
  if((z - miu_hat) * Betahat >= 0){
    T2[s] = 1
  }
  else{
    T2[s] = 2
  }


  }
correctness_1_n1 = 0
correctness_2_n1 = 0
for(s in 1:nrep){
  if(T1[s] == true[s]){
    correctness_1_n1 = correctness_1_n1 + 1
  }
  if(T2[s] == true[s]){
    correctness_2_n1 = correctness_2_n1 + 1
  }
}
correctness_1_n1 = correctness_1_n1/ nrep;correctness_1_n1
correctness_2_n1 = correctness_2_n1/ nrep;correctness_2_n1


## Model 2 normal -distribution
for (s in 1:nrep){
  set.seed(1004145944+s)
  a = sample(1:2, 1, replace=TRUE)
  if(a == 1){
    z = rnorm(1, miu_1, segma_2)
    true[s] = 1
  }
  else if(a ==2){
    z =  rnorm(1, miu_2, segma_2)
    true[s] = 2
  }
  ## LDA
  if(t(z - miu)%*%omega_2%*%se >= 0){
    T1[s] = 1
  }
  else{
    T1[s] = 2
  }
  
  ##LPD
  X = rnorm(n_1, miu_1, segma_2)
  Y = rnorm(n_2, miu_2, segma_2)
  X_BAR = 1/n_1 * (sum(X))
  Y_BAR = 1/n_2 * (sum(Y))
  se_Bar = X_BAR - Y_BAR
  miu_hat = (X_BAR+Y_BAR)/2
  segma_x_hat = 1/n_1 * (sum((X - X_BAR) * t(X - X_BAR)))
  segma_y_hat = 1/n_2 * (sum((Y - Y_BAR) * t(Y - Y_BAR)))
  segma_hat = 1/(n_1+n_2) * (n_1 * segma_x_hat + n_2 * segma_y_hat)
  delta_p =  t(se)%*%omega_2%*%se
  lambda_n = sqrt(delta_p * log(p)/(n_1 = n_2))
  Betahat <- ((- lambda_n + (X_BAR - Y_BAR))/ segma_hat)
  if((z - miu_hat) * Betahat >= 0){
    T2[s] = 1
  }
  else{
    T2[s] = 2
  }
  
  
}
correctness_1_n2 = 0
correctness_2_n2 = 0
for(s in 1:nrep){
  if(T1[s] == true[s]){
    correctness_1_n2 = correctness_1_n2 + 1
  }
  if(T2[s] == true[s]){
    correctness_2_n2 = correctness_2_n2 + 1
  }
}
correctness_1_n2 = correctness_1_n2/ nrep;correctness_1_n2
correctness_2_n2 = correctness_2_n2/ nrep;correctness_2_n2
