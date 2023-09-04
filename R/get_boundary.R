pocock = function(alpha, t){
  f = NULL
  for(i in 1:length(t)){
    f = append(f,alpha*min(log(1 + (exp(1) - 1) * t[i]), 1))
  }
  return(f)
}


obf = function(alpha, t){
  f = NULL
  for(i in 1:length(t)){
    if(t[i] != 0){
      f = append(f, min(2 * (1 - pnorm(qnorm(1 - alpha / 2) / sqrt(t[i]))), alpha))
    }
    else{
      f = append(f,0)
    }
  }
  return(f)
}

gamma_family = function(alpha, gamma, t){
  f = NULL
  if(gamma == 0){
    for(i in 1:length(t)){
      f = append(f, alpha * min(t[i], 1))
    }
  }
  else{
    for(i in 1:length(t)){
      f = append(f, alpha * min(1, (1 - exp(-gamma * t[i])) / (1 - exp(-gamma)), 1))
    }
  }
  return(f)
}



get_boundary = function(K, alpha = 0.05, t = NULL, func= 'obf', gamma = 0){
  if(is.null(t)){
    t = seq(1/K, 1, 1/K)
  }
  if(length(t) != K){
    stop('Length of t must be equal to K.')
  }
  
  c = NULL
  S = matrix(rep(0, K * K), nrow = K) 
  for(i in 1:K){
    for(j in 1:K){
      S[i,j] = min(t[i], t[j]) / max(t[i], t[j])
    }
  }
  mean = rep(0, K)
  sigma = S
  n = 1e+4
  set.seed(20230810)
  simu = MASS::mvrnorm(n, mean, sigma) 
  
  if(func == 'gamma_family'){
    a = gamma_family(alpha, gamma, t)
  }else if(func == 'obf'){
    a = obf(alpha,t)
  }else{
    a = pocock(alpha, t)
  }
    for(k in 1:length(t)){ 
    if(k == 1){
      c0 = qnorm(1 - a[k])
      c = append(c, c0)
    }
    else{
      delta = a[k] - a[k - 1]
      itr = seq(c[k - 1], 0, -0.01) 
      for(d in itr){
        s = 0 
        for(t in 1:n){
          flag = 1
          for(j in 1:(k-1)){
            if(simu[t,j] >= c[j]){
              flag = 0
              break
            }
          }
          if(flag == 1){
            if(simu[t,k] <= d){
              flag = 0
            }else{
              s = s + 1
            }
          }
        }
        prob = s / n
        if(prob >= delta){
          c = append(c, d)
          break
        }
      }
    }
  }
  
  return(c)
}