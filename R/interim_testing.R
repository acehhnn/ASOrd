#' Use proportional odds model approach in the interim analysis
#' 
#' @param K number of tests including the final analysis
#' @param k order of the test among all tests
#' @param y1 outcomes of currently involved members in the new treatment group
#' @param y0 outcomes of currently involved members in the control group
#' @param alpha overall Type I error
#' @param t a sequence of fractions in the interim analysis
#' @param func choosing an alpha spending function
#' @param gamma parameter used in gamma family alpha spending function
#' @return Results of hypothesis testing
#' @export
interim_odds = function(K, k, y1, y0, alpha = 0.05, t = NULL, func = 'obf', gamma = 0){
  if(is.null(t)){
    t = seq(1/K, 1, 1/K)
  }
  if(length(t) != K){
    stop('Length of t must be equal to K.')
  }
  
  if(gamma == 0){
    c = get_boundary(K = K, alpha = alpha, t = t, func = func)
  }else{
    c = get_boundary(K = K, alpha = alpha, t = t, func = func, gamma = gamma)
  }
  
  n1 = length(y1)
  n0 = length(y0)
  
  df = data.frame(
    y = c(y1, y0),
    x = c(rep(1, n1), rep(0, n0))
  )
  df$y = as.factor(df$y)
  
  fit = MASS::polr(y ~ x, df, method = "logistic", Hess = TRUE) 
  beta = fit$coefficients 
  sigma = (solve(fit$Hessian))[1,1] 
  z = beta / sqrt(sigma)
  if(z >= c[k]){
    p_value = 1 - pnorm(z)
    print(paste('Reject H0 with p =', p_value))
  }else{
    print('Not enough evidence to reject H0.')
  }
}



#' Use rank sum statistics in the interim analysis
#' 
#' @param K number of tests including the final analysis
#' @param k order of the test among all tests
#' @param y1 outcomes of currently involved members in the new treatment group
#' @param y0 outcomes of currently involved members in the control group
#' @param A number of levels of the outcomes
#' @param alpha overall Type I error
#' @param t a sequence of fractions in the interim analysis
#' @param func choosing an alpha spending function
#' @param gamma parameter used in gamma family alpha spending function
#' @return Results of hypothesis testing
#' @export
interim_rank = function(K, k, y1, y0, A, alpha = 0.05, t = NULL, func = 'obf', gamma = 0){
  if(is.null(t)){
    t = seq(1/K, 1, 1/K)
  }
  if(length(t) != K){
    stop('Length of t must be equal to K.')
  }
  
  if(gamma == 0){
    c = get_boundary(K = K, alpha = alpha, t = t, func = func)
  }else{
    c = get_boundary(K = K, alpha = alpha, t = t, func = func, gamma = gamma)
  }
  
  n1 = length(y1)
  n0 = length(y0)
  
  y = c(y1, y0)
  
  d = NULL # number of patients in each level
  for(a in 1:A){
    d = append(d, sum(y == a))
  }
  
  r = rank(y)
  w = sum(r[1:n1])
  mu = n1 * (n1 + n0 +1) / 2
  var = n1 * n0 * (n1 + n0 + 1) / 12 - n1 * n0 * sum(d^3 - d) / (12 * (n1 + n0) * (n1 + n0 - 1))
  z = (w - mu) / sqrt(var)
  if(z >= c[k]){
    p_value = 1 - pnorm(z)
    print(paste('Reject H0 with p =', p_value))
  }else{
    print('Not enough evidence to reject H0.')
  }
}