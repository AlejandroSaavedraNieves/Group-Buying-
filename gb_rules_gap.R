group_buying_rules_gap <- function(q, a, b, eta) {
  
  if (any(q <= 0)) {stop("All quantities q_i must be strictly positive.")}
  
  if (a < 0) {stop("a must be positive")}
  
  if (b <= 0) {stop("b must be positive")}
  
  if (eta<=0 ||eta>=1) {stop("eta should satisfy 0<eta<1")}
  
  
  # Schotanus price function
  P <- function(x){a + b / x^eta}
  
  n <- length(q)
  q_N <- sum(q)
  C_GB <- P(q_N) * q_N
  individual_cost <- P(q) * q
  
  q_minus_i <- q_N - q
  C_minus_i <- P(q_minus_i) * q_minus_i
  SC <- C_GB - C_minus_i
  
  NSC <- C_GB - sum(SC)
  
  Equal_Price <- P(q_N) * q
  
  ENSC <- SC + NSC / n
  
  surplus_gap <- individual_cost - SC
  
  denominator <- sum(surplus_gap)
  
  if (abs(denominator) < .Machine$double.eps) {
    SCRB <- rep(NA_real_, n)
    warning("The denominator of the SCRB rule is zero; SCRB is undefined.")
  } else {
    SCRB <- SC + (surplus_gap / denominator) * NSC
  }
  
  
  # Allocation results
  allocation <- data.frame(agent = seq_len(n),quantity = q,
                           individual_cost = individual_cost,separable_cost = SC,
                           equal_price = Equal_Price,ENSC = ENSC,SCRB = SCRB)
  
  return(list(parameters = list(a = a,b = b,eta = eta),
              total_quantity = q_N,
              group_price = P(q_N),
              total_cost = C_GB,
              nonseparable_cost = NSC,
              separable_cost = SC,
              allocation = allocation
  ))
}

result <- group_buying_rules_gap(q = c(2, 3, 5),
                                 a = 10, b = 100, eta = 0.5)
result$allocation