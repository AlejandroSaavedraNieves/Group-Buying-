bankruptcy_cost_allocations <- function(q, a, b, eta) {
  
  n <- length(q)
  q_N <- sum(q)
  P <- function(x) {a + b / x^eta}
  
  individual_cost <- P(q) * q
  q_minus_i <- q_N - q
  E <- b*(sum(q^(1 - eta)) - q_N^(1 - eta))
  d <- b*(q^(1 - eta) + q_minus_i^(1 - eta) - q_N^(1 - eta))
  
  PROP <- E * d / sum(d)
  PROP_cost <- individual_cost -  PROP
  
  cea_allocation <- function(claims, estate) {
    
    if (estate >= sum(claims)) {return(claims)}
    sorted_claims <- sort(claims)
    lambda <- uniroot(
      function(x) sum(pmin(claims, x)) - estate,
      interval = c(0, max(claims))
    )$root
    return(pmin(claims, lambda))
  }
  
  CEA <- cea_allocation(d, E)
  CEA_cost <- individual_cost -  CEA
  
  
  cel_allocation <- function(claims, estate) {
    
    if (estate >= sum(claims)) {return(claims)}
    lambda <- uniroot(
      function(x) sum(pmax(claims - x, 0)) - estate,
      interval = c(0, max(claims))
    )$root
    return(pmax(claims - lambda, 0))
  }
  
  CEL <- cel_allocation(d, E)
  CEL_cost <- individual_cost -  CEL
  
  
  
  TAL <- numeric(n)
  if (E <= sum(d) / 2) {
    lambda <- uniroot(
      function(x) sum(pmin(d / 2, x)) - E,
      interval = c(0, max(d / 2))
    )$root
    TAL <- pmin(d / 2, lambda)
  } else {
    losses <- sum(d) - E
    lambda <- uniroot(
      function(x) sum(pmin(d / 2, x)) - losses,
      interval = c(0, max(d / 2))
    )$root
    TAL <- d - pmin(d / 2, lambda)
  }
  TAL_cost <- individual_cost - TAL
  
  
  allocations <- data.frame(agent = seq_len(n),
                            quantity = q,
                            individual_cost = individual_cost,
                            claim = d,
                            PROP = PROP_cost,
                            CEA = CEA_cost,
                            CEL = CEL_cost,
                            TAL = TAL_cost
  )
  
  rel_savings<- data.frame(agent = seq_len(n),
                           quantity = q,
                           individual_cost = individual_cost,
                           claim = d,
                           PROPs = (individual_cost-PROP_cost)/individual_cost,
                           CEAs = (individual_cost-CEA_cost)/individual_cost,
                           CELs = (individual_cost-CEL_cost)/individual_cost,
                           TALs = (individual_cost-TAL_cost)/individual_cost
  )
  
  return(list(estate = E,claims = d,
              bankruptcy_allocations = data.frame(agent = seq_len(n),
                                                  PROP = PROP,CEA = CEA,CEL = CEL,TAL = TAL),
              cost_allocations = allocations,
              relative_savings = rel_savings
  ))
}

result <- bankruptcy_cost_allocations(q = c(2, 3, 5),
                                      a = 10, b = 100, eta = 0.5)
result$cost_allocations
