
#####################################
# b^PROP
#####################################

b_PROP <- function(q, a, b, eta) {
  
  n <- length(q)
  qN <- sum(q)
  qi_pow <- q^(1 - eta)
  qN_pow <- qN^(1 - eta)
  q_minus_i <- qN - q
  q_minus_i_pow <- q_minus_i^(1 - eta)
  
  d <- b*(qi_pow +q_minus_i_pow -qN_pow)
  E <- b*(sum(qi_pow)-qN_pow)
  
  sum_d <- sum(qi_pow + q_minus_i_pow - qN_pow)
  
  PROP <- b * (qi_pow + q_minus_i_pow - qN_pow)/sum_d*(sum(qi_pow) - qN_pow)
  
  b_PROP <- b*(1-(qi_pow + q_minus_i_pow - qN_pow)/sum_d*(sum(qi_pow)-qN_pow)/qi_pow)
  
  cost <- a * q + b * qi_pow - PROP
  
  price <- a + b_PROP * q^(-eta)
  
  data.frame(i = 1:n,q = q,d_GB = d, PROP = PROP, b_PROP = b_PROP, cost = cost, price = price)
}

#####################################
# b^EP
#####################################

b_EP <- function(q, a, b, eta) {
  
  n <- length(q)
  qN <- sum(q)
  
  # Equal Price cost allocation
  EP <- a * q + b * q / qN^eta
  
  # Individual coefficient b^EP
  b_EP <- b * (q / qN)^eta
  
  # Individual cost under the Schotanus-type price function
  cost <- a * q + b * q^(1 - eta)
  
  # Effective individual price under EP
  price <- EP / q
  
  data.frame(
    i = 1:n,
    q = q,
    EP = EP,
    b_EP = b_EP,
    cost = cost,
    price = price
  )
}

#####################################
# b^CEA
#####################################
b_CEA <- function(q, a, b, eta) {
  n <- length(q)
  qN <- sum(q)
  qN_pow <- qN^(1 - eta)
  qi_pow <- q^(1 - eta)
  
  q_minus_i <- qN - q
  q_minus_i_pow <- q_minus_i^(1 - eta)
  d <- b * (qi_pow + q_minus_i_pow - qN_pow)
  E <- b * (sum(qi_pow) - qN_pow)
  
  if (min(d) >= E / n) {
    b_CEA <- b * (1 - (sum(qi_pow) - qN_pow) / (n * qi_pow))
    return(list(caso = "a",k = 0,E = E,d = d,b_CEA = b_CEA))
  }
  
  
  orden <- order(d)
  
  d_ord <- d[orden]
  qi_pow_ord <- qi_pow[orden]
  q_minus_i_pow_ord <- q_minus_i_pow[orden]
  
  k <- NA
  E_k <- NA
  
  for (kk in 1:(n - 1)) {
    
    E_kk <- b / (n - kk) * (
      (sum(qi_pow_ord) - qN_pow) -
        sum(
          qi_pow_ord[1:kk] +
            q_minus_i_pow_ord[1:kk] -
            qN_pow
        )
    )
    
    if (d_ord[kk] <= E_kk && E_kk < d_ord[kk + 1]) {
      k <- kk
      E_k <- E_kk
      break
    }
  }
  
  if (is.na(k)) {
    stop("No k satisfying the CEA condition has been found.")
  }
  
  
  b_CEA_ord <- numeric(n)
  b_CEA_ord[1:k] <- b * ((qN_pow - q_minus_i_pow_ord[1:k]) / qi_pow_ord[1:k])
  
  
  
  termino_comun <- (sum(qi_pow_ord[(k + 1):n]) + sum(qN_pow - q_minus_i_pow_ord[1:k]) - qN_pow)
  
  b_CEA_ord[(k + 1):n] <- b * (1 - termino_comun/(qi_pow_ord[(k + 1):n] * (n - k)))
  
  
  b_CEA <- numeric(n)
  b_CEA[orden] <- b_CEA_ord
  
  
  return(list(caso = "b",k = k,E = E,d = d,orden = orden,E_k = E_k,b_CEA = b_CEA))
}
#####################################
# b^CEL
#####################################

b_CEL <- function(q, a, b, eta) {
  
  n <- length(q)
  qN <- sum(q)
  qi_pow <- q^(1 - eta)
  qN_pow <- qN^(1 - eta)
  q_minus_i <- qN - q
  q_minus_i_pow <- q_minus_i^(1 - eta)
  
  d <- b * (qi_pow + q_minus_i_pow - qN_pow)
  E <- b * (sum(qi_pow) - qN_pow)
  
  L <- sum(d) - E
  
  orden <- order(d)
  
  d_ord <- d[orden]
  qi_pow_ord <- qi_pow[orden]
  q_minus_i_pow_ord <- q_minus_i_pow[orden]
  
  
  if (min(d) >= L / n) {
    ajuste_comun <- (qN_pow + sum(q_minus_i_pow - qN_pow)) / n
    
    b_CEL <- b * (qN_pow - q_minus_i_pow + ajuste_comun) / qi_pow
    
    CEL <- d - L / n
    cost <- a * q + b * qi_pow - CEL
    
    # Precio
    price <- a + b_CEL * q^(-eta)
    
    
    return(data.frame(i = 1:n,q = q,d_GB = d,CEL = CEL,b_CEL = b_CEL,cost = cost,
                      price = price,caso = "a",k = 0))
  }
  
  
  k <- NA
  L_k <- NA
  
  for (kk in 1:(n - 1)) {
    L_kk <- b / (n - kk) * (qN_pow + sum(q_minus_i_pow_ord - qN_pow) -
                              sum(qi_pow_ord[1:kk] + q_minus_i_pow_ord[1:kk] - qN_pow))
    
    if (d_ord[kk] <= L_kk && L_kk < d_ord[kk + 1]) {
      k <- kk
      L_k <- L_kk
      break
    }
  }
  
  
  # Comprobación
  if (is.na(k)) {
    stop(
      "No k satisfying the CEL condition has been found."
    )
  }
  
  
  b_CEL_ord <- numeric(n)
  b_CEL_ord[1:k] <- b
  
  suma_k <- sum(qi_pow_ord[1:k])
  
  for (r in (k + 1):n) {
    
    b_CEL_ord[r] <- b * (qN_pow - suma_k + sum(q_minus_i_pow_ord[(k + 1):n] - q_minus_i_pow_ord[r])) / (
      (n - k) *
        qi_pow_ord[r]
    )
  }
  
  
  b_CEL <- numeric(n)
  b_CEL[orden] <- b_CEL_ord
  
  
  
  CEL_ord <- numeric(n)
  
  # i <= k
  CEL_ord[1:k] <- 0
  
  # i > k
  CEL_ord[(k + 1):n] <- d_ord[(k + 1):n] -
    L_k
  
  CEL <- numeric(n)
  CEL[orden] <- CEL_ord
  
  
  cost <- a * q +
    b * qi_pow -
    CEL
  
  
  price <- a +
    b_CEL * q^(-eta)
  
  
  data.frame(
    i = 1:n,
    q = q,
    d_GB = d,
    CEL = CEL,
    b_CEL = b_CEL,
    cost = cost,
    price = price,
    caso = "b",
    k = k
  )
}

#####################################
# b^TAL
#####################################
b_TAL <- function(q, a, b, eta) {
  
  n <- length(q)
  qN <- sum(q)
  
  qi_pow <- q^(1 - eta)
  qN_pow <- qN^(1 - eta)
  
  q_minus_i <- qN - q
  q_minus_i_pow <- q_minus_i^(1 - eta)
  
  
  d <- b * (qi_pow + q_minus_i_pow - qN_pow)
  E <- b * (sum(qi_pow) - qN_pow)
  
  sum_d <- sum(d)
  
  # Order agents by their claims
  orden <- order(d)
  
  d_ord <- d[orden]
  qi_pow_ord <- qi_pow[orden]
  q_minus_i_pow_ord <- q_minus_i_pow[orden]
  
  
  # ============================================================
  # LOW ESTATE
  # ============================================================
  
  if (E <= sum_d / 2) {
    
    k <- NA
    TAL_k <- NA
    
    for (kk in 0:(n - 1)) {
      
      if (kk == 0) {
        
        termino <- E / n
        
        limite_inferior <- 0
        limite_superior <- d_ord[1] / 2
        
      } else {
        
        suma_k <- sum(
          qi_pow_ord[1:kk] +
            q_minus_i_pow_ord[1:kk] -
            qN_pow
        )
        
        termino <- (E - (b / 2) * suma_k) / (n - kk)
        
        limite_inferior <- d_ord[kk] / 2
        limite_superior <- d_ord[kk + 1] / 2
      }
      
      if (limite_inferior <= termino &&
          termino < limite_superior) {
        
        k <- kk
        TAL_k <- termino
        break
      }
    }
    
    if (is.na(k)) {
      stop(
        "No k satisfying the Talmud condition in low estate "
        "has been found."
      )
    }
    
    
    # ------------------------------------------------------------
    # b_TAL
    # ------------------------------------------------------------
    
    b_TAL_ord <- numeric(n)
    
    # i <= k
    if (k > 0) {
      
      b_TAL_ord[1:k] <-
        (b / 2) *
        (
          1 +
            (qN_pow - q_minus_i_pow_ord[1:k]) /
            qi_pow_ord[1:k]
        )
    }
    
    
    # i > k
    if (k < n) {
      
      if (k == 0) {
        suma_k <- 0
      } else {
        suma_k <- sum(
          qi_pow_ord[1:k] +
            q_minus_i_pow_ord[1:k] -
            qN_pow
        )
      }
      
      termino_comun <-
        sum(qi_pow_ord) -
        qN_pow -
        suma_k / 2
      
      b_TAL_ord[(k + 1):n] <-
        b * (
          1 -
            termino_comun /
            ((n - k) * qi_pow_ord[(k + 1):n])
        )
    }
    
    
    # ------------------------------------------------------------
    # TAL allocation
    # ------------------------------------------------------------
    
    TAL_ord <- numeric(n)
    
    if (k > 0) {
      TAL_ord[1:k] <- d_ord[1:k] / 2
    }
    
    if (k < n) {
      TAL_ord[(k + 1):n] <- TAL_k
    }
    
    
    # Return to original ordering
    b_TAL <- numeric(n)
    TAL <- numeric(n)
    
    b_TAL[orden] <- b_TAL_ord
    TAL[orden] <- TAL_ord
    
    
    # Cost and price
    cost <- a * q + b * qi_pow - TAL
    price <- a + b_TAL * q^(-eta)
    
    
    return(
      data.frame(
        i = 1:n,
        q = q,
        d_GB = d,
        TAL = TAL,
        b_TAL = b_TAL,
        cost = cost,
        price = price,
        regime = "low_estate",
        k = k
      )
    )
    
    
  } else {
    
    
    # ============================================================
    # HIGH ESTATE
    # ============================================================
    
    k <- NA
    TAL_k <- NA
    
    for (kk in 0:(n - 1)) {
      
      if (kk == 0) {
        
        termino <-
          b / n *
          (
            sum(q_minus_i_pow - qN_pow) +
              qN_pow
          )
        
        limite_inferior <- 0
        limite_superior <- d_ord[1] / 2
        
      } else {
        
        suma_k <- sum(
          qi_pow_ord[1:kk] +
            q_minus_i_pow_ord[1:kk] -
            qN_pow
        )
        
        termino <-
          b / (n - kk) *
          (
            sum(q_minus_i_pow_ord - qN_pow) +
              qN_pow -
              suma_k / 2
          )
        
        limite_inferior <- d_ord[kk] / 2
        limite_superior <- d_ord[kk + 1] / 2
      }
      
      if (limite_inferior <= termino &&
          termino < limite_superior) {
        
        k <- kk
        TAL_k <- termino
        break
      }
    }
    
    if (is.na(k)) {
      stop(
        "No k satisfying the Talmud condition in high estate "
        "has been found."
      )
    }
    
    
    # ------------------------------------------------------------
    # Common term A
    #
    # A = sum(d) - E - (1/2) sum_{j=1}^k d_j
    # ------------------------------------------------------------
    
    if (k == 0) {
      suma_d_k <- 0
    } else {
      suma_d_k <- sum(d_ord[1:k])
    }
    
    A <- sum_d - E - (1 / 2) * suma_d_k
    
    
    # ------------------------------------------------------------
    # b_TAL
    # ------------------------------------------------------------
    
    b_TAL_ord <- numeric(n)
    
    
    # i <= k
    if (k > 0) {
      
      b_TAL_ord[1:k] <- (b / 2) *
        (1 + 
           (qN_pow - q_minus_i_pow_ord[1:k]) /
           qi_pow_ord[1:k])
    }
    
    
    # i > k
    if (k < n) {
      b_TAL_ord[(k + 1):n] <-
        b -
        d_ord[(k + 1):n] /
        qi_pow_ord[(k + 1):n] +
        A /
        (
          (n - k) *
            qi_pow_ord[(k + 1):n]
        )
    }
    
    
    # ------------------------------------------------------------
    # TAL allocation
    # ------------------------------------------------------------
    
    TAL_ord <- numeric(n)
    
    
    # i <= k
    if (k > 0) {TAL_ord[1:k] <- d_ord[1:k] / 2}
    
    
    # i > k
    if (k < n) {
      TAL_ord[(k + 1):n] <-
        d_ord[(k + 1):n] -
        A / (n - k)
    }
    
    
    # Return to original ordering
    b_TAL <- numeric(n)
    TAL <- numeric(n)
    
    b_TAL[orden] <- b_TAL_ord
    TAL[orden] <- TAL_ord
    
    
    # Cost and price
    cost <- a * q + b * qi_pow - TAL
    price <- a + b_TAL * q^(-eta)
    
    
    return(
      data.frame(
        i = 1:n,
        q = q,
        d_GB = d,
        TAL = TAL,
        b_TAL = b_TAL,
        cost = cost,
        price = price,
        regime = "high_estate",
        k = k
      )
    )
  }
}