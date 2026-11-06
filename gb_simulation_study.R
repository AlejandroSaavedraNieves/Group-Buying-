#################################################
# checking monotonicity
################################################
is.monotonic<-function(x){
  cond<-0
  if (all(diff(x) < 0)) cond<-1 # Decreasing
  if (all(diff(x) > 0)) cond<-2 # Increasing
  return(cond)}

#################################################
# simulation study
################################################
n=10
set.seed(1)

nsim<-1000
simulation_results<-matrix(0,nrow=nsim,ncol=5)
simulation_results_savings<-matrix(0,nrow=nsim,ncol=5)
simulation_results_allocations<-matrix(0,nrow=nsim,ncol=5)
sim_dif<-rep(0,nrow=nsim)

for (i in 1:nsim){
  set.seed(i)
  
  x0i<-sample(1000:3000,size=n,replace=F)
  a0i<-round(runif(1,0,1),3)*10
  b0i<-round(runif(1,0.01,1),3)*10
  eta0<-round(runif(1,0.01,0.5),3)
  q0i<-sort(x0i,decreasing=F)
  
  ENSC<-group_buying_rules_gap(q0i, a=a0i, b=b0i, eta=eta0)$allocation[,6]
  
  bepaux<-b_EP(q0i, a=a0i, b=b0i, eta=eta0)
  bep<-bepaux[,4]
  bpi<-b_PROP(q0i, a=a0i, b=b0i, eta=eta0)[,5]
  bceai<-b_CEA(q0i, a=a0i, b=b0i, eta=eta0)$b_CEA
  bceli<-b_CEL(q0i, a=a0i, b=b0i, eta=eta0)[,5]
  btali<-b_TAL(q0i, a=a0i, b=b0i, eta=eta0)[,5]
  
  res<-bankruptcy_cost_allocations(q0i, a=a0i, b=b0i, eta=eta0)
  res_sav<-res$relative_savings
  res_all<-res$cost_allocations
  
  simulation_results[i,1]<-is.monotonic(bpi)
  simulation_results[i,2]<-is.monotonic(bpi)
  simulation_results[i,3]<-is.monotonic(bceai)
  simulation_results[i,4]<-is.monotonic(bceli)
  simulation_results[i,5]<-is.monotonic(btali)
  
  simulation_results_allocations[i,1]<-is.monotonic(bepaux[,3])
  simulation_results_allocations[i,2]<-is.monotonic(res_all[,5])
  simulation_results_allocations[i,3]<-is.monotonic(res_all[,6])
  simulation_results_allocations[i,4]<-is.monotonic(res_all[,7])
  simulation_results_allocations[i,5]<-is.monotonic(res_all[,8])
}


#0: ND; 1:D; 2:I
# Allocations
table(simulation_results_allocations[,1])
table(simulation_results_allocations[,2])
table(simulation_results_allocations[,3])
table(simulation_results_allocations[,4])
table(simulation_results_allocations[,5])

#b
table(simulation_results[,1])
table(simulation_results[,2])
table(simulation_results[,3])
table(simulation_results[,4])
table(simulation_results[,5])
