
beta_OR_transform <- function(beta, se){
  OR = exp(beta)
  LL_OR = exp(beta-1.96*se)
  UL_OR = exp(beta+1.96*se)
  output <- cbind(OR, LL_OR, UL_OR)
  return(output)
}