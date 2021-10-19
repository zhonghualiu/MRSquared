#' Performs multiply robust Mendelian randomization analysis with many invalid instruments
#'
#' @param Y length n continuous outcome vector
#' @param A length n continuous exposure vector
#' @param G.mat n by K binary SNP matrix
#' @param k 1<=k<=K, the number of valid IVs
#'
#' @return Point estimate of causal effect of A on Y and standard error
#'
#' @export
MRSq = function(Y,A,G.mat,k=5){
  G.hat<- apply(G.mat,2,mean) ## MAF
  ## Generate Z, the new IV matrix
  K = dim(G.mat)[2]
  n = dim(G.mat)[1]
  Z<- numeric();
  for (L in  (K-k+1):K) {
    print("Creating the Z matrix")
    if (L==(K-k+1)) {
      S<- t(combn(1:K,L))
      Z.L<- matrix(nrow=n,ncol=dim(S)[1])
      for (j in 1:dim(S)[1]) {
        if (length(S[j,])==1) {
          Z.L[,j]<- G.mat[,S[j,]]-G.hat[S[j,]];
        } else {
          Z.L[,j]<- apply(sweep(G.mat[,S[j,]],2,G.hat[S[j,]]),1,prod);
        }
      }
    } else {
      S<- t(combn(1:K,L))
      Z.L<- matrix(nrow=n,ncol=dim(S)[1])
      for (j in 1:dim(S)[1]) {
        k.L  <- S[j,];
        if (K==k) {
          P.L <- 1;
        } else if ((K-k)==1) {
          P.L  <- G.mat[,k.L[1:(K-k)]]-G.hat[k.L[1:(K-k)]];
        } else {
          P.L  <- apply(sweep(G.mat[,k.L[1:(K-k)]],2,G.hat[k.L[1:(K-k)]]),1,prod)
        }
        Z.L[,j]<- P.L*(apply(G.mat[,k.L[(K-k+1):L]],1,prod)-prod(G.hat[k.L[(K-k+1):L]]));
      }
    }
    Z<- cbind(Z,Z.L);
  }
  print("The Z matrix created!, running regressions")
  print(dim(Z))
  tsls.out <- summary(ivreg(Y ~ A | Z), vcov = sandwich::sandwich, diagnostics = F)
  print(tsls.out)
  est <- c(tsls.out$coef[2,1], tsls.out$coef[2,2])
  return(est)
}
