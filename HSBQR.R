HSBQR <- function(y,X,quant,nsave,nburn,thin,iter){
  # HSBQR is a function that applies the Horseshoe prior to the Bayesian
  # Quantile Regression. Please cite the authors (David Kohns and Tibor
  # Szendrei) paper (found in zipped file) in your work.The function
  # generates a vector of betas which is an average for all draws.

  
  # This code is free to use for academic purposes only, provided that the 
  # paper is cited as:
  #
  # Kohns, D.E. and Szendrei, T. (2020). Horseshoe Prior Bayesian Quantile Regression, arXiv preprint arXiv:2006.07655
  #
  # This code comes without technical support of any kind. It is expected to
  # reproduce the results reported in the paper. Under no circumstances will
  # the authors be held responsible for any use (or misuse) of this code in
  # any way.
  
    
  # The input arguments are the following:
  # 1. y=LHS variable                                                    
  # 2. X=Matrix of explanatory variables
  # 3. quantile=set the quantile (between 0 and 1)
  # 4. nsave=Number of simulations (after burn-in)
  # 5. nburn=Number of burn-in runs
  # 6. thin=thining
  # 7. iter= print every ith iteration
  
  #Load required packages
  require(pracma)
  
  
  #Preamble
  ntot <- nsave + nburn# Number of total draws
  
  y <- as.matrix(y)
  x <- X
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  # Define Priors
  V = 9*diag(p)
  Vinv = solve(V)
  
  N=nburn+nsave
  effsamp=(N-nburn)/thin
  
  #Output
  betaout=matrix(0,p,effsamp)
  lambdaout=matrix(0,p,effsamp)
  tauout=matrix(0,effsamp,1)
  sigmaSqout=matrix(0,effsamp,1)
  
  #Matrices
  I_n=diag(n)
  l=matrix(1,n,1)
  
  #Prior for sigma2 ~ IG(a0,b0)
  a0 = 0.1
  b0 = 0.1
  
  #Initialise vectors
  n_q=length(quant)
  beta=matrix(0,p,n_q)
  z = matrix(1,n,n_q)
  sigma = matrix(1,1,n_q)
  theta = matrix(0,1,n_q)
  tau_sq = matrix(0,1,n_q)
  
  #Parameters
  Beta=matrix(0,p,1)
  plambda=matrix(1,p,n_q)
  ptau=matrix(1,1,n_q)
  
  #Storage Matrices
  beta_draws = array(0,c(p,n_q,nsave))
  
  for (irep in 1:ntot){
    if (irep%%iter==0){
      print(irep)
    }
    for (q in 1:n_q){
      tau_sq[,q] = 2/(quant[q]*(1-quant[q]))
      theta[,q] = (1-2*quant[q])/(quant[q]*(1-quant[q]))
      lambda=as.matrix(plambda[,q])
      tau=ptau[,q]
      
      #HS prior implementation
      #Sample regression coefficients beta
      U = diag( 1/(sqrt(sigma[1,q])*tau_sq[,q]*z[,q]) )
      
      y_tilde = Chol_of_U(U)%*%(y-theta[,q]*z[,q])
      X_tilde = Chol_of_U(U)%*%x
      
      lambda_star=tau*lambda
      
      U_bar=(repmat(lambda_star^2,1,n)*t(X_tilde))
      
      u=as.matrix(rnorm(length(lambda_star),0,lambda_star))
      v=X_tilde%*%u+rnorm(n,0,l)
      
      v_star=mldivide((X_tilde%*%U_bar+I_n),(y_tilde-v))
      beta[,q]=(u+U_bar%*%v_star)
      
      #Update lambda_j's in a block using slice sampling
      eta=1/lambda^2
      upsi=as.matrix(runif(length(eta),0,1/(1+eta)))
      tempps=as.matrix(beta[,q]^2)/(2*tau^2)
      ub=(1-upsi)/upsi
      
      Fub=1-exp(-tempps*ub)
      Fub[Fub<(1e-4)]=1e-4
      up=as.matrix(runif(length(Fub),0,Fub))
      eta=-log(1-up)/tempps
      lambda=1/sqrt(eta)
      
      #Update tau
      tempt=sum((beta[,q]/lambda)^2)/2
      et=1/tau^2
      utau=runif(length(et),0,1/(1+et))
      ubt=(1-utau)/utau
      Fubt=pgamma(ubt,(p+1)/2,tempt)
      Fubt=max(Fubt,1e-8)
      ut=runif(length(Fubt),0,Fubt)
      et=qgamma(ut,(p+1)/2,tempt)
      tau=1/sqrt(et)
      
      #Sample regression variance sigma2
      a1=a0+3*n/2
      sse=(y-x%*%as.matrix(beta[,q])-as.matrix(theta[,q]*z[,q]))^2
      a2=b0+sum(sse/(2*z[,q]*tau_sq[,q]))+sum(z[,q])
      sigma[1,q]=1/rgamma(length(a1),a1,a2)
      
      #Sample latent variables z_t
      for (t in 1:n){
        k1=sqrt(theta[,q]^2+2*tau_sq[,q])/abs(y[t,]-x[t,]%*%beta[,q])
        k2=(theta[,q]^2+2*tau_sq[,q])/(sigma[1,q]*tau_sq[,q])
        z[t,q]=max(1/Draw_IG(k1,k2),1e-4)
      }
      plambda[,q]=lambda
      ptau[,q]=tau
    }
    if (irep>nburn){
      beta_draws[,,irep-nburn]=beta
    }
  }
  beta_out=apply(beta_draws,c(1,2),mean)
  return(beta_out)
}

Chol_of_U <-function(U){
  tryCatch({
    chol(U)
    return(chol(U))
  },error=function(e){
    return(chol(nearest_spd(U)))
  }
  )
}

Draw_IG <- function(mu,lambda){
  v0=rnorm(1,0,1)^2
  x1=mu+(.5*(mu^2)*v0)/lambda - (.5*mu/lambda)*sqrt(4*mu*lambda*v0 + (mu^2)*(v0^2))
  x2=(mu^2)/x1
  p1_v0=mu/(mu+x1)
  rand=runif(1,0,1)
  if (rand>p1_v0){
    y=x1
  }else{
    y=x2
  }
  return(y)
}