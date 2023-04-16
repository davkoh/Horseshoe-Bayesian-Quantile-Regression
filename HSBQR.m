function [beta_out]=HSBQR(y,X,quantile,nsave,nburn,thin,iter)
% HSBQR is a function that applies the Horseshoe prior to the Bayesian
% Quantile Regression. Please cite the authors (David Kohns and Tibor
% Szendrei) paper (found in zipped file) in your work.The function
% generates a vector of betas which is an average for all draws.
%
% The input arguments are the following:
% 1. y=LHS variable
% 2. X=Matrix of explanatory variables
% 3. quantile=set the quantile (between 0 and 1)
% 4. nsave=Number of simulations (after burn-in)
% 5. nburn=Number of burn-in runs
% 6. thin=thining 
% 7. iter= print every ith iteration


ntot = nsave + nburn;   % Number of total draws
T=size(X,1);
y = y;
x = X;
q=quantile;
[n,p] = size(x);

% ==============| Define priors
V = 9*eye(p);
Vinv = inv(V);

N=nburn+nsave;
effsamp=(N-nburn)/thin;

%% output %%
betaout=zeros(p,effsamp);
lambdaout=zeros(p,effsamp);
tauout=zeros(effsamp,1);
sigmaSqout=zeros(effsamp,1);

%% matrices %%
I_n=eye(n);
l=ones(n,1);

% prior for sigma2 ~ IG(a0,b0)
a0 = 0.1;
b0 = 0.1;

% ==============| Initialize vectors
beta = zeros(p,1);
z = ones(n,1);
sigma = ones(1);
theta = zeros(1);
tau_sq = zeros(1);

%% paramters %%
Beta=zeros(p,1);
lambda=ones(p,1);
tau=1;

% ==============| Storage matrices
beta_draws = zeros(p,nsave);

for irep = 1:ntot
    % Print every "iter" iterations on the screen
    if mod(irep,iter)==0
        disp(irep)
    end

    tau_sq = 2/(q*(1-q));
    theta = (1-2*q)/(q*(1-q));
    
    
    %% HS prior Implementation
    
    % Sample regression coefficients beta
    U = diag(1./(sqrt(sigma).*tau_sq.*z));
    try
        y_tilde = chol(U)*(y-theta*z);
        X_tilde = chol(U)*x;
    catch
        y_tilde = chol(nearestSPD(U))*(y-theta*z);
        X_tilde = chol(nearestSPD(U))*x;
    end
    
    lambda_star=tau*lambda;
    U_bar=bsxfun(@times,(lambda_star.^2),X_tilde');
    
    u=normrnd(0,lambda_star);
    v=X_tilde*u+normrnd(0,l);
    
    v_star=(X_tilde*U_bar+I_n)\(y_tilde-v);
    beta=(u+U_bar*v_star);
    
    %% update lambda_j's in a block using slice sampling %%
    eta = 1./(lambda.^2);
    upsi = unifrnd(0,1./(1+eta));
    tempps = beta.^2/(2*tau^2);
    ub = (1-upsi)./upsi;
    
    % now sample eta from exp(tempv) truncated between 0 & upsi/(1-upsi)
    Fub = 1 - exp(-tempps.*ub); % exp cdf at ub
    Fub(Fub < (1e-4)) = 1e-4;  % for numerical stability
    up = unifrnd(0,Fub);
    eta = -log(1-up)./tempps;
    lambda = 1./sqrt(eta);
    
    %% update tau %%
    tempt = sum((beta./lambda).^2)/(2);
    et = 1/tau^2;
    utau = unifrnd(0,1/(1+et));
    ubt = (1-utau)/utau;
    Fubt = gamcdf(ubt,(p+1)/2,1/tempt);
    Fubt = max(Fubt,1e-8); % for numerical stability
    ut = unifrnd(0,Fubt);
    et = gaminv(ut,(p+1)/2,1/tempt);
    tau = 1/sqrt(et);
    
    % Sample regression variance sigma2
    a1 = a0 + 3*n/2;
    sse = (y-x*beta - theta*z).^2;
    a2 = b0 + sum(sse./(2*z*tau_sq)) + sum(z);
    sigma = 1./gamrnd(a1,1./a2);
    
    % Sample latent variables z_{t}
    for t = 1:n
        k1 = sqrt(theta.^2 + 2*tau_sq)/abs(y(t,:)-x(t,:)*beta);
        k2 = (theta.^2 + 2*tau_sq)/(sigma*tau_sq);
        z(t) = max(1./Draw_IG(k1,k2),1e-4);
    end
    %end
    
    if irep > nburn
        beta_draws(:,irep-nburn) = beta;
    end
end
beta_out=mean(beta_draws,2);
end