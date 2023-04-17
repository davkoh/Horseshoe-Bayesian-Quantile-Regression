function [eps,beta, X]=MC_data_indepvar(T,beta_sparse,rho,error_struct)

% Support function for horseshoe prior Bayesian quantile Regression (2020)
%
% See:
% Kohns, D.E. and Szendrei, T. (2020). Horseshoe Prior Bayesian Quantile Regression, arXiv preprint arXiv:2006.07655


% MC_data_indepvar Monte Carlo independent data generation (primarily) for 
% quantile regression
%
%   The follwing OUTPUT arguments are generated:
%
%   1. y is the vector of dependent variable
%
%   2. X is the matrix of explanatory variables
%
%   The following INPUT arguments are used for the function:
%
%   1. T is the sample size
%
%
%   2. beta_sparse determines beta sparsity
%       -'sparse': sparse data matrix
%       -'dense': dense data matrix
%       -'block': block sparse data matrix
%
%   3. rho determines the "correlation" used to generate the Sigma for
%   generating the X matrix
%
%   4. X_struct
%       -'not_het': Positive and negative values for X
%       -'het':Positive values for X only
%
%   5. error_struct determines the epsilon structure
%       -'normal': N(0,1)
%       -'student': Student-t(3)
%
%
%Note: In its current setup all diagonal elements of Sigma are 1 and the
%means of the variables are 0.

%% -------------------- Generate Beta
if strcmp(beta_sparse,'sparse')==1
    beta=[1;1;1/2;1/3;1/4;1/5;zeros(T*2,1)];
elseif strcmp(beta_sparse,'dense')==1
    beta=[1;0.85*ones(T,1)];
elseif strcmp(beta_sparse,'block')==1
    beta=[1;0.85*ones(T,1);zeros(T,1);0.85*ones(T,1)];
else
    print('Please specify a valid sparsity: sparse, dense, block');
    pkill
end

K=size(beta,1)-1;

%% -------------------- Generate Variables
Sigma=zeros(K,K);
mu=zeros(K,1);

for i=1:K
    for j=1:K
        Sigma(i,j)=rho^(abs(i-j));
    end
end

X=mvnrnd(mu,Sigma,T);

X=[ones(T,1) X];


%% -------------------- Generate epsilon

if strcmp(error_struct,'normal')==1
    eps=randn(T,1);
elseif strcmp(error_struct,'student')==1
    eps=trnd(3,[T,1]);
else
    print('Please specify a valid error structure: normal, student')
    pkill
end
end