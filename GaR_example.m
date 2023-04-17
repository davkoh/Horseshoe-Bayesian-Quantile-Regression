% This is the main run file for the empirical application of the HS-BQR
% Model "Horseshoe Prior Bayesian Quantile Regression" by D. Kohns and T.
% Szendrei (2020)
%
% This code is free to use for academic purposes only, provided that the 
% paper is cited as:
%
% Kohns, D.E. and Szendrei, T. (2020). Horseshoe Prior Bayesian Quantile Regression, arXiv preprint arXiv:2006.07655
%
% This code comes without technical support of any kind. It is expected to
% reproduce the results reported in the paper. Under no circumstances will
% the authors be held responsible for any use (or misuse) of this code in
% any way.



clear all
clc

currentdir=cd;
folder = strcat(currentdir,'/functions'); % Needs some chaning
folder2= strcat(folder,'/functions/NearestSymmetricPositiveDefinite'); % Needs some changing
folder3=strcat(currentdir,'/GaR application'); % Needs some changing
addpath(folder,folder2,folder3);

run('dataforHSBQR.m')

tau=[0.1 0.3 0.5 0.7 0.9]';

%% Parellelised Loop HSBQR - wideDB

% -- Forecast Preliminaries --
Tmax=size(wideDB,1); 
tin=50; 
tf=Tmax-tin; % out-of-sample period
% Number of rolling forecasts
experiments=tf;
nfor = 1; % forecast horizons

% -- Transform and Shift Forward -- 
y=wideDB(:,1);
yf = zeros(size(y,1),1);
    for t = 1: size(y,1)-nfor
        yf(t)=sum(y(t+1:t+nfor));
    end
% Divide by nfor unless level is being forecast
yf = yf/nfor;
% Now correct observations after taking lagged values
y = yf(1:end-nfor,:);
Yraw = y;

% -- Shift X Matrix Back --
Xraw = wideDB(:,2:end);
Xraw = Xraw(1:end-nfor,:);
X=Xraw;

% Storage Matrix
pbeta_HSBQR_wide=zeros(size(wideDB(1:tin,2:end),2),size(tau,1),experiments-1);


% Execute HSBQR Function
parfor w = 1:experiments-nfor
    [HSBQR_wide]=HSBQR(y(1:tin+w,:)*100,X(1:tin+w,:),tau,5000,5000,1,1000);
    pbeta_HSBQR_wide(:,:,w)=HSBQR_wide;
end

%% Forecast Evaluation Wide
tin=50;
n_q=size(tau,1);
% HS-BQR
res=zeros(experiments,n_q);
fit_wide=zeros(experiments,n_q);
yf=[];
%Yraw=wideDB(:,1);
Xraw=wideDB(:,2:end);

for i=1:experiments-nfor-1
    testb=pbeta_HSBQR_wide(:,:,i);
    y=Yraw(1:tin+i+1,:);
    X=Xraw(1:tin+i+1,:);
    %X=[ones(size(y,1),1) data(1:tin+i,:)];
    fitt=zeros(1,n_q);
    rest=zeros(1,n_q);
    for j=1:n_q
        betat=testb(:,j);
        fitt=X*betat;
        fit_wide(i,j)=fitt(end,:);
        et=y(end,:)-fitt(end,:);
        rest(:,j)=et;
    end
    res(i,:)=rest;
    yf=[yf;y(end,:)];
end
fit_wide=fit_wide/100;
rmsfe_hsbqr=sqrt(sum(res.^2)./experiments);
fit_hsbqr_wide=fit_wide;