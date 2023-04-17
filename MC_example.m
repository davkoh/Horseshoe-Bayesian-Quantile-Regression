% This is an example run file for the simulations of the HS-BQR
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

mainloc=cd;
folder1 = mainloc+"\functions\";
folder2 = folder1+"\NearestSymmetricPositiveDefinite\";
folder3 = mainloc+"\MC\";
addpath(folder1,folder2,folder3);
clear folder*

%% Generate a random sample using the MC
Sparse_type='sparse' %Sparse_type='sparse','dense','block'
T=100 %Needs to be larger than 100 and divisible by 100.
[data_het,data_nohet,eps,beta]=datagen_func(Sparse_type,T);

%% Select an MC to run
Design_type='y1' %Design_type='y1','y2','y3','y4'
Experiment_num=36 %Experiment_num=1-100

if strcmp(Design_type,'y1')==1
    y=squeeze(data_nohet(1:end,1,Experiment_num));
    X=squeeze(data_nohet(1:end,3:size(data_nohet,2),Experiment_num));
elseif strcmp(Design_type,'y2')==1
    y=squeeze(data_nohet(1:end,2,Experiment_num));
    X=squeeze(data_nohet(1:end,3:size(data_nohet,2),Experiment_num));
elseif strcmp(Design_type,'y3')==1
    y=squeeze(data_het(1:end,2,Experiment_num));
    X=squeeze(data_het(1:end,6:size(data_het,2),Experiment_num));
elseif strcmp(Design_type,'y4')==1
    y=squeeze(data_het(1:end,3,Experiment_num));
    X=squeeze(data_het(1:end,6:size(data_het,2),Experiment_num));
end

%Clear unused variables from memory
clear data_* eps

%% Run HS-BQR
quantile=0.5; %Can be anything between (not including) 0 and 1.
y_fit=y(1:T,1:end);
X_fit=X(1:T,1:end);
%Last 100 observations always used for testing.

Draw_Save=2000;
Burn_in=4000;
thining=1;
disp_iter=1000; %Display when it reaches every disp_iter iteration.

disp("Running HS-BQR");
[HSBQR_temp]=HSBQR(y_fit,X_fit,quantile,Draw_Save,Burn_in,thining,disp_iter);