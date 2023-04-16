function [data_het,data_nohet,eps,beta]=datagen_func(Sparse_type,datasize)
%Sparse_type='sparse','dense','block'

%% Monte Carlo Settings
%Dimensions of betas: (varibs,quantiles,y varib used,iterations)

%DON'T CHANGE THESE PARAMETERS!
rng(1994,'twister')
T=200; %Sample size of each Monte-Carlo.
rho=0.5; %Correlation of the X matrix
iters=1000; %Number of Monte Carlo iterations
ynum=7; %This is the number of type of y
%% Monte Carlo
data_het=[];
data_nohet=[];
eps=[];


for j=1:iters
    [eps_normal]=MC_data_indepvar(T,Sparse_type,rho,'normal');
    [eps_student,beta,X]=MC_data_indepvar(T,Sparse_type,rho,'student');
    eps_uniform=rand(T,1)*2;
    
    X_het=X;
    X_nohet=X;
    X_het(:,2)=exp(X(:,2));
    X_het(:,7)=exp(X(:,7));
    
    y=zeros(T,ynum);
    y(:,1)=X_nohet*beta+eps_normal; %y1
    y(:,2)=X_nohet*beta+eps_student; %y2
    y(:,3)=X_het*beta+eps_normal;
    y(:,4)=X_het*beta+(1+X_het(:,2)).*eps_normal; %y3
    y(:,5)=X_het*beta+eps_normal+(X_het(:,2)).*eps_uniform; %y4
    y(:,6)=X_het*beta+(1+X_het(:,2)).*eps_normal+X_het(:,7).*eps_uniform;
    y(:,7)=X_het*beta+eps_student+X_het(:,2).*eps_uniform;
    
    y_het=[y(:,3), y(:,4), y(:,5), y(:,6), y(:,7)];
    y_nohet=[y(:,1), y(:,2)];
    
    data_het_iter=[y_het, X_het];
    data_nohet_iter=[y_nohet, X_nohet];
    
    eps_iter=[eps_normal,eps_student,eps_uniform];
    
    data_het(:,:,j)=data_het_iter;
    data_nohet(:,:,j)=data_nohet_iter;
    eps(:,:,j)=eps_iter;
end

%Reshape array for longer runs (Inefficient way to do things. Only reason it's
%done this way is because T=100 MC's were done first and we wanted to keep
%runs comparable. For people reading this code: I am deeply sorry for doing it like this)
data_het = datareshape(data_het,datasize);
data_nohet = datareshape(data_nohet,datasize);
eps = datareshape(eps,datasize);

end

