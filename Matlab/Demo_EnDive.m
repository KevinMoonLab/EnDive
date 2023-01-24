% A demo for the EnDive estimator. The Renyi-alpha divergence integral between two
% truncated Gaussian distributions is estimated multiple times and averaged. Requires
% the functions TruncatedGaussian.m and rmvnrnd.m to be in the filepath. 

% This file contains a demonstration using the default settings, the second
% estimator, a different kernel, and a different distance function.

%% Setup
% Dimension
d=10;

% Gaussian means
mu1=.7*ones(d,1);
mu2=.3*ones(d,1);

% Covariance matrix
sig=.4*eye(d);

% Truncate samples to the unit cube
xlo=zeros(d,1);
xhi=1*ones(d,1);

% Sample size
N=1000;

% Number of trials
S=50;

% Alpha parameter
a=.5;

% Compute the true value of the divergence integral
invf1=inv(sig);
invf2=inv(sig);
s1=a*invf1+(1-a)*invf2;
s2=(a*mu1'*invf1+(1-a)*mu2'*invf2)';
s3=a*mu1'*invf1*mu1+(1-a)*mu2'*invf2*mu2;
A=det(invf1)^(a/2)*det(invf2)^((1-a)/2)/sqrt(det(s1));
A=A*exp(-.5*s3);
B=exp(-.5*s2'*s1^-1*s2);
divgauss=A/B;

lambda=inv(s1);
gcdf=mvncdf(xlo',xhi',mu1',sig);
fcdf=mvncdf(xlo',xhi',mu2',sig);

divtrunc=divgauss*mvncdf(xlo',xhi',(s1\s2)',lambda)/(gcdf^a*fcdf^(1-a));

divtrue=divtrunc;

%% Divergence integral estimation with default settings

tic
div_est=zeros(S,1);

disp 'Running simulation with default settings'
% Change this to a for loop if parallel processing isn't an option for you
parfor s=1:S
%     Generate the data
    Y=rmvnrnd(mu1,sig,N,[eye(d); -eye(d)],[xhi; -xlo]);
    X=rmvnrnd(mu2,sig,N,[eye(d); -eye(d)],[xhi; -xlo]);

%     Run the estimator using default parameters
    div_est(s)=EnDive(X,Y,'type','renyi','alpha',a,'quiet');
    
end

disp 'Simulation completed!'

% Mean estimate
div_mean=mean(div_est);

% Estimated MSE
div_mse=mean((div_est-divtrue).^2);

% Print results
fprintf('Mean estimate: %0.4f \n',div_mean)
fprintf('Estimated MSE: %0.4f \n',div_mse)

toc

%% Divergence integral estimation using the second estimator

tic
div_est=zeros(S,1);

disp 'Running simulation using the second estimator'
% Change this to a for loop if parallel processing isn't an option for you
parfor s=1:S
%     Generate the data
    Y=rmvnrnd(mu1,sig,N,[eye(d); -eye(d)],[xhi; -xlo]);
    X=rmvnrnd(mu2,sig,N,[eye(d); -eye(d)],[xhi; -xlo]);

%     Run the estimator using default parameters
    div_est(s)=EnDive(X,Y,'type','renyi','alpha',a,'quiet','est',2);
    
end

disp 'Simulation completed!'

% Mean estimate
div_mean=mean(div_est);

% Estimated MSE
div_mse=mean((div_est-divtrue).^2);

% Print results
fprintf('Mean estimate: %0.4f \n',div_mean)
fprintf('Estimated MSE: %0.4f \n',div_mse)

toc

%% Divergence integral estimation with uniform kernel

tic
div_est=zeros(S,1);

disp 'Running simulation with the uniform kernel'
% Change this to a for loop if parallel processing isn't an option for you
parfor s=1:S
%     Generate the data
    Y=rmvnrnd(mu1,sig,N,[eye(d); -eye(d)],[xhi; -xlo]);
    X=rmvnrnd(mu2,sig,N,[eye(d); -eye(d)],[xhi; -xlo]);

%     Run the estimator using default parameters
    div_est(s)=EnDive(X,Y,'type','renyi','alpha',a,'quiet','kernel','uniform');
    
end

disp 'Simulation completed!'

% Mean estimate
div_mean=mean(div_est);

% Estimated MSE
div_mse=mean((div_est-divtrue).^2);

% Print results
fprintf('Mean estimate: %0.4f \n',div_mean)
fprintf('Estimated MSE: %0.4f \n',div_mse)

toc

%% Divergence integral estimation using the uniform kernel with the Chebychev distance

tic
div_est=zeros(S,1);

disp 'Running simulation with the uniform kernel and the Chebychev distance'
% Change this to a for loop if parallel processing isn't an option for you
parfor s=1:S
%     Generate the data
    Y=rmvnrnd(mu1,sig,N,[eye(d); -eye(d)],[xhi; -xlo]);
    X=rmvnrnd(mu2,sig,N,[eye(d); -eye(d)],[xhi; -xlo]);

%     Run the estimator using default parameters
    div_est(s)=EnDive(X,Y,'type','renyi','alpha',a,'quiet','kernel','uniform','distance','chebychev');
    
end

disp 'Simulation completed!'

% Mean estimate
div_mean=mean(div_est);

% Estimated MSE
div_mse=mean((div_est-divtrue).^2);

% Print results
fprintf('Mean estimate: %0.4f \n',div_mean)
fprintf('Estimated MSE: %0.4f \n',div_mse)

toc
