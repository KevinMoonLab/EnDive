function [w0, epsreturn]=calculateWeightsKDE_cvx(N,d,h, varargin)

% Compute the optimal weight for divergence/entropy estimation using CVX
% and a KDE plug-in estimator. Note that CVX must be installed and the
% paths must be set up prior to using this function.

% % % Inputs
%   N = sample size
%   d = dimension
%   h = vector of bandwidths
% varargin:
%   'eta' (default = 1)
%       Determines a tradeoff between bias and variance. Larger eta allows
%       for larger variance and smaller bias.
%   'est' (default = 1)
%       Determines which version of EnDive is used. A value of 1 gives the
%       weights for the estimator given in the main paper below. A value of 2
%       gives the modified estimator in Appendix B. In theory, this second
%       estimator has better convergence rates, although the range of
%       possible bandwidths is more limited.
%   'delta' (default = 0.2)
%       Only used when est == 2. The delta parameter is used to determine
%       the exponent relating the bandwidth to the basis functions for the
%       weights. See Appendix B in the paper below.

% % % Outputs
% w0 = optimal weight
% epsreturn = corresponding value of epsilon in the optimization problem
% 
% Written by Kevin Moon, November 2022
% 
% Relevant Citation: K. Moon, K. Sricharan, K. Greenewald, A.O. Hero III, 
% "Ensemble Estimation of Information Divergence," Entropy (Special Issue 
% on Information Theory in Machine Learning and Data Science), vol. 20, 
% no. 8, pp. 560, July 2018. 

% Set up default parameters
eta=1;
est=1;
delta=.2;

% Get input parameters
for i=1:length(varargin)
%     Eta parameter
    if(strcmp(varargin{i},'eta'))
        eta=varargin{i+1};
    end
%     Which estimator to use.
    if(strcmp(varargin{i},'est'))
        est=varargin{i+1};
        if est~=1&&est~=2
            est=1;
            disp 'Invalid choice of estimator. Using the default value.'
        end
    end
    if(strcmp(varargin{i},'delta'))
        delta=varargin{i+1};
    end
end

% The parameter for the basis function
L=length(h);
if est==1
    lparam=h*N^(1/(2*d));
else
    lparam=h*N^(1/(d+delta));
end

% Basis functions
if est==1    
    psi=zeros(L,d);
    for r=1:d
        psi(:,r)=lparam.^r; 
    end
    Nfunc=N.^(-(1:d)/(2*d));
else
    psi=zeros(L,1);
    Nfunc=0;
    r=1;
    for j=0:floor((d+delta)/2)
        for q=0:floor((d+delta)/2)
            if j+q*delta>0&&j+q*delta<(d+delta)/2
                psi(:,r)=lparam.^(j-d*q);
                Nfunc(r)=N^(-(j+q*delta)/(d+delta));
                r=r+1;
            end
        end
    end
end

% Find the optimal weights. See (6) in the paper.
cvx_begin quiet
    variables w(L) t(1); %t=epsilon
    minimize t
    subject to
      sum(w)==1;
%       norm(w,2)<=eta;
      norm(w,2)<=eta*t;
      max(N^.5*Nfunc.*abs(w'*psi))<=t;
cvx_end

% Check if the weights are not nans
if isnan(sum(w))
    w0=ones(L,1)/L;
    disp 'Optimization returned NaNs. Setting w0=1/L.'
else
    w0=w;
end
epsreturn=t;
