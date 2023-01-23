function [w0, epsreturn]=calculateWeightsKDE(N,d,h,varargin)

% Compute the optimal weight for divergence/entropy estimation using
% a KDE plug-in estimator.

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
%     The delta parameter associated with the second estimator
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

% The sum-to-one equality constraint
Aeq=ones(1,L);
Aeq(end+1)=0;
Beq=1;

% The inequality constraints with the basis functions
if est==1
    A=zeros(d,L);
    for r=1:d
        A(r,:)=lparam.^r*N^(.5-r/(2*d));
    end
    A=[A;-A];
    A(:,L+1)=-ones(2*d,1);
    B=zeros(2*d,1);
else
    A=zeros(1,L);
    r=1;
    for j=0:floor((d+delta)/2)
        for q=0:floor((d+delta)/2)
            if j+q*delta>0&&j+q*delta<(d+delta)/2
                A(r,:)=lparam.^(j-d*q)*N^(.5-(j+q*delta)/(d+delta));
                r=r+1;
            end
        end
    end
    A=[A;-A];
    [numParams,~]=size(A);
    A(:,L+1)=-ones(numParams,1);
    B=zeros(numParams,1);
end

% Initialization
winit=ones(L,1)/L;
winit(end+1)=100;

% fmincon options
options=optimset('Algorithm','interior-point','MaxFunEvals',100000,'MaxIter',10000,'Display','off');

% Optimization
w0=fmincon(@(w) w(end),winit,A,B,Aeq,Beq,[],[],@(w) wnorm(w,eta),options);
epsreturn=w0(end);
w0(end)=[];
w0=w0/sum(w0);

% Function used for the constraints on the l2 norm
function [c, ceq]=wnorm(w,eta)
c=norm(w(1:end-1))-eta*w(end);
ceq=[];
