function [div_est,w,h]=EnDive(X,Y,varargin)
% Estimates the f-divergence functional from points in X and Y, i.e. the
% quantity \int{f(q(x)/p(x))*p(x)} using the EnDive estimator. Points in X 
% are drawn from p and points in Y are drawn from q. The estimator is an 
% optimally weighted ensemble estimator using KDEs. Requires the function
% calculateWeightsKDE.m or calculateWeightsKDE_cvx.m as well.

% % Inputs
%   X = points from p and is a Nxd matrix
%   Y = points from q and is a Mxd matrix
% 
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
%       weights. See Appendix B in the paper below. Must be greater than
%       zero.
%   'h' (default = chosen based on k-nn distances)
%       A length L vector containing the set of bandwidths to use in the
%       ensemble of estimators. If not given, then the default is to select
%       the bandwidths based on the k-nn distances of the data.
%   'kmax' (default = chosen as a function of N)
%       Used for automatically computing the bandwidths. This is the 
%       k-value used to get the maximum value of the bandwidth range. See
%       'kmin' below for the minimum. Must be an integer greater than zero.
%   'kmin' (default = chosen as a function of N)
%       Used for automatically computing the bandwidths. This is the 
%       k-value used to get the minimum value of the bandwidth range. See
%       'kmax' above for the maximum. Must be an integer greater than zero.
%   'L' (default = 40)
%       Only used when automatically computing the bandwidths. This gives
%       the number of bandwidth values (and thus base estimators) in
%       between the minimum and maximum values. See 'kmin' and 'kmax'
%       above. Must be an integer greater than zero.
%   'prc_thresh' (default = 95)
%       Only used when automatically computing the bandwidths. Using 'kmin'
%       and 'kmax', the k-nn distances within the data are computed. The
%       minimum value of the bandwidths 'h' is selected to be the percentile
%       given by 'prc_thresh' of all 'kmin'-nn distances. The maximum value
%       of 'h' is chosen analogously using 'kmax'. Must be a number between
%       0 and 100.
%   'Distance' (default = 'Euclidean')
%       The distance type used for computing the KDEs. Options include the
%       standard options for functions such as knnsearch, pdist, pdist2,
%       etc.
%   'CVX' 
%       If this option is included, then the optimal weights for the
%       ensemble estimator are computing using the CVX package. CVX must be
%       installed and set up before this option can be used. The default is
%       to not use CVX.
%   'Kernel' (default = 'RBF')
%       A string or function handle specifying which kernel function to use 
%       to compute the KDEs. The value can be one of the following:
%       'gaussian'  - A Gaussian kernel with the form exp(-.5*x.^2)
%       'rbf'       - Same as a Gaussian kernel
%       'uniform'   - The uniform kernel x<=1;
%       'custom'    - The next argument must be a function that takes a
%                   single input which consists of the distances divided by 
%                   the bandwidth. In other words, the kernel function must
%                   satisfy the definition of a radial kernel.
%   'type' (default = 'KL')
%       A string or function handle specifying the type of divergence
%       functional. The value can be one of the following:
%       'KL'        - The Kullback-Leibler divergence.
%       'Renyi'     - The Renyi-divergence integral. Note that this does
%                   not take the log of the integral and divide by alpha-1.
%       'HP'        - This is the Henze-Penrose divergence. See K. Moon et
%                   al., "Meta learning of bounds on the Bayes classifier
%                   error", IEEE Signal Processing Workshop, Aug. 2015.
%       'Dp'        - Same as the HP divergence
%       'custom'    - The next argument must be a function that takes a
%                   single input that is intended for the likelihood ratio
%                   in the f-divergence functional.
%   'alpha' (default = 0.5)
%       A numeric value between 0 and 1 that gives the exponent value for
%       the Renyi divergence functional. Only used if the 'Renyi' option is
%       selected. 
%   'quiet' 
%       Adding this option suppresses the output commentary.

% Written by Kevin Moon, January 2023
% 
% Relevant Citation: K. Moon, K. Sricharan, K. Greenewald, A.O. Hero III, 
% "Ensemble Estimation of Information Divergence," Entropy (Special Issue 
% on Information Theory in Machine Learning and Data Science), vol. 20, 
% no. 8, pp. 560, July 2018. 


% Set up default parameters
eta=1;
est=1;
delta=.2;
hcustom=0;
kmax_custom=0;
kmin_custom=0;
prc_thresh=95;
L=40;
dist_type='euclidean';
cvx_w=0;
kernel='RBF';
div_type='KL';
a=0.5;
quiet=0;

% Get input parameters
for i=1:length(varargin)
%     Eta parameter
    if(strcmpi(varargin{i},'eta'))
        eta=varargin{i+1};
    end
%     Which estimator to use.
    if(strcmpi(varargin{i},'est'))
        est=varargin{i+1};
        if est~=1&&est~=2
            est=1;
            disp 'Invalid choice of estimator. Using the default value.'
        end
    end
%     Delta parameter for the modified estimator
    if(strcmpi(varargin{i},'delta'))
        delta=varargin{i+1};
        if delta<=0
            delta=.2;
            disp 'Invalid choice of delta. Using the default value.'
        end
    end
%     User-defined choice of bandwidths h
    if(strcmpi(varargin{i},'h'))
        h=varargin{i+1};
        hcustom=1;
    end
%     Parameters used to automatically select h
    if(strcmpi(varargin{i},'kmax'))
        kmax=varargin{i+1};
        kmax_custom=1;
    end
    if(strcmpi(varargin{i},'kmin'))
        kmin=varargin{i+1};
        kmin_custom=1;
    end
    if(strcmp(varargin{i},'L'))
        L=varargin{i+1};
    end
    if(strcmpi(varargin{i},'prc_thresh'))
        prc_thresh=varargin{i+1};
    end
    if(strcmpi(varargin{i},'distance'))
        dist_type=varargin{i+1};
    end
    
%     Determine if CVX should be used to solve for the weights.
    if(strcmpi(varargin{i},'cvx'))
        cvx_w=1;
    end
    
%     Choice of kernel function
    if(strcmpi(varargin{i},'kernel'))
        kernel=varargin{i+1};
        if(strcmpi(kernel,'custom'))
            kernel_func=varargin{i+2};
        end
    end
    
%     Choice of divergence type
    if(strcmpi(varargin{i},'type'))
        div_type=varargin{i+1};
        if(strcmpi(div_type,'custom'))
            div_func=varargin{i+2};
        end
        
    end
    
%     Alpha parameter for Renyi divergence
    if(strcmpi(varargin{i},'alpha'))
        a=varargin{i+1};
    end
    if(strcmpi(varargin{i},'quiet'))
        quiet=1;
    end
end

% Obtain the kernel function
switch lower(kernel)
    case 'rbf'
        kernel_func=@(x) exp(-.5*x.^2);
    case 'gaussian'
        kernel_func=@(x) exp(-.5*x.^2);
    case 'uniform'
        kernel_func=@(x) x<=1;
    case 'custom'
        disp 'Using a custom kernel function'
    otherwise
        disp 'Invalid kernel name. Using a Gaussian kernel instead.'
        kernel_func=@(x) exp(-.5*x.^2);
end

% Get sizes of the data
[N,dx]=size(X);
[M,dy]=size(Y);

if(dx~=dy)
    error('Dimensions of X and Y do not match')
end
d=dx;

% Choose the minimum value for the sample size.
N=min(N,M);

% Obtain the desired divergence functional
switch lower(div_type)
    case 'kl'
        div_func=@(x) -log(x);
    case 'renyi'
        div_func=@(x) x.^a;
    case 'hp'
        qx=N/(N+M);
        qy=1-qx;
        div_func=@(x) 1-4*qx*qy*x./(qy*x+qx);
    case 'dp'
%         This is the same as the HP divergence given above.
        qx=N/(N+M);
        qy=1-qx;
        div_func=@(x) 1-4*qx*qy*x./(qy*x+qx);
    case 'custom'
        disp 'Using a custom divergence functional'
    otherwise
        disp 'Invalid divergence functional. Computing the KL divergence'
        div_func=@(x) -log(x);
end




% Automatically select bandwidths h if not given

if ~hcustom
%     Use knn to select the minimum and maximum values
    kl=2;
    km=3;
    if ~kmin_custom
        if est==1
            kmin=round(kl*sqrt(N));
        else
            kmin=round(kl*N^(1/(d+delta)));
        end
    end
    if ~kmax_custom
        if est==1
            kmax=round(km*sqrt(N));
        else
           kmax=round(km*N^(1/(d+delta)));
        end
    end

    if kmax<=kmin
        disp 'kmax <= kmin. Setting kmax = kmin + 1.'
        kmax=kmin+1;
    end
    if ~quiet
        disp 'Automatically selecting h via k-nearest neighbors.'
        fprintf('kmin: %i \nkmax: %i \n',kmin,kmax)
    end
    
    [~,kdists]=knnsearch([X; Y],[X; Y],'k',kmax+1,'Distance',dist_type);
    if ~quiet
        fprintf('Choosing the %ith percentile distances among the data. \n',prc_thresh)
    end
    hmin=prctile(kdists(:,kmin+1),prc_thresh);
    hmax=prctile(kdists(:,kmax+1),prc_thresh);
    h=linspace(hmin,hmax,L);
else
    L=length(h);
    
end

% Compute the optimal weights

if ~quiet
    disp 'Computing the optimal weights'
end
if cvx_w
    w=calculateWeightsKDE_cvx(N,d,h,'eta',eta,'est',est,'delta',delta);
else
    w=calculateWeightsKDE(N,d,h,'eta',eta,'est',est,'delta',delta);
end

% Compute pairwise distances

if ~quiet
    disp 'Computing the pairwise distances' 
end
Dxx=squareform(pdist(X,dist_type));
Dxy=pdist2(Y,X,dist_type);

% Estimating the plug-in estimators for each bandwidth value
if ~quiet
    disp 'Computing the plug-in estimators'
end
temp_kern=zeros(L,1);
for t=1:L
    Qx=1/(M*h(t)^d)*sum(kernel_func(Dxy/h(t)));
    Px=1/((N-1)*h(t)^d)*sum(kernel_func(Dxx/h(t)));
    Ltemp=Qx./Px;
%     Check for bad values
    Ltemp(isnan(Ltemp)|Ltemp==Inf)=[];
%     Average across samples
    temp_kern(t)=mean(div_func(Ltemp));
    
end

% Take the weighted sum of the plug-in estimators
div_est=w'*temp_kern;
if ~quiet
    disp 'Ensemble estimation complete'
end

if div_est<0
    fprintf( 'Warning! Ensemble estimator returned a value less than zero. \nThis may not be appropriate for the chosen divergence functional. \nConsider changing the bandwidth range and rerunning. \n')
end
