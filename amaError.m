function [E,Eall] = amaError(errorType,X,ctgInd,pp,ppAll) 

% function [E Eall] = amaError(errorType,X,ctgInd,pp,ppAll)
% 
% compute cost for different error types (cost functions) given a 
% posterior probability distribution over the categories
% 
% errorType:   cost function type
%              'MAP' -> maximum a posteriori estimator   (KL-divergence)
%              'MSE' -> mean squared error cost function (squared error)
% X:           category values                           [ 1    x nCtg ]
% ctgInd:      category index for each stimulus          [ nStm x   1  ]
% pp:          posterior probability at correct category [ nStm x   1  ]
% ppAll:       posterior probability at  each   category [ nStm x nCtg ]
%                      [pp,ppAll]=AMAengine(.);
% %%%%%%%%%%%%%%%
% E:           error, averaged across all stimuli        [ 1    x  1   ]
% Eall:        error, for each stimulus                  [ nStm x  1   ]

if     strcmp(errorType,'MAP') || strcmp(errorType,'ARE') || strcmp(errorType,'L0N') 
    Eall = -log(pp);              % (ARE) (A)verage (R)elative (E)ntropy between empirical and ideal posteriors 
    E    = mean(Eall);             % average relative entropy
elseif strcmp(errorType,'MSE') || strcmp(errorType,'L2N') || strcmp(errorType,'Est') 
    Xhat =  ppAll*X';              % expected value of posterior 
    Eall = (Xhat - X(ctgInd)').^2; % MSE between expected value of posterior and true value  
    E    = mean(Eall);              
elseif strcmp(errorType,'MED') || strcmp(errorType,'L1N')
    Xhat  = zeros(size(ppAll,1),1);
    for i = 1:size(ppAll,1)        % INTERPOLATE MEDIAN
    Xhat(i,1) = interp1(cumsum(ppAll(i,:)),X,0.5,'linear',min(X));
    end
    Eall = abs(Xhat - X(ctgInd)');
    E    = mean(Eall);
    disp(['amaError: WARNING! untested code for errorType: ' errorType]);
elseif strcmp(errorType,'circMSE') 
    Xhat = circ_mean(X', ppAll')';        % expected value of response 
    Eall = circ_dist(X(ctgInd)',Xhat).^2; % circular distance squared between expected value of response and true value  
    E    = mean(Eall);       % 
else
    error(['amaError: WARNING! unhandled error type: ' errorType '. WRITE CODE!']);
end