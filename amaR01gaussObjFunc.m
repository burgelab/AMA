function [E,Eall,L] = amaR01gaussObjFunc(f,fLrn,s,ctgInd,X,rMax,fano,var0,errorType)

% function [E,Eall,L] = amaR01gaussObjFunc(f,fLrn,s,ctgInd,X,rMax,fano,var0,errorType)
% 
% objective function for AMA assuming Gaussian-distributed measurement distributions
% 
% s:         stimuli. vector magnitude of each stimulus must equal 1  [nDim x nStims ]       
% f:         filters  vector magnitude of each  filter  must equal 1  [nDim x nf]          
% fLrn:      filters  that have already been learned                  [nDim x nfLrn]
% ctgInd:    index of category for each stimulus                      [nDim x 1]   
% rMax:      response maximum (on average)
% fano:      response fano factor (UNUSED IN THIS FUNCTION)
% var0:      baseline variance
% errorType: error type
%           'MAP' -> maximum aposteriori
%           'MSE' -> mean squared error 
% %%%%%%%%%%%%
% E:         objective to minimize
% Eall:      cost for each stimulus
% L:         likelihood of correct category for each stimulus

% FILTER RESPONSE MEAN AND 
fAll       = [fLrn f];                       % combine all filters into a single filter matrix   [nDims x nFall]    
nF         = size(fAll,2);
r          = stim2resp( s,fAll,rMax);        % mean response from filter weights and stimuli     [nStim x nFall]     
sigma      = resp2sigma(r,fano,var0);        % sigma from mean response                          [nStim x nFall]      
% MATCH ADDITIVE NOISE VARIANCE TO AVERAGE MULTIPLICATIVE NOISE VARIANCE
varAvg     = mean(sigma.^2);
% UNCORRELATED INTERNAL NEURAL NOISE
COVint     = diag( varAvg.*ones(1,nF) ); 

% CORRELATED INTERNAL NEURAL NOISE (debugging... if response noise corr important, add to param list )
% rho = 0.0; % IF 0 < RHO <= 1.0, INTERNAL FILTER RESPONSE NOISE IS CORRELATED 
% COVint(logical(triu(ones(nF), 1))) = corr2cov(varAvg,rho);
% COVint(logical(tril(ones(nF),-1))) = corr2cov(varAvg,rho);

% COMPUTE POSTERIOR PROBABILITY (ASSUMING THAT EXPECTED COST IS APPRX EQUAL TO COST OF EXPECTED RESPONSE )
[pp,ppAll,L,Lall] =  AMAengineGSS(r,r,COVint,ctgInd,X);  % posterior probability of correct X value  (pp)     
                                                %                   and across all X values (ppAll)                                       
% AVERAGE LIKELIHOOD (WARNING! UNTESTED)
L = mean(log(L));        

% COMPUTE ERROR
[E,Eall] = amaError(errorType,X,ctgInd,pp,ppAll); 

