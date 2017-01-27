% function [pp ppAll] = AMAengine(R,r,sigma,ctgInd);
%
%   example call: % PP OF X GIVEN EXPECTED RESPONSE FOR FULL MODEL
%                   [pp ppAll] = AMAengine(r,r,sigma,ctgInd);
%
%                 % PP OF X GIVEN NOISY    RESPONSE FOR FULL MODEL
%                   [pp ppAll] = AMAengine(R,r,sigma,ctgInd);
%
% compute the posterior probability (pp) of each stimulus' category from
% the mean filter response(s) to each stimulus in the training set, 
% 
% R:      noisy filter response(s) for each stimulus in training set   [ nStm x nF   ]  
% r:      mean  filter response(s) for each stimulus in training set   [ nStm x nF   ]  
% sigma:  SD of filter response(s) for each stimulus in training set   [ nStm x nF   ]                                     
% ctgInd: category index of each stimulus in training set              [ nStm x 1    ]       
% %%%%%%%%%%%%%%%%%%%%%%
% pp:     posterior probability of correct category, for each stimulus [ nStm x 1    ]
% ppAll:  posterior probability of all categories,   for each stimulus [ nStm x nCtg ]