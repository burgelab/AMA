% function [fgrd,pp,ppAll] = AMAengineGradientMAP(R,r,sigma,ctgInd,s,f,fano,nFfix)
%
%   example call: [fgrd,pp,ppAll] = AMAengineGradientMAP(R,r,sigma,ctgInd,s,f,fano,nFfix)
%
% computes euclidean gradient of MAP cost for simultaneous learning routine
% 
% R:      noisy filter response to stim kl              [ nStm x   q  ]
% r:      mean  filter response to stimuli              [ nStm x   q  ]
% sigma:  filter response stddev to stimuli             [ nStm x   q  ]
% ctgInd: category indices                              [ nStm x   1  ]
% s:      stimuli                                       [ nDim x nStm ]  
% f:      filters                                       [ nDim x   q  ]  
% fano:   fano factor for internal noise                [  1   x   1  ]
% nFfix:  number of fixed filter shapes                 [  1   x   1  ]
% %%%%%%%%%%%%%%%%%%%%%%
% fgrd:   euclidean gradient of filters                 [ nDim x   q  ] 
% pp:     posterior probability at correct category     [ nStm x   q  ] 
% ppAll:  posterior probability at all the categories   [ nStm x nCtg ] 
