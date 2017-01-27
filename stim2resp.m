function r = stim2resp(s,f,rMax,fType)

% function [r] = stim2resp(s,f,rMax,fType)
%
%   example call: 
% 
% s:     stimulus (usually contrast normalized)  [  d   x nStm ]
% f:     filter weights                          [  d   x nF   ]
% rMax:  maximum response
% fType: filter type
%        space   domain filters take dot product with stimulus
%        fourier domain filters take dot product with stimulus fft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('rMax', 'var') || isempty(rMax),  rMax =      1; end
if ~exist('fType','var') || isempty(fType), fType = 'NPM'; end
tol = 1e-8;
% if sum(abs(mean(s,1))>tol)>0            error('stim2resp: WARNING! each column of s must be zero mean');           end
% if sum( abs(sqrt(sum(s.^2,1))-1)>tol)>0 error('stim2resp: WARNING! each column of s must have an L2norm of zero'); end

%%%  RESPONSE FUNCTION  %%%
r = rMax.*(s'*f);
