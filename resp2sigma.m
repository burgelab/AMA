function [sigma] = resp2sigma(r,fano,var0)

% function [sigma] = resp2sigma(r,fano,var0)
%
%   example call:
%
% SD of response given mean response
%
% r:        mean response
% fano:     fano factor
% var0:     baseline variance

sigma = sqrt( fano.*abs(r) + var0 );