function [ f0 ] = initRandomFilters( f0,s,nF )

% function [f0] = initRandomFilters( f0,s,nF )
%
%   example call: % INIT 2 RANDOM FILTERS
%                   f0 = initRandomFilters([],s,2);
%
% initializes random filters. if some filters are fixed, initializes
% sufficient number of random filters to make nF filters 
%
% f0:      filters that are input
% s:       stimulus matrix          [ d x nStm ]
% nF:      total number of filters
%%%%%%%%%%%%%%%%%
% f0:      random initial filters

% RANDOM FILTER WEIGHTS
f0 = [f0(:,1:min([nF size(f0,2)])) randn(size(s,1),nF-size(f0,2))];
% NORMALIZED TO L2 NORM OF 1.0
f0 = bsxfun(@rdivide,f0,sqrt(sum(f0.^2,1)));

