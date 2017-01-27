function [fInitNew, Ebtch] = updateSGD(errorType, fFix,fInit, stpSz, fgrd, ind,sBtch,ctgIndBtch,X,rMax,fano,var0, Ebtch, b, j)

% [fInitNew, Ebtch] = updateSGD(errorType, fFix,fInit, stpSz, fgrd, ind,sBtch,ctgIndBtch,X,rMax,fano,var0, Ebtch, b, j)
%
%   example call : [fInit, f, Ebtch] = updateSGD(orderType, errorType, fInit, f, stpSz, fgrd, ind,sBtch,ctgIndBtch,X,rMax,fano,var0, Ebtch, b, j)
% 
%  updates filters by taking a step in direction of gradient
%
% errorType:    cost function type
%               'MAP' -> maximum a posteriori estimator
%               'MSE' -> mean squared error cost function
% fFix:         fixed filters for the iteration
% fInit:        initial filters for the iteration
% stpSz:          step size for the current update step
% fgrd:         gradient in tangent plane to hypersphere for descent
% ind:          current filter number  being learnt
% rBch:         response for the stimuli in the current batch
% sBtch:        stimuli in the current batch under consideration
% ctgIndBtch:   category index of the stimuli in the current batch
% X:            category values                        [ 1 x nCtg ]
% rMax:         response maximum (on average)
% fano:         response fano factor
% var0:         baseline variance
% Ebtch:            error associated with each iteration
% b:            counter for current batch in consideration
% j:            counter for iterations over the whole dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fInitNew:     updated filters (will serve as initial point for next iter)
% f:            filters after update
% Ebtch:        error associated with each iteration
%

% % COST W. CURRENT BATCH & FILTERS
% Ebtch(b,j)    = amaR01fullObjFunc(fInit,fFix,sBtch,ctgIndBtch,X,rMax,fano,var0,errorType);
                
% STEP ALONG GRADIENT W. SPECIFIED SIZE
fInitNew     = fInit - stpSz.*fgrd(:,ind);
% NORMALIZE FILTERS
fInitNew     = bsxfun(@rdivide,fInitNew,sqrt(sum(fInitNew.^2)));
% ERROR W. CURRENT BATCH AND NEW FILTERS
Ebtch(b+1,j) = amaR01fullObjFunc(fInitNew,fFix,sBtch,ctgIndBtch,X,rMax,fano,var0,errorType);
% RETURN INPUT FILTERS ONLY IF STEP REDUCES COST
if Ebtch(b+1,j) > Ebtch(b,j)
    fInitNew = fInit;
end

