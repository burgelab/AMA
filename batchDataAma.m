function [sBtch, ctgIndBtch, rBtch, sigmaBtch] = batchDataAma(btchSz,b,f,sBtchAll,ctgIndBtchAll,rMax,fano,var0)

% function [ind, sBtch, ctgIndBtch, rBtch, sigmaBtch] = batchDataAma(btchSz,b,f,sBtchAll,ctgIndBtchAll,rMax,fano,var0)
%
%   example call : [sBtch, ctgIndBtch, rBtch, sigmaBtch] = batchDataAma(btchSz,b,f,sBtchAll,ctgIndBtchAll,rMax,fano,var0)
%
% evaluates parameters from the data required for gradient descent for the current batch of permuted data 
%
% btchSz:        number of stimuli per batch
% b:             iteration over different batches
% f :            current set of filters
% sBtchAll:      stimuli for the permuted dataset(all)
% ctgIndBtchAll: category indices of the stimuli in the permuted dataset(all)
% rMax:          response maximum (on average)
% fano:          response fano factor
% var0:          baseline variance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sBtch:         stimuli in the current batch under consideration
% ctgIndBtch:    category index of the stimuli in the current batch 
% rBtch:         response for the stimuli in the current batch
% sigmaBtch:     sigma for the current batch under consideration
%
% INDICES OF CURRENT BATCH
ind = [1:btchSz] + (b-1).*btchSz;
% STIMULI AND CTG INDEX FOR BATCH
sBtch      = sBtchAll(:,ind);
ctgIndBtch = ctgIndBtchAll(ind);
% COMPUTE RESPONSE
rBtch      = stim2resp(sBtch,f,rMax);
% COMPUTE RESPONSE SD
sigmaBtch  = resp2sigma(rBtch,fano,var0);