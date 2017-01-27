function [E,Eall] = amaR01fullObjFunc(fInit,fFix,s,ctgInd,X,rMax,fano,var0,errorType,bGPU)

% function [E,Eall] = amaR01fullObjFunc(fInit,fFix,s,ctgInd,X,rMax,fano,var0,errorType,bGPU)
% 
%   example call: rMax = 7.3; fano = 0.5; var0 = 0.23; errorType = 'MAP';
%                 for i = 1:size(f,2), [E,Eall] = amaR01fullObjFunc([],AMA.f,AMA.s,AMA.ctgInd,AMA.X,rMax,fano,var0,errorType); end
%   
% objective function for AMA (revision R01)
% 
% fInit:     filters  vector magnitude of each  filter  must equal 1  [ d x nF    ]          
% fFix:      filters  that are fixed (i.e. have already been learned) [ d x nFfix ]
% s:         stimuli. vector magnitude of each stimulus must equal 1  [ d x nStm  ]       
% ctgInd:    index of category for each stimulus                      [ d x 1     ]   
% rMax:      response maximum (on average)
% fano:      response fano factor
% var0:      baseline variance
% errorType: type help amaError.m
% bGPU:      boolean on whether to attempt to use GPU
%            1 -> try to use GPU
%            0 -> don't (default)
% %%%%%%%%%%%%
% E:         cost of objective function
% Eall:      cost of each stimulus in training set

if ~exist('bGPU','var') || isempty(bGPU) bGPU = 0; end

% FILTERS
f          = [fFix fInit];                                            % combine filters into a single filter matrix, f    [nDims x nFall]    
% FILTER RESPONSE MEAN AND STANDARD DEVIATION
r          =  stim2resp(s,f,rMax);                                    % mean response from filter weights and stimuli     [nStim x nFall]     
sigma      =  resp2sigma(r,fano,var0);                                % sigma from mean response                          [nStim x nFall]    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE POSTERIOR PROBABILITY %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bGPU == 0
        [pp,ppAll] =  AMAengine(        r,r,sigma,ctgInd            );  % posterior probability of correct X value  (pp) and across all X values (ppAll)    
elseif bGPU == 1
  try   [pp,ppAll] =  AMAengineGPU_comp(r,r,sigma,ctgInd,max(ctgInd));  % posterior probability of correct X value  (pp) and across all X values (ppAll)    
  catch [pp,ppAll] =  AMAengine(        r,r,sigma,ctgInd            );  % posterior probability of correct X value  (pp) and across all X values (ppAll)    
        disp(['amaR01fullObjFunc: WARNING! bGPU = 1, but GPU is not being utilized!!! Quit and type gpuDeviceReset() at command line!']); 
  end
end
%%%%%%%%%%%%%%%%%
% COMPUTE ERROR %
%%%%%%%%%%%%%%%%%
[E,Eall]   = amaError(errorType,X,ctgInd,pp,ppAll); 
