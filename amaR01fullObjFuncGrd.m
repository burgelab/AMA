function [E,fgrd,Eall] = amaR01fullObjFuncGrd(f,fFix,s,ctgInd,X,rMax,fano,var0,errorType,bGPU)
                         
% function [E,fgrd,Eall] = amaR01fullObjFuncGrd(f,fFix,s,ctgInd,X,rMax,fano,var0,errorType,bGPU)
% 
%   example call: amaR01fullObjFuncGrd(AMA.f,[],AMA.s,AMA.ctgInd,AMA.X,AMA.paramRSP.rMax,AMA.paramRSP.fano,AMA.paramRSP.var0,AMA.errorType)
%
% objective function for AMA (revision R01) w. gradient computation
% 
% f:         filters  vector magnitude of each  filter  must equal 1  [ Nd x nF    ]          
% fFix:      filters  that are fixed (i.e. have already been learned) [ Nd x nFfix ]
% s:         stimuli. vector magnitude of each stimulus must equal 1  [ Nd x nStm  ]       
% ctgInd:    index of category for each stimulus                      [ Nd x 1     ]   
% rMax:      response maximum (on average)
% fano:      response fano factor
% var0:      baseline variance
% errorType: type >> help amaError.m
% bGPU:      boolean on whether to attempt to use GPU
%            1 -> try to use GPU
%            0 -> don't (default)
% %%%%%%%%%%%%
% E:         objective to minimize
% fgrd:      gradient of cost with respect to filters
% Eall:      cost of each stimulus


% INPUT HANDLING
% f = reshape(f,[],nF);
if ~strcmp(class(s),   'double') error(['amaR01fullObjFuncGrd: WARNING! input argument ''s''    is not of type ''double''']); end
if ~strcmp(class(f),   'double') error(['amaR01fullObjFuncGrd: WARNING! input argument ''f''    is not of type ''double''']); end
if ~strcmp(class(fFix),'double') error(['amaR01fullObjFuncGrd: WARNING! input argument ''fFix'' is not of type ''double''']); end
if ~exist('bGPU','var') || isempty(bGPU) bGPU = 0; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POSTERIOR PROBABILITIES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NUMBER OF FIXED FILTERS
nFfix      = size(fFix,2);
fAll       = [fFix f];                       % combine all filters into a single filter matrix   [nDims x nFall]    
r          =  stim2resp(s,fAll,rMax);        % mean response from filter weights and stimuli     [nStim x nFall]     
sigma      =  resp2sigma(r,fano,var0);       % sigma from mean response                          [nStim x nFall]      

%%%%%%%%%%%%
% GRADIENT %
%%%%%%%%%%%%
% COMPUTE EUCLIDEAN GRADIENT & POSTERIOR PROBABILITIES
if strcmp(errorType,'MAP')
    % CHANGE INPUT PARAM LIST?;                           (r,r,sigma,ctgInd,fAll,nFfix,s,rMax,fano,length(X))
    if bGPU == 0
        [fgrd,pp,ppAll] = AMAengineGradientMAP(        r,r,sigma,ctgInd,rMax.*s ,fAll,fano,nFfix);
    elseif bGPU == 1
        try   [fgrd,pp,ppAll] = AMAengineGradientMAPGPU_comp(r,r,sigma,ctgInd,rMax.*s',fAll,fano,nFfix,length(X));
        catch [fgrd,pp,ppAll] = AMAengineGradientMAP(        r,r,sigma,ctgInd,rMax.*s ,fAll,fano,nFfix);
              disp(['amaR01fullObjFuncGrd: WARNING! bGPU = 1, but GPU code for MAP gradient is not being utilized!!! Quit and type gpuDeviceReset() at command line!']); 
        end
    end
elseif strcmp(errorType,'MSE')
    % CHANGE INPUT PARAM LIST?;                         (r,r,sigma,ctgInd,fAll,nFfix,s,rMax,fano,X)
    if bGPU == 0
        [fgrd,pp,ppAll] = AMAengineGradientMSE(        r,r,sigma,ctgInd,rMax.*s ,fAll,fano,nFfix,X);
    elseif bGPU == 1
        try   [fgrd,pp,ppAll] = AMAengineGradientMSEGPU_comp(r,r,sigma,ctgInd,rMax.*s',fAll,fano,nFfix,X);
        catch [fgrd,pp,ppAll] = AMAengineGradientMSE(        r,r,sigma,ctgInd,rMax.*s ,fAll,fano,nFfix,X);
              disp(['amaR01fullObjFuncGrd: WARNING! bGPU = 1, but GPU code for MSE gradient is not being utilized!!! Quit and type gpuDeviceReset() at command line!']); 
        end
    end
else
    disp(['amaR01fullObjFuncGrd: WARNING! unhandled errorType ' errorType ' cost function... No function for computing fgrd' ]);
end
fgrd = fgrd(:,(nFfix+1):size(fgrd,2));

% COMPUTE ERROR
[E,Eall] = amaError(errorType,X,ctgInd,pp,ppAll); 
