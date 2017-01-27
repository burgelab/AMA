function [fgrd,pp,ppAll] = gradSGD(errorType,r,sigma,ctgInd,s,nFfix,nFset,f,X,rMax,fano,bGPU)

% function [fgrd,pp,ppAll] = gradSGD(errorType,r,sigma,ctgInd,s,nFfix,nFset,f,X,rMax,fano,bGPU)
%
%   example call: [fgrd,pp,ppAll] = gradSGD('MAP',rBtch,sigmaBtch,ctgIndBtch,sBtch,nFfix,nFset,f,X,fano,bGPU)
%
% evaluates the euclidean gradient for current filters and current stim batch 
% AND projects the euclidean gradient on the tangent plane of hypersphere
% 
% NOTE! the gradient code expects batches of r, sigma, ctgInd, and s...
%       otherwise there is little benefit to using this routine
%
% errorType:  cost function type
%             'MAP' -> maximum a posteriori estimator
%             'MSE' -> mean squared error cost function
% r:          response for the stimuli in the current batch
% sigma:      sigma for the current batch under consideration
% ctgInd:     category index of the stimuli in the current batch
% s:          stimuli in the current batch under consideration
% nFfix:      number of fixed filters
% nFset:      number of filters per set
% f :         current set of filters                    [ d  x  nF  ]
% X:          category values                           [ 1  x nCtg ]
% fano:       fano factor
% bGPU:       boolean on whether to attempt to use GPU
%             1 -> try to use GPU
%             0 -> don't (default)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fgrd:       gradient on tangent plane of hypersphere
%             (i.e. the dimensionality is [ d x nF ] where nF = size(f,2)
% pp:         posterior probability at correct category for all stimuli
% ppAll:      posterior probability distributions       for all stimuli

if ~strcmp(class(r),'double') error(['gradSGD: WARNING! input argument ''r'' is not of type ''double''']); end
if ~strcmp(class(f),'double') error(['gradSGD: WARNING! input argument ''f'' is not of type ''double''']); end
if ~exist('bGPU','var') || isempty(bGPU) bGPU =  0; end

% FIND INDICES OF FILTERS THAT ARE NOT FIXED
ind = (nFfix+1):min(nFfix+nFset,size(f,2));

if strcmp(errorType, 'MAP') || strcmp(errorType, 'L0N')
    % RETURN GRADIENT OF SIZE [ d x nF ] ... where nF = size(f,2);
    % R = r + sigma.*randn(size(r));
    if bGPU == 0
        [fgrd_euclid,pp,ppAll] = AMAengineGradientMAP(        r,r,sigma,ctgInd,rMax.*s ,f,fano,nFfix);
    elseif bGPU == 1
        try   [fgrd_euclid,pp,ppAll] = AMAengineGradientMAPGPU_comp(r,r,sigma,ctgInd,rMax*s',f,fano,nFfix,length(X));  if sum(isnan(fgrd_euclid(:))) ==0, disp('Yippee!!! MAPgradientGPU works!'); else disp('Yikes!!! MAPgradientGPU don''t work!');  end
        catch [fgrd_euclid,pp,ppAll] = AMAengineGradientMAP(        r,r,sigma,ctgInd,rMax*s ,f,fano,nFfix); 
              disp(['gradSGD: WARNING! bGPU = 1, but GPU is not being used w. AMAengineGradientMAP!!! Quit, type gpuDeviceReset(), and retry!']);
        end     
    end
    killer = 1;
elseif strcmp(errorType, 'MSE') || strcmp(errorType, 'L2N')
    % RETURN GRADIENT OF SIZE [ d x nF ] ... where nF = size(f,2);
    if length(X)~=length(unique(ctgInd)), disp(['gradSGD: WARNING! length(unique(ctgInd))=' num2str(length(unique(ctgInd))) ' does not match length(X)=' num2str(length(X))]); 
                                          disp(['         Increase the batch size so sampled batches have at least one sample per ctg ']); 
                                          disp(['         Code likely to error because of non-robust code in mex file.']); 
                                          disp(['         See notes in AMAengineGradientMSE.cpp when ready to fix!']); 
    end
    % EUCLIDEAN GRADIENT 
    % R = r + sigma.*randn(size(r));
    if bGPU == 0
        [fgrd_euclid,pp,ppAll] = AMAengineGradientMSE(        r,r,sigma,ctgInd,rMax.*s ,f,fano,nFfix,X(unique(ctgInd))); 
    elseif bGPU == 1
        try   [fgrd_euclid,pp,ppAll] = AMAengineGradientMSEGPU_comp(r,r,sigma,ctgInd,rMax.*s',f,fano,nFfix,X(unique(ctgInd))); if sum(isnan(fgrd_euclid(:))) ==0, disp('Yippee!!! MSEgradientGPU works!');  end
        catch [fgrd_euclid,pp,ppAll] = AMAengineGradientMSE(        r,r,sigma,ctgInd,rMax.*s ,f,fano,nFfix,X(unique(ctgInd))); 
              disp(['gradSGD: WARNING! bGPU = 1, but GPU is not being used w. AMAengineGradientMSE!!! Quit, type gpuDeviceReset(), and retry!']);
        end
    end
end
% PROJECT GRADIENT ON TANGENT PLANE OF HYPERSPHERE AND NORMALIZE... see ../ama/Notes/SGD_Overview.pdf
for i = ind
    % PROJECT ONTO TANGENT PLANE
    fgrd(:,i) = fgrd_euclid(:,i) - (f(:,i)*f(:,i)')*fgrd_euclid(:,i);
    % NORMALIZE
    fgrd(:,i) = bsxfun(@rdivide,fgrd(:,i),sqrt(sum(fgrd(:,i).^2)));
end   