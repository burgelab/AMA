function [f Ebtch] = amaR01sgdObjFunc( fInit,fFix,s,ctgInd,X,rMax,fano,var0,errorType,nF,nFfix,nFset,btchSz,nIterMax,stpSzMax,stpSzMin,stpSzEta,bGPU)

% function [f Ebtch] = amaR01sgdObjFunc( fInit,fFix,s,ctgInd,X,rMax,fano,var0,errorType,nF,nFfix,nFset,btchSz,nIterMax,stpSzMax,stpSzMin,stpSzEta,bGPU)
%
%   example call:
%
% implements stochastic gradient descent routine
%
% fInit:      initial filter values                                    [ d x <= nF ] 
% s:          stimuli. vector magnitude of each stimulus must equal 1  [ d x nStm  ]
% ctgInd:     index of category for each stimulus                      [ d x 1     ]
% rMax:       response maximum (on average)
% fano:       response fano factor
% var0:       baseline variance
% errorType:  cost function type
%             'MAP' -> maximum a posteriori estimator
%             'MSE' -> mean squared error cost function
% nF:         number of filters to learn
% nFfix:      number of fixed filters
% nFset:      number of filters in each set
% btchSz:     batch size
% nIterMax:   number of iterations (NOTE: an iteration loops over all batches) 
% stpSzMax:   initial (i.e. max) step size
% stpSzMin:   final   (i.e. min) step size
% stpSzEta:   factor by which stpSz is reduced on each batch
% bGPU:       boolean on whether to attempt to use GPU
%             1 -> try to use GPU
%             0 -> don't (default)
%%%%%%%%%%%%%%%%%%%%%
% f:          learned filters 
% Ebtch:      error for each batch

if btchSz > size(s,2), error(['amaR01sgdObjFunc: WARNING! requested btchSz=' num2str(btchSz) ' larger than nStm=' num2str(size(s,2)) '! Fix it...']); end
if ~exist('bGPU','var')     || isempty(bGPU)     bGPU     =  0;       end

% TOTAL NUMBER OF BATCHES
nBtch = floor(size(ctgInd,1)./btchSz);

% STARTING (MAXIMUM) STEP SIZE
stpSz = stpSzMax;
% INDICES OF FILTERS CURRENTLY BEING LEARNED
ind   = (nFfix+1):min(nFfix+nFset,nF);
for j = 1:nIterMax
    % PERMUTE THE DATA
    [sBtchAll,ctgIndBtchAll] = dataPermute(s,ctgInd,size(ctgInd,1));
    for b = 1:nBtch
        % SAMPLE BATCHES FROM RELEVANT PARAMETERS
        [sBtch, ctgIndBtch, rBtch, sigmaBtch] = batchDataAma(btchSz,b,[fFix fInit],sBtchAll,ctgIndBtchAll,rMax,fano,var0);
        % GRADIENT OF OBJECTIVE FUNCTION ON TANGENT PLANE OF HYPERSPHERE
        [fgrd,pp,ppAll] = gradSGD(errorType, rBtch, sigmaBtch, ctgIndBtch, sBtch,nFfix,nFset,[fInit fFix], X, rMax,fano,bGPU);
        % COST W. CURRENT BATCH & FILTERS
        Ebtch(b,j)      = amaError(errorType,X,ctgIndBtch,pp,ppAll); 
        % UPDATE FILTERS 
        [fInit,Ebtch]   = updateSGD(errorType,fFix,fInit,stpSz, fgrd, ind,sBtch,ctgIndBtch,X,rMax,fano,var0, Ebtch, b, j);     
        % UPDATE STEP-SIZE FOR NEXT ITERATION
        stpSz           = max(stpSzMin,stpSz - stpSzEta.*stpSz);
        % PLOT FILTER EVOLUTION
        bPLOT = 1;
        if bPLOT  && b == 1
        plotSGDcostNfilters(fInit,nBtch,Ebtch,stpSz,b,j,nIterMax,sBtch,ctgIndBtch,rMax,fano,var0);
        end
    end
    % PROGRESS REPORT
	progressreport(j,2,nIterMax)
end
% SAVE FILTERS
f = fInit;

function bQuit = plotSGDcostNfilters(fInit,nBtch,Ebtch,stpSz,b,j,nIterMax,s,ctgInd,rMax,fano,var0)

% function bQuit = plotSGDcostNfilters(fInit,nBtch,Ebtch,stpSz,b,j,nIterMax,s,ctgInd,rMax,fano,var0)
%
% plot sgd routine data

% PLOT FILTER
clf; % PLOT ERROR
global Eavg
if b==1 && j == 1 Eavg = []; end
if b==nBtch       Eavg = mean(Ebtch,1); end

subplot(1,3,1); hold on;
indPlt = 1:(nBtch*(j-1) + b);
plot(Ebtch(indPlt),'k');
if j >= 1 && exist('Eavg','var') && ~isempty(Eavg)
    
    indPlt2 = ceil(nBtch/2):nBtch:(ceil(nBtch/2)+(nBtch*(j)));
    indPlt2 = intersect(indPlt,indPlt2);
    [ b j];
    plot( indPlt2(1:length(Eavg)), Eavg(1:end),'ro-','linewidth',3);
    
end 
axis square;
xlim([0 nIterMax*nBtch+1]); ylim([floor(min(Ebtch(indPlt))) ceil(max(Ebtch(indPlt)))]);
formatFigure('Number of Batches','Cost',['Cost=' num2str(Ebtch(b+1,j),'%.4f')]);

subplot(1,3,2); hold on;
plot(fInit); plot(fInit,'--'); axis square
formatFigure('Filter Position','Filter Values',['stpSz=' num2str(stpSz,'%.4f')]);
xlim([0 size(fInit,1)+1]); ylim([-.5 .5]);


try 
    subplot(1,3,3); 
    %%
    hold on; cla;
    r     = stim2resp(s,fInit,rMax);
    sigma = resp2sigma(r,fano,var0);
    ctgIndUnq = unique(ctgInd);
   
    colors  = repmat(get(gca,'colorOrder'),[5 1]);
    axis(rMax*[-1 1 -1 1]); 
    formatFigure('R1','R2');
    numCtg = length(unique(ctgInd));
    axis square;
    for c = 1:length(ctgIndUnq)
        ind  = find(ctgInd == ctgIndUnq(c)); 
        plot(r(ind,1),r(ind,2),'.','color',colors(c,:))
        for i = sort(randsample(1:length(ind),1 ))
        plotEllipse(r(ind(i),:),diag(sigma(ind(i),:).^2),5,[],0.5,colors(c,:));
        end
    end
    
end
pause(.00001);