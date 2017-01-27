function plotFilters_Disparity(f,fSmp,figh,nCol)

% function plotFilters_Disparity(f,fSmp,figh,nCol)
%
%   example call: plotFilters_Disparity(f,fSmp,nCol)
%
% plot filters for estimating disparity
%
% f:          filter weights
% fSmp:       values at which filter weights are sampled
% figh:       figure handle

if ~exist('fSmp','var') || isempty(fSmp) fSmp = [0:(size(f,1)/2-1)]; fSmp = fSmp-max(fSmp)/2-1;  end
if ~exist('figh','var') || isempty(figh) figh = []; end

% HANDLE FIGURE
if isempty(figh); figure;
else              figure(figh);
end
set(gcf,'position',[52 581 1629 375]);

% PLOT
for j = 1:size(f,3)
    % SUBPLOT ORGANIZATION: <= 4 PLOTS PER ROW
    if ~exist('nCol','var') || isempty(nCol) nCol = 4; end 
    for i = 1:size(f,2)
        if size(f,3) > 1
        subplot(ceil(numel(f)./(nCol.*size(f,1))),nCol,j+(i-1)*nCol)
        else
        subplot(ceil(numel(f)./(nCol.*size(f,1))),nCol,i)
        end
        % EYE-SPECIFIC FILTER INDEXES
        indLE = 1:size(f,1)/2;
        indRE = indLE + max(indLE);
        % LEFT AND RIGHT EYE FILTERS
        fLE(:,i) = f(indLE,i,j);
        fRE(:,i) = f(indRE,i,j);
        hold on;
        % PLOT
        plot(fSmp.*60,fLE(:,i),'ko-' ,'markerface','k')
        plot(fSmp.*60,fRE(:,i),'ko--','markerface','w')
        axis xy; axis square; 
        xlim(max(abs(fSmp.*60)).*[-1 1]);
        ylim([-0.5 0.5]);
        formatFigure('Position (arcmin)','Weight',['RF=' num2str(i)])
        
    end
end