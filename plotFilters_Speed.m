function plotFilters_Speed(f,fSmp,figh)

% function plotFilters_Speed(f,fSmp,figh)
%
%   example call: plotFilters_Speed(f,fSmps)
%
% plot filters for speed estimation
%
% f:          filter weights
% fSmp:       values at which filter weights are sampled
% figh:       figure handle

if ~exist('fSmp','var') || isempty(fSmp) fSmp = []; end
if ~exist('figh','var') || isempty(figh) figh = []; end

% HANDLE FIGURE
if isempty(figh); figure;
else              figure(figh);
end
set(gcf,'position',[991 973 920 502]);

%% PLOT
for i = 1:size(f,2)
    % SUBPLOT ORGANIZATION: <= 4 PLOTS PER ROW
    subplot(ceil(size(f,2)./4),4,i)
    % PLOT FILTERS
    imagesc(fSmp,fSmp,reshape(f(:,i),sqrt(size(f,1)),[])'); 
    % MAKE PRETTY
    axis xy; axis image; % caxis([-.2 .2]);
    formatFigure('Position','Time',['Filter ' num2str(i)]);
    colormap gray
end