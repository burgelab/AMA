function plotFilters(taskType,f,fSmp,figh)

% function plotFilters(taskType,f,fSmp,figh)
%
%   example call: plotFilters('Speed',f,fSmp)
%
% plot filters for a number of different tasks
%
% taskType: type of task the filters are optimized for
%           'Speed'
%           'Speed3D'
%           'Disparity'
%           'Disparity2D'
%           'FigureGround'
% f:        filter values
% fSmp:     values at which filter weights are sampled
% figh:     set figure handle


if ~exist('fSmp','var')  fSmp   = []; end
if ~exist('figh','var')  figh   = []; end

if strcmp(taskType,'Defocus')
    plotFilters_Defocus(f,fSmp,figh);
elseif strcmp(taskType,'Defocus16') || strcmp(taskType,'Defocus32')
    plotFilters_Defocus(f,fSmp,figh);
elseif strcmp(taskType,'Disparity')
    plotFilters_Disparity(f,fSmp,figh);
elseif strcmp(taskType,'Disparity1D')
    plotFilters_Disparity1D(f,fSmp,figh);
elseif strcmp(taskType,'Disparity2D')
    plotFilters_Disparity2D(f,fSmp,figh);
elseif strcmp(taskType,'FigureGround')
    plotFilters_FigureGround(f,fSmp,figh);
elseif strcmp(taskType,'Spectral')
    plotFilters_Spectral(f,fSmp,figh);
elseif strcmp(taskType,'Speed')
    plotFilters_Speed(f,fSmp,figh);
elseif strcmp(taskType,'Speed3D')
    plotFilters_Speed3D(f,fSmp,figh);
else
   error(['plotFilters: WARNING! unhandled taskType: ' num2str(taskType) '. WRITE CODE!']);
end
