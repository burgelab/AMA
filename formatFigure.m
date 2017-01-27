function formatFigure(xlbl,ylbl,tlt,bLogX,bLogY,fsLbl,fsTck,xtck,ytck,xtckSD,ytckSD)

% function formatFigure(xlbl,ylbl,tlt,bLogX,bLogY,fsLbl,fsTck,xtck,ytck,xtckSD,ytckSD)
%
%   example call: formatFigure('X','Y','Data Plot',1,1,22,18,xtck,ytck)
%
% automatically format figure with axis labels, title, log/linear scaling, and specified fontsize
%
% xlbl:     x-axis label
% ylbl:     y-axis label
% bLogX:    log scale x-axis
%           0 -> linear scale (default)
%           1 -> log    scale
% bLogX:    log scale y-axis
%           0 -> linear scale (default)
%           1 -> log    scale
% fsLbl: fontsize for labels
% fsTck:  fontsize for ticks
% xtck:    x tick values
% ytck:    y tick values
% xtckSD:  number of significant digits for xtck values
% ytckSD:  number of significant digits for ytck values

if ~exist('tlt','var')     || isempty(tlt),    tlt = [];    end
if ~exist('h','var')       || isempty(h)       h = gca;     end
if ~exist('bLogX','var')   || isempty(bLogX)   bLogX = 0;   end
if ~exist('bLogY','var')   || isempty(bLogY)   bLogY = 0;   end
if ~exist('fsLbl','var')   || isempty(fsLbl)   fsLbl = 22;  end
if ~exist('fsTck','var')   || isempty(fsTck)   fsTck = 18;  end


xlabel(xlbl,'fontsize',fsLbl);
ylabel(ylbl,'fontsize',fsLbl);
title(tlt,'fontsize',fsLbl);
if bLogX == 1
    set(h,'xscale','log');
end
if bLogY == 1
    set(h,'yscale','log');
end
set(gca,'fontsize',fsTck);
set(gca,'fontWeight','normal');
set(gca,'XColor','k');
set(gca,'YColor','k');
set(gcf,'color','w');
box on

% SET TICKS AND SINGIFICANT DIGITS
if exist('xtck','var')   & ~isempty(xtck)   set(gca,'xtick',[xtck]); end
if exist('ytck','var')   & ~isempty(ytck)   set(gca,'ytick',[ytck]); end
if exist('xtckSD','var') & ~isempty(xtckSD) set(gca,'xticklabel',num2str(get(gca,'xtick')',['%.' num2str(xtckSD) 'f'] )); end
if exist('ytckSD','var') & ~isempty(ytckSD) set(gca,'yticklabel',num2str(get(gca,'ytick')',['%.' num2str(ytckSD) 'f'] )); end

% UNBOLD TITLE DEFAULTS FROM MATLABv2015 or later
v = version('-release');
if str2num(v(1:4)) >= 2015
    set(gca,'TitleFontWeight','normal');
end