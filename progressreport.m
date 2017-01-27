function progressreport(i,modIter,totIter)

% function progressreport(i,modIter,totIter)
%
%   example call: for i = 1:100, progressreport(i,10,100); end
%
% print progress report every mod(iter,modIter) number of iterations
%
% i:        iteration number
% modIter:   mod number of iterations
% totIter:   total number of iterations

if mod(i,modIter) == 0
    disp([num2str(i) ' of ' num2str(totIter) '. ' num2str(100.*i./totIter,'%.1f') '% complete...']); 
end