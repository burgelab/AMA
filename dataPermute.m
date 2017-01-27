function [sPermute, ctgIndPermute] = dataPermute(s,ctgInd,btchSz)

% function [s, ctgInd] = dataPermute(s,ctgInd,btchSz)
% 
%   example call: [a b]= dataPermute([repmat(reshape(1:15,[3 5])),1:5,5)
% 
% resamples data without replacement
% 
% s:             stimuli             [ d    x nStm ]
% ctgInd:        category indices    [ nStm x 1    ]
% btchSz:        number of elements in each batch
%%%%%%%%%%%%%%%%%%%%
% sPermute:      permuted stimuli
% ctgIndPermute: permuted category indices

[ctgIndPermute, ind] = datasample(ctgInd,btchSz,'Replace',false);
sPermute = s(:,ind);
