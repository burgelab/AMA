function [pp,ppAll,L,Lall] =  AMAengineGSS(R,r,LAMBDA,ctgInd,X)

% function [pp,ppAll,L,Lall] =  AMAengineGSS(R,r,LAMBDA,ctgInd,X)
%
% compute posterior probability distribution assuming Gaussian-distributed
% conditional response distributions
%
% R:      noisy filter response for each stimulus           [ nStm x  nF  ]
% r:      mean  filter response for each stimulus           [ nStm x  nF  ]
% LAMBDA: covariance of internal response noise             [  nF  x  nF  ]
%         if matrix is diagonal, noise is uncorreleated
% ctgInd: category index                                    [ nStm x  1   ]
% X:      category values                                   [  1   x nCtg ]
% %%%%%%%%%%%%%%%%%
% pp:     posterior prob of correct category, for each stim [ nStm x 1    ]
% ppAll:  posterior prob of all   categories, for each stim [ nStm x nCtg ]
% L:      likelihood     of correct category, for each stim [ nStm x nCtg ]  
% Lall:   likelihood     of all   categories, for each stim [ nStm x nCtg ]

for c = 1:length(X)
    % CATEGORY INDICES
    bInd       = ctgInd==c;
    % MEAN VECTOR      : RESPONSE MEAN DUE TO SIGNAL, FOR EACH X LEVEL
    MU(c,:)    = mean(r(bInd,:));
    % COVARIANCE MATRIX: RESPONSE COVARIANCE FOR EACH X LEVEL (Ctot = Cext + Cint) 
    COV(:,:,c) = cov( r(bInd,:)) + LAMBDA;
    % LIKELIHOOD FUNCTION VALUE FOR EACH STIMULUS AT EACH X LEVEL (ctg)
    Lall(:,c)   = mvnpdf(R,MU(c,:),COV(:,:,c));
    % LIKELIHOOD FUNCTION VALUE AT CORRECT CATEGORY
    L(bInd,1)   = Lall(bInd,c);
end
% POSTERIOR PROBABILITY DISTRIBUTION
ppAll = bsxfun(@rdivide,Lall,sum(Lall,2)); % posterior probability of correct X value  (pp)   
% POSTERIOR PROBABILITY AT CORRECT CATEGORY
pp    = bsxfun(@rdivide,L   ,sum(Lall,2)); %                   and across all X values (ppAll) 