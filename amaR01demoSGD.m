%%%%%%%%%%%%%%%%%%%%%
%%% AMA- SGD DEMO %%%
%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1 : LOAD DESIRED DATA SET %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATASET 
dataType = 'Disparity';
% dataType = 'Speed';

% LOAD DATA: FOR DETAILS ON OUTPUTS, TYPE >> help loadAMAdata.m 
[s, ctgInd, X] = loadAMAdata(dataType); 

%% %%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2 : SET PARAMETERS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIX GENERAL PARAMETER VALUES
errorType = 'MAP'; %  cost function type
nF        = 2;     %  number of filters to learn
nFfix     = 0;     %  number of filters to fix
nFset     = 2;     %  number of filters in each set
f0        = [];    %  filter weights to initialize filter learning
rMax      = 5.7;   %  maximum mean response
fano      = 0.5;   %  fano factor
var0      = 0.23;  %  baseline variance
rndSd     = 3;     %  random seed

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: AMA-SGD TO LEARN & PLOT FILTERS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET SGD-SPECIFIC PARAMETERS
amaType  = 'SGD';  % ama type
btchSz   = 570;     % batch size for SGD 
nIterMax = 15;    % number of iterations (NOTE: an iteration loops over all batches)
stpSzMax = 0.1;   % initial (i.e. max) step size 
stpSzMin = 0.001; % final   (i.e. min) step size
stpSzEta = 0.01;   % factor by which stpSz is reduced on each iteration

% LEARN FILTERS: FOR DETAILS ON INPUTS AND OUTPUTS, TYPE >> help amaR01.m 
[fSGD, ESGD, minTimeSecSGD, f0SGD, E0SGD, EbtchSGD,AMASGD] = amaR01(amaType,errorType,nF,nFfix,nFset,f0,s,ctgInd,X,rMax,fano,var0,rndSd,btchSz,nIterMax,stpSzMax,stpSzMin,stpSzEta);

% PLOT FILTERS FROM AMA-SGD: FOR DETAILS ON INPUTS TYPE >> help plotFilters.m     
plotFilters('Disparity',fSGD,[]);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 4: AMA-FLL TO LEARN & PLOT FILTERS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET FLL-SPECIFIC PARAMETERS
amaType = 'FLL'; 

% LEARN FILTERS: FOR DETAILS ON INPUTS AND OUTPUTS, TYPE >> help amaR01.m 
[fFLL, EFLL, minTimeSecFLL, f0FLL, E0FLL,    ~    ,AMAFLL] = amaR01(amaType,errorType,nF,nFfix,nFset,f0,s,ctgInd,X,rMax,fano,var0,rndSd);               

% PLOT FILTERS FROM AMA-FLL: FOR DETAILS ON INPUTS TYPE >> help plotFilters.m  
plotFilters('Disparity',fFLL,[]);
