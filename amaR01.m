function [f, E, minTimeSec, f0, E0, Ebtch,AMA] = amaR01(amaType,errorType,nF,nFfix,nFset,f0,s,ctgInd,X,rMax,fano,var0,rndSd,btchSz,nIterMax,stpSzMax,stpSzMin,stpSzEta,bGPU,bUseGrd)

% function f = amaR01(amaType,errorType,nF,nFfix,nFset,f0,s,ctgInd,X,rMax,fano,var0,rndSd,btchSz,nIterMax,stpSzMax,stpSzMin,stpSzEta,bGPU,bUseGrd)
%
%   example calls: % LEARN 2 FILTERS SEQUENTIALLY W. FULL MODEL
%                    [s ctgInd     X] = loadAMAdevData('Disparity','Big');
%                    [f E minTimeSec] = amaR01('FLL','MAP',2,0,1,[],s,ctgInd,X,5.7,1.36,0.23,[]);   
%                    
%                  % LEARN 2 FILTERS SIMULATANEOUSLY
%                    [f E minTimeSec] = amaR01('FLL','MAP',2,0,2,[],s,ctgInd,X,5.7,1.36,0.23,[]); 
%
%                  % LEARN 2 FILTERS SIMULATANEOUSLY W. GAUSSIAN MODEL
%                    [f E minTimeSec] = amaR01('GSS','MAP',2,0,2,[],s,ctgInd,X,5.7,0,1.36,[]);   
%
%                  % LEARN 2 FILTERS SIMULATANEOUSLY W. STOCHASTIC GRADIENT DESCENT             
%                    [f E minTimeSec] = amaR01('SGD','MAP',2,0,2,[],s,ctgInd,X,5.7,1.36,0.23,1,570,15,.1,.001,.01);
%
%                  % LEARN 2 FILTERS SIMULATANEOUSLY W. STOCHASTIC GRADIENT DESCENT ON THE GPU             
%                    [f E minTimeSec] = amaR01('SGD','MAP',2,0,2,[],s,ctgInd,X,5.7,1.36,0.23,1,570,15,.1,.001,.01,1);
%    
% runs all amaR01* models from a single command line interface
%
% amaType:     AMA model type to run
%              'FLL' -> full     model
%              'SGD' -> full     model with stochastic gradient descent
%              'GSS' -> gaussian model with class-conditional gaussian assumption    
% errorType:   cost function type
%              'MAP' -> maximum a posteriori estimator...   optimal given KL divergence cost function
%              'MSE' -> mean squared error cost function... optimal given squared error cost function 
% nF:          number of filters to learn
%                1   -> learn 1 filter  total
%                2   -> learn 2 filters total
%               etc.
% nFfix:       number of filters to fix
% nFset:       number of filters in each set
% f0:          filter weights to initialize learning    [  nDims   x   nF  ]
%              default: []
% s:           noisy contrast normalized stimuli        [  nDims   x  nStm ]
% ctgInd:      index of category membership             [  nStm    x   1   ] 
% X:           category values                          [    1     x length(unique(ctgInd)) ] 
% rMax:        maximum mean response
% fano:        fano factor
% var0:        baseline variance
% rndSd:       random seed 
% btchSz:      batch size
% nIterMax:    number of iterations (NOTE: an iteration loops over all batches) 
% stpSzMax:    initial (i.e. max) step size
% stpSzMin:    final   (i.e. min) step size
% stpSzEta:    factor by which stpSz is reduced on each iteration
% bGPU:        boolean on whether to attempt to use GPU
%              1 -> try to use GPU
%              0 -> don't (default)
% bUseGrd:     use gradient
%              1 -> compute gradient analytically 
%              0 -> compute gradient  numerically (default)
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% f:           AMA filters                              [ nDims x   nF  ]
% E:           error (i.e. cost)                        [   1   x nDims ]
% minTimeSec:  minimization time in seconds
% f0:          AMA filters, initial values
% E0:          error (i.e. cost) for initial filter values
% Ebtch:       error for stochastic gradient descent routine
% AMA:         AMA structure with relevant outputs etc...

% INPUT HANDLING
if ~exist('nFfix'   ,'var') || isempty(nFfix)    nFfix    = 0;        end
if ~exist('nFset'   ,'var') || isempty(nFset)    nFset    = 1;        end
if ~exist('f0'      ,'var') || isempty(f0)       f0       = [];       end
if ~exist('rndSd'   ,'var') || isempty(rndSd)    rndSd = randi(2^16); end
if ~exist('bPLOT'   ,'var') || isempty(bPLOT)    bPLOT    = 1;        end
if ~exist('btchSz'  ,'var') || isempty(btchSz)   btchSz   = [];       end
if ~exist('nIterMax','var') || isempty(nIterMax) nIterMax = [];       end
if ~exist('stpSzMax','var') || isempty(stpSzMax) stpSzMax = [];       end
if ~exist('stpSzMin','var') || isempty(stpSzMin) stpSzMin = [];       end
if ~exist('stpSzEta','var') || isempty(stpSzEta) stpSzEta = [];       end
if ~exist('bGPU','var')     || isempty(bGPU)     bGPU     =  0;       end
if ~exist('bUseGrd','var')  || isempty(bUseGrd)  bUseGrd  =  0;       end
if nFset < 1 || nFset > nF,     error(['amaR01: WARNING! nFset=' num2str(nFset) '. Must be greater than 0 and less than nF= ' num2str(nF)]); end
if size(ctgInd,1) ~= size(s,2), error(['amaR01: WARNING! size(ctgInd)=[' num2str(size(ctgInd)) '] must be column vector of size=[' num2str(size(s,2)) ' 1]']); end
tol = 1e-9;
% if  min(            mean(s,1)      ) <  -tol || max(            mean(s,1)      ) >  +tol, error(['amaR01: WARNING! the columns of s must have mean   of 0.0. Check input!']); end
% if  min(unique( sqrt(sum(s.^2,1)) )) < 1-tol || max(unique( sqrt(sum(s.^2,1)) )) > 1+tol, error(['amaR01: WARNING! the columns of s must have L2norm of 1.0. Check input!']); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEED RANDOM NUMBER GENERATOR %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng('default'); rng(rndSd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE FILTER VALUES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RANDOM FILTER WEIGHTS (IF NECESSARY)
f0 = initRandomFilters(f0,s,nF);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE MINIMIZATION ROUTINE OPTIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minimizationRoutine = 'fmincon';
opts = optimset(minimizationRoutine);
opts.Algorithm      = 'active-set';
opts.LargeScale     = 'off';
opts.UseParallel    = 'always';
opts.Display        = 'iter';
opts.MaxFunEvals    = size(s,1)*800;
opts.MaxIter        = 25;
opts.TolFun         = 1e-3;
opts.TolX           = 1e-3;
opts.TolCon         = 1e-3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USE GRADIENT WITH fmincon.m OR NOT? %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bGPU == 1
    opts.UseParallel = 'never';
end
if bUseGrd == 1
    if strcmp(amaType,'FLL') && bGPU == 1, 
        bUseGrd = 0; 
        disp(['amaR01: WARNING! bUseGrd = 1 AND bGPU = 1. Setting bUseGrd = 0 to protect GPU! Update call...']);
    end
    if strcmp(amaType,'SGD') || strcmp(amaType,'GSS')
        bUseGrd = 0;
        disp(['amaR01: WARNING! bUseGrd = 1 AND amaType=' amaType '. Setting bUseGrd = 0 b/c it is non-functional in this case! Update call...']);
    end
end
if bUseGrd == 1
    opts.GradObj    = 'on';
    opts.GradConstr = 'on';
end

% GPU MESSAGES: PRINT TO SCREEN
if     bGPU == 1, disp(['amaR01: bGPU = 1. Attempting to use GPU for computations...']);
elseif bGPU == 0, disp(['amaR01: bGPU = 0. GPU use will not be attempted...']);
end

tic
f = f0;             % INITIALIZE FILTERS
nFfixTMP = nFfix;   % TEMP VARIABLE TO CONTROL NUM FIXED FILTERS
for i = 1:nF % NUMBER OF FILTERS
    if i > nFfixTMP
        for j = 1:5, % MINIMIZE 5x TO INCREASE LIKELIHOOD OF CONVERGENCE
            % SET INDICES & PARAMETERS
            ind   = (nFfixTMP+1):min(nFfixTMP+nFset,nF); % FILTER INDICES TO BE LEARNED (NOT TO EXCEED nF)
            fFix  = f(:,1:nFfixTMP);                     % FILTERS TO BE FIXED
            fInit = f(:,ind);                            % FILTERS TO BE LEARNED, INIT VALUES FOR THIS ITERATION
            % PROGRESS REPORT
            disp(['amaR01: Learning Fs ' num2str(nFfixTMP+1) '-' num2str(max(ind)) '... (iter' num2str(j) ') ' num2str(toc,'%.1f') ' seconds elapsed...']); 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            % MINIMIZE OBJECTIVE FUNCTION % 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if strcmp(amaType,'FLL')
                % FULL AMA MODEL (WITH fmincon.m)
                if bUseGrd == 0
                [f(:,ind), fval,exitflag,output,lambda,grad] = fmincon(@(rfs)  amaR01fullObjFunc(   rfs,fFix,s,ctgInd,X,rMax,fano,var0,errorType,bGPU ),fInit,[],[],[],[],[],[],@(rfs) wConEq(rfs),opts);
                elseif bUseGrd == 1
                [f(:,ind), fval,exitflag,output,lambda,grad] = fmincon(@(rfs)  amaR01fullObjFuncGrd(rfs,fFix,s,ctgInd,X,rMax,fano,var0,errorType,bGPU ),fInit,[],[],[],[],[],[],@(rfs) wConEq(rfs),opts);
                end
                killer = 1;
            elseif strcmp(amaType,'SGD')
                % STOCHASTIC GRADIENT DESCENT
                [f(:,ind), Ebtch(:,:,i)] = amaR01sgdObjFunc(fInit,fFix,s,ctgInd,X,rMax,fano,var0,errorType,nF,nFfixTMP,nFset,btchSz,nIterMax,stpSzMax,stpSzMin,stpSzEta,bGPU);         
                break;
            elseif strcmp(amaType,'GSS')
                % GAUSSIAN DISTRIBUTIONS OVER CATEGORIES (WITH fmincon.m)
                [f(:,ind), fval,exitflag,output,lambda,grad] = fmincon(@(rfs)  amaR01gaussObjFunc(rfs,fFix,s,ctgInd,X,rMax,fano,var0,errorType),fInit,[],[],[],[],[],[],@(rfs) wConEq(rfs),opts);
                killer = 1;
            end
        end
        % PRINT MINIMIZATION TIME %%% 
        minTimeSec(i) = toc; disp(['Learned Fs ' num2str(nFfixTMP+1) '-' num2str(max(ind)) ' of ' num2str(nF) ' in ' num2str(minTimeSec(i)) ' seconds']);     disp([' ']);
        % UPDATE NUMBER OF FIXED FILTERS SPECIFIC INSTRUCTIONS
        nFfixTMP = nFfixTMP + nFset; 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE RESPONSES & RESPONSE STATS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEAN RESPONSE
r     = stim2resp(s,f,rMax);
% RESPONSE STANDARD DEVIATION
sigma = resp2sigma(r,fano,var0);

% CONDITIONAL (MEAN) RESPONSE DISTRIBUTIONS 
modelCR = 'gaussian';
try
[MU,COV] = fitCondRespDistribution(modelCR,s,f,r,ctgInd);
catch
MU = []; COV = []; disp(['amaR01: WARNING! could not fit MU and COV to data']);    
end

%%%%%%%%%%%%%%%
% PCA FILTERS %
%%%%%%%%%%%%%%%
fpca = pca(s'); 
fpca = fpca(:,1:nF);

% COMPUTE ERROR AS FUNCTION OF FILTERS
for i = 1:nF
    tic;
    % ERROR FOR  FINAL  FILTERS
    E(i)    = amaR01fullObjFunc(   f(:,1:i),[],s,ctgInd,X,rMax,fano,var0,errorType);
    % ERROR FOR INITIAL FILTERS
    E0(i)   = amaR01fullObjFunc(  f0(:,1:i),[],s,ctgInd,X,rMax,fano,var0,errorType);
    % ERROR FOR  PCA FILTERS
    Epca(i) = amaR01fullObjFunc(fpca(:,1:i),[],s,ctgInd,X,rMax,fano,var0,errorType);
    disp(['amaR01: Computed error for ' num2str(i) ' of ' num2str(nF) ' initial and final filters for ' amaType ' model: E=' num2str(E(i),3) '; Time elapsed = ' num2str(toc,'%.1f') ' sec...']);
end
if exist('Ebtch','var') && ~isempty(Ebtch)
    Ebtch = Ebtch(:,:,1:nFset:end);
else
    Ebtch = [];
end

% PARAM STRUCTS
paramRSP = struct('rMax',rMax,'fano',fano,'var0',var0);
paramSGD = struct('btchSz',btchSz,'nIterMax',nIterMax,'stpSzMin',stpSzMin,'stpSzMax',stpSzMax,'stpSzEta',stpSzEta,'Ebtch',Ebtch);

% LOADUP AMA STRUCT
AMA      = struct('errorType',errorType,'amaType',amaType,'rndSd',rndSd, ...
                  'nF',nF,'nFfix',nFfix,'nFset',nFset, ...
                  'bGPU',bGPU,'bUseGrd',bUseGrd, ...
                  'X',X,'ctgInd',ctgInd, ...
                  's',s,'f',f,'f0',f0, ...
                  'r',r,'sigma',sigma, ...
                  'paramRSP',paramRSP,'paramSGD',paramSGD, ...
                  'E',E,'E0',E0,'fpca',fpca,'Epca',Epca, ...
                  'modelCR',modelCR,'MU',MU,'COV',COV, ...
                  'minTimeSec',minTimeSec);
  
%%% CONSTRAINT FUNCTION %%%
function [c,ceq,cGrd,ceqGrd] = wConEq(f)

c = []; 
% FILTER VECTOR MAGNITUDE (i.e. NORM) MUST EQUAL 1.0
ceq = sum(f.^2)-1; 
% IF GRADIENT IS PASSED
if nargout > 2
    nF = size(f,2); 
    K = kron(eye(nF),2.*f); 
    cGrd = [];
    ceqGrd = K(:,1:(nF+1):end); % GRADIENT OF CONSTRAINT FUNCTION                      
end

%%%  RESPONSE FUNCTION  %%%
function [r] = stim2resp(s,f,rMax)
r = rMax.*(s'*f);

%%%   SIGMA   FUNCTION  %%% 
function [sigma] = resp2sigma(r,fano,var0)
sigma = sqrt( fano.*abs(r) + var0 );
