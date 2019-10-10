#include "AMAengineGradientMSE.h"
#include "mex.h"
#include <set>
        
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
     // CHECK THE NUMBER OF INPUTS
    if (nrhs != 9)
        mexErrMsgTxt ("AMAengineGradientMSE: WARNING! wrong number of inputs");
    
    // GET NOISY FILTER RESPONSE, R (1st RH INPUT)  
    const mwSize *dims = mxGetDimensions (prhs[0]);   // get dimensions
    mwSize nStm = dims[0];                            // R is     [ nStm x q ]
    int q = dims[1];
	double *R = mxGetPr (prhs[0]);                 // get pointer to elements

    // GET MEAN FILTER RESPONSE, r (2nd RH INPUT)
    dims = mxGetDimensions (prhs[1]);              // actual dimensions
    mwSize rRow = dims[0];                            // r is     [ nStm x q ]
    mwSize rCol = dims[1];
    if (rRow != nStm || rCol != q)
        mexErrMsgTxt ("R and r must be the same size");
	double *r = mxGetPr (prhs[1]);                 // get pointer to elements
   
    // GET RESPONSE STDDEV, sigma (3nd RH INPUT)
    dims = mxGetDimensions (prhs[2]);              // actual dimensions
    int sigmaRow = dims[0];                        // sigma is [ nStm x q ]
    int sigmaCol = dims[1];
    if (sigmaRow != nStm || sigmaCol != q)
        mexErrMsgTxt ("r and sigma must be the same size");
	double *sigma = mxGetPr (prhs[2]);             // get pointer to elements

    // GET CATEGORY INDEX, ctgInd (4rd RH INPUT)
    dims = mxGetDimensions (prhs[3]);              // actual dimensions
    int ctgIndRow = dims[0];                       // ctgInd is [ nStm x 1 ]
    int ctgIndCol = dims[1];
    if (ctgIndRow != nStm || ctgIndCol != 1)
        mexErrMsgTxt ("ctgInd must be jx1");
	double *ctgInd = mxGetPr (prhs[3]);            // get pointer to elements
    
    // GET STIMULI, s             (5th RH INPUT)
    dims = mxGetDimensions (prhs[4]);              // actual dimensions
    int d = dims[0];                               // s        [ d x nStm ]
    int sCol = dims[1];
    if (sCol != nStm)
        mexErrMsgTxt ("s must be d x nStm");
	double *s = mxGetPr (prhs[4]);                 // get pointer to elements
    
    // GET FILTERS, f             (6th RH INPUT)
    dims = mxGetDimensions (prhs[5]);              // actual dimensions
    mwSize fRow = dims[0];                            // f        [ d x q ]
    mwSize fCol = dims[1];
    if (fRow != d || fCol != q)
        mexErrMsgTxt ("f must be d x q");
	double *f = mxGetPr (prhs[5]);                 // get pointer to elements

    // GET FANO FACTOR, fano      (7th RH INPUT)
    dims = mxGetDimensions (prhs[6]);              // actual dimensions
    int fanoRow = dims[0];                            // f        [ d x q ]
    int fanoCol = dims[1];
    if (fanoRow != 1 || fanoCol != 1)
        mexErrMsgTxt ("fano must be a scalar");
	double *fanoPr = mxGetPr (prhs[6]);            // get pointer to elements
    const double fano = fanoPr[0];
    
    // GET NUMBER OF FIXED FILTERS (8th RH INPUT)
    dims = mxGetDimensions (prhs[7]);              // actual dimensions
    int nFfixRow = dims[0];                        // nFfix    [ 1 x 1 ]
    int nFfixCol = dims[1];
    if (nFfixRow != 1 && nFfixCol != 1)
        mexErrMsgTxt ("nFfix must be a scalar");
	double *nFfixPr = mxGetPr (prhs[7]); 
    const int nFfix = int(nFfixPr[0]);
    
    // GET CTG VALUES, X            (9th RH INPUT)
    dims = mxGetDimensions (prhs[8]);         // actual dimensions
    int Xrow = dims[0];                       // X        [ 1 x nCtg ]
    int Xcol = dims[1];
    // NOT ROBUST (BREAKS UNLESS AT LEAST ONE STIM SAMPLED PER CATEGORY) 
    mwSize nCtg = *std::max_element (&ctgInd[0], &ctgInd[0] + nStm);
    
    // CHECK IF ctgInd IS OF TYPE SINGLE
    if (mxGetClassID(prhs[3]) == mxSINGLE_CLASS){ 
        std::cout << std::endl;
        std::cout << "AMAengineGradientMSE: WARNING! ctgInd is of type single instead of double" << std::endl;
        std::cout << "                      nStm=" << nStm << "; nCtg=" << nCtg << ". Function will break!!!" << std::endl;
        std::cout << "                      QUITTING!!!" << std::endl;
        return;
    }

    // MORE ROBUST: GET NUMBER OF UNIQUE ELEMENTS (CRASHES OUT)
    // TRY 1 START
    //     vector<double> ctgIndCnt(nStm);                             // construct stl vector (with bounds checking, iterators, etc.)
    //     memcpy(&ctgIndCnt[0],mxGetPr(prhs[3]),nStm*sizeof(double)); // copy data from MATLAB array into vector
    //     set<double> ctgIndUnq(ctgIndCnt.begin(),ctgIndCnt.end());      // same as before, but with iterators, bounds-checking, etc.
    //     int nCtg = int(ctgIndUnq.size());
    // TRY 1 END
    // TRY 2 START
    // 	std::set<double> ctgIndUnq(&ctgInd[0], &ctgInd[0] + nStm);
    //     int nCtg = int( ctgIndUnq.size() );    
    // TRY 2 END
    
    double *X = mxGetPr (prhs[8]);            // get pointer to elements
    if (Xrow != 1 || Xcol != nCtg)            
       mexErrMsgTxt ("X must be 1 x nCtg");
    /*
    // ALLOCATE MEMORY FOR pp OUTPUT: nStm x q MATRIX OF POSTERIOR PROB OF CORRECT CTG
    int output1dims[2] = { fRow, fCol };
    plhs[0] = mxCreateNumericArray (2, output1dims, mxDOUBLE_CLASS, mxREAL);
    // save address of allocated memory into the return value
	double *grd = mxGetPr (plhs[0]);
    
    // ALLOCATE MEMORY FOR ppAll OUTPUT: nStm x nCtg MATRIX OF POSTERIOR PROB
    int output2dims[2] = { nStm, nCtg };
    plhs[1] = mxCreateNumericArray (2, output2dims, mxDOUBLE_CLASS, mxREAL);
    // save address of allocated memory into the return value
	double *ppAll = mxGetPr (plhs[1]);
    */
    
    // ALLOCATE MEMORY FOR grd   OUTPUT: nStm x q MATRIX OF POSTERIOR PROB OF CORRECT CTG
    mwSize output1dims[2] = { fRow, fCol };
    plhs[0] = mxCreateNumericArray (2, output1dims, mxDOUBLE_CLASS, mxREAL);
    // save address of allocated memory into the return value
	double *grd = mxGetPr (plhs[0]);
    
	// ALLOCATE MEMORY FOR pp    OUTPUT: nStm x 1 MATRIX OF POSTERIOR PROB OF CORRECT CTG
    mwSize output2dims[2] = { nStm, 1 };
    plhs[1] = mxCreateNumericArray (2, output2dims, mxDOUBLE_CLASS, mxREAL);
    // save address of allocated memory into the return value
	double *pp = mxGetPr (plhs[1]);
    
    // DETERMINE NUMBER OF CATEGORIES
    // const int nCtg = *std::max_element (&ctgInd[0], &ctgInd[0] + nStm);
    
    // ALLOCATE MEMORY FOR ppAll OUTPUT: nStm x nCtg MATRIX OF POSTERIOR PROB ALL CTGs
    mwSize output3dims[2] = { nStm, nCtg };
    plhs[2] = mxCreateNumericArray (2, output3dims, mxDOUBLE_CLASS, mxREAL);
    // save address of allocated memory into the return value
	double *ppAll = mxGetPr (plhs[2]);
    
    // GRADIENT FOR AMA
    AMAgradientMSE (nStm, q, d, R, r, sigma, ctgInd, s, f, fano, nFfix, X, grd, pp, ppAll );
    
}
