#include "AMAengine.h"
#include "mex.h"
        
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // CHECK THE NUMBER OF INPUTS
    if (nrhs != 4)
        mexErrMsgTxt ("AMAengine: WARNING! wrong number of inputs");

    // GET NOISY FILTER RESPONSE, R (1st RH INPUT)
    const int *dims = mxGetDimensions (prhs[0]);   // get dimensions
    int nStm = dims[0];                            // R is     [ nStm x q ]
    int q = dims[1];
	double *R = mxGetPr (prhs[0]);                 // get pointer to elements
    
    // GET MEAN  FILTER RESPONSE, r (2nd RH INPUT)
    dims = mxGetDimensions (prhs[1]);  // actual dimensions
    int rRow  = dims[0];                          // r     is  [ nStm x q ]
    int rCol     = dims[1];
    if (rRow != nStm || rCol != q)
        mexErrMsgTxt ("R and r must be the same size");
	double *r = mxGetPr (prhs[1]);                // pointer to elements

    // GET RESPONSE STDDEV, sigma (3nd RH INPUT)
    dims = mxGetDimensions (prhs[2]);             // actual dimensions
    int sigmaRow  = dims[0];                      // sigma is  [ d x nStm ]
    int sigmaCol  = dims[1];
    if (sigmaRow != nStm || sigmaCol != q)
        mexErrMsgTxt ("r and sigma must be the same size");
	double *sigma = mxGetPr (prhs[2]);            // get pointer to elements

    // GET CATEGORY INDEX, ctgInd (4rd RH INPUT)
    dims = mxGetDimensions (prhs[3]);             // actual dimensions
    int ctgIndRow = dims[0];                      // ctgInd is [ nStm x 1 ]
    int ctgIndCol = dims[1];
    if (ctgIndRow != nStm || ctgIndCol != 1)
        mexErrMsgTxt ("ctgInd must be jx1");
	double *ctgInd = mxGetPr (prhs[3]);           //  pointer to elements

    // ALLOCATE MEMORY FOR pp OUTPUT: nStm x 1 MATRIX OF POSTERIOR PROB OF CORRECT CTG
    int output1dims[2] = { nStm, 1 };
    plhs[0] = mxCreateNumericArray (2, output1dims, mxDOUBLE_CLASS, mxREAL);
    // save address of allocated memory into the return value
	double *pp = mxGetPr (plhs[0]);
    
    // DETERMINE NUMBER OF CATEGORIES
    const int nCtg = *std::max_element (&ctgInd[0], &ctgInd[0] + nStm);
    
    ///////////////////////////////////////
    // CHECK IF ctgInd IS OF TYPE SINGLE //
    ///////////////////////////////////////
    if (mxGetClassID(prhs[3]) == mxSINGLE_CLASS){ 
        std::cout << std::endl;
        std::cout << "AMAengine:            WARNING! ctgInd is of type single instead of double" << std::endl;
        std::cout << "                      nStm=" << nStm << "; nCtg=" << nCtg << ". Function will break!!!" << std::endl;
        std::cout << "                      QUITTING!!!" << std::endl;
        return;
    }
    
    // ALLOCATE MEMORY FOR ppAll OUTPUT: nStm x nCtg MATRIX OF POSTERIOR PROB ALL CTGs
    int output2dims[2] = { nStm, nCtg };
    plhs[1] = mxCreateNumericArray (2, output2dims, mxDOUBLE_CLASS, mxREAL);
    // save address of allocated memory into the return value
	double *ppAll = mxGetPr (plhs[1]);

    // POSTERIOR PROB OF CORRECT (pp) & ALL (ppAll) CATEGORIES FOR EACH STIMULUS
    AMAposteriors (nStm, q, R, r, sigma, ctgInd, pp, ppAll);

}
