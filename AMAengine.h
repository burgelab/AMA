#ifndef AMA_ENGINE
#define AMA_ENGINE

#include <algorithm>
#include <cmath>
#include <iostream> 
#include <vector>

#if (_MSC_VER>= 1600)       // checks whether we're using the msvc 2010 compiler
   #include<yvals.h>        // includes the original definition of char16_t in yvals.h
   #define __STDC_UTF_16__  // tells matrix.h (in preprocessor) not to redefine char16_t
#endif

// COMPUTE AMA POSTERIOR PROBABILITIES OF CORRECT CATEGORY AND ALL CATEGORIES
//
// nStm:   number of stimuli (rows of r, sigma, ctgInd, pp, ppAll)
// q:      number of filters (cols of r, sigma                   )
// R:      noisy filter response to stim kl                 [ nStm x   q  ]
// r:      mean  filter response to stimuli                 [ nStm x   q  ]
// sigma:  filter response stddev to stimuli                [ nStm x   q  ]
// ctgInd: category indices                                 [ nStm x   1  ]
//////////////////////////////////////
// pp:     posterior probability of correct category        [ nStm x   1  ]
// ppAll:  posterior probability across all categories      [ nStm x nCtg ]

void AMAposteriors (int nStm, int q, double *R, double *r, double *sigma, double *ctgInd, double *pp, double *ppAll){ // , double *Y, double *Z){
    // GET TOTAL NUMBER OF CATEGORIES
    const int nCtg = *std::max_element (&ctgInd[0], &ctgInd[0] + nStm); // call std library to get pointer to max element in a range
    // FOR EACH CANDIDATE STIMULUS
    for (int l = 0; l < nStm; ++l){   // l = kl
        std::vector<double> Y(nCtg);  // same as 'const double Y[nCtg];' but with better compatibility with old C compilers
        double              Z = 0.0;
        // Z[l]=0.0;
        
        // FOR EACH CATEGORY
        for (int i = 0; i < nCtg; ++i){ // LOOP OVER CATEGORIES: i = ij
            Y[i] = 0.0;
            // Y[l+i*nStm] = 0.0;
            // FOR EACH STIMULUS
            for (int j = 0; j < nStm; ++j){ // LOOP OVER STIMS IN CATEGORY
                int ctg = ctgInd[j]-1;      // ctgInd stores ctg from 1 to nCtg. category must be indexed from 0 to nCtg-1  
                if (ctg == i){
                    double sigPrdInv  = 1.0f;
                    double DELTA      = 0.0f;
                    // FOR EACH FILTER
                    for (int t = 0; t < q; ++t){
                        sigPrdInv  *= (1/sigma[j+t*nStm]);
                        DELTA      += pow(( (R[l+t*nStm]-r[j+t*nStm]) / (sigma[j+t*nStm]) ),2);
                    }
                    DELTA = -0.5*DELTA;
                    double prod = sigPrdInv*exp(DELTA);
                    // NUMERATOR AND DENOMINATOR
                    Y[ctg]   += prod;
                    Z        += prod;
                }
            }
        }
        
        // INDEX OF CORRECT CATEGORY
        int ctg     = ctgInd[l]-1;     // ctgInd stores ctg from 1 to q. category must be indexed from 0 to q-1  
        // POSTERIOR PROBABILITY OF CORRECT CATEGORY
        pp[l]       = Y[ctg] / Z;     
        // POSTERIOR PROBABILITY OF ALL CATEGORIES
        for (int i = 0; i < nCtg; ++i){ 
            ppAll[l+i*nStm] = Y[i] / Z;   
        } 
        
        /*
        // POSTERIOR PROBABILITY OF CORRECT CATEGORY
        int ctg     = ctgInd[l]-1;     // ctgInd stores ctg from 1 to q. category must be indexed from 0 to q-1  
        pp[l]       = Y[l+ctg*nStm] / Z[l];     
        // POSTERIOR PROBABILITY OF ALL CATEGORIES
        for (int i = 0; i < nCtg; ++i){ 
            ppAll[l+i*nStm] = Y[l+i*nStm] / Z[l];   
        }
        */

    }
}

#endif
