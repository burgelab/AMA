#ifndef AMA_ENGINE_GRADIENT_MSE
#define AMA_ENGINE_GRADIENT_MSE

#include <algorithm>
#include <cstring>
#include <cmath>
#include <iostream> 
#include <vector>

using namespace std;

#if (_MSC_VER>= 1600)       // checks whether we're using the msvc 2010 compiler
   #include<yvals.h>        // includes the original definition of char16_t in yvals.h
   #define __STDC_UTF_16__  // tells matrix.h (in preprocessor) not to redefine char16_t
#endif

// COMPUTE AMA POSTERIOR PROBABILITIES OF CORRECT CATEGORY AND ALL CATEGORIES FOR MSE COST FUNCTION
//
// nStm:   number of stimuli (rows of R & r)
// q:      number of filters (cols of R & r)
// d:      filter dimensionality ( d )
// R:      noisy filter response to stim kl                 [ nStm x   q  ]
// r:      mean  filter response to stimuli                 [ nStm x   q  ]
// sigma:  filter response stddev to stimuli                [ nStm x   q  ]
// ctgInd: category indices                                 [ nStm x   1  ]
// s:      stimuli                                          [  d   x nStm ]  
// f:      filters                                          [  d   x   q  ]  
// fano:   fano factor for internal noise                   [  1   x   1  ]
// nFfix:  number of fixed filter shapes                    [  1   x   1  ] 
//////////////////////////////////////
// grd:    gradient                                         [  d   x   q  ]
// pp:     destination of posteriors                        [ nStm x   1  ]
// ppAll:  posterior probability across all categories      [ nStm x nCtg ]

void AMAgradientMSE (const int nStm, const int q, const int d, double *R, double *r, double *sigma, double *ctgInd, double *s, double *f, double fano, const int nFfix, double *X, double *grd, double *pp, double *ppAll){    
    
    // GET TOTAL NUMBER OF CATEGORIES
    const int nCtg = *std::max_element (&ctgInd[0], &ctgInd[0] + nStm); // call std library to get pointer to max element in a range
    // FOR EACH CANDIDATE STIMULUS
    for (int l = 0; l < nStm; ++l){             // l = kl
        // NUMERATOR (Y) & DENOMINATOR (Z) OF POSTERIOR PROBABILITY
        std::vector<double>                               Y(nCtg);    
        double                                            Z;
        // GRADIENTS OF NUMERATOR & DENOMINATOR
        std::vector<std::vector<std::vector<double> > >   grdY(    d,std::vector<std::vector<double> >(q,std::vector<double>(nCtg,0)));
        std::vector<std::vector<double> >                 grdZ(    d,std::vector<double> ( q , 0 ));        // GRADIENT OF DENOMINATOR
        // OTHER ATOMIC VARIABLES
        std::vector< double >                             Rkl_minus_rij(q, 0);
        std::vector< double >                             VARij(q, 0);
        std::vector<std::vector<double> >                 grdDELTA(d,std::vector<double> ( q , 0 ));        // GRADIENT OF SUM IN EXPONENTIAL  (ij)
        std::vector<std::vector<double> >                 grdVAR(  d,std::vector<double> ( q , 0 ));        // GRADIENT OF VARIANCE            (ij)
        // GRADIENT OF OPTIMAL ESTIMATE
        std::vector<std::vector<double> >                 grdXopt(d,std::vector<double> ( q , 0 ));
        
        // INITIALIZE Z TO ZERO
        Z = 0.0;
        // FOR EACH CATEGORY
        for (int i = 0; i < nCtg; ++i){
            // INITIALIZE Y TO ZERO
            Y[i] = 0.0;                         // same as 'const double Y[nCtg];' but with better compatibility with old C compilers
            // FOR EACH STIMULUS
            for (int j = 0; j < nStm; ++j){ // NOTE: j = ij
                int ctg = ctgInd[j]-1;      // NOTE: ctgInd stores ctg from 1 to q. category must be indexed from 0 to q-1  
                // IF CURRENT RESPONSE r(i,j) IS FROM STIM IN CURRENT CATEGORY
                if (ctg == i){
                    double sigInvPrd = 1.0f;                                     // PRODUCT OF SD^-1 ACROSS FILTERS (ij) 
                    double DELTA     = 0.0f;                                     // SUM IN EXPONENTIAL              (ij)    

                    //////////////////////////
                    // COMPUTE USEFUL TERMS //
                    //////////////////////////
                    for (int t = 0; t < q; ++t){
                        // RESPONSE DIFFERENCE
                        Rkl_minus_rij[t] = R[l+t*nStm]-r[j+t*nStm];
                        // RESPONSE VARIANCE
                        VARij[t]         = pow(sigma[j+t*nStm],2);
                        // PRODUCT OF INVERSE SIGMAS
                        sigInvPrd  *= (1/sigma[j+t*nStm]); 
                        // EXPONENTIAL SUM: SEE EQ A2 IN SGD APPENDIX: DELTA(kl,ij)... SUM WITHIN EXPONENTIAL FOR Y(k,l) & Z(k,l)
                        DELTA      +=      pow( Rkl_minus_rij[t],2) / VARij[t]; // / (sigma[j+t*nStm]) ),2);
                    }
                    DELTA = -0.5*DELTA;
                    // EXPONENTIAL OF SUM
                    double prod = sigInvPrd * exp(DELTA);
                    //////////////////////////////////////////////////////
                    // NUMERATOR & DENOMINATOR OF POSTERIOR PROBABILITY //
                    //////////////////////////////////////////////////////
                    Y[ctg] += prod; // Y(kl)
                    Z      += prod;
                    
                    ////////////////////////////////////////////
                    // COMPUTE TERMS IN GRADIENT OF NUMERATOR //
                    ////////////////////////////////////////////
                    for (int t = nFfix; t < q; ++t){ // loop over filters
                        for (int dd = 0; dd < d; ++dd){ // loop over filter dimension
                            // EQ A10 (1ST term) IN AMA SGD APPENDIX GRD DELTA W.O. MULTIPLICATIVE NOISE TERM 
                            grdDELTA[dd][t] = -( Rkl_minus_rij[t] )*( s[dd+l*d]-s[dd+j*d] ) / VARij[t]; 
                            if (fano != 0){
                                // EQ A9b IN AMA SGD APPENDIX: GRADIENT OF THE VARIANCE TERM 
                                double rSgn; // = r[j+t*nStm]/abs(r[j+t*nStm]);
                                if (r[j+t*nStm] > 0){  rSgn =  1;}
                                else                {  rSgn = -1;}
                                grdVAR[dd][t]        = rSgn*fano*s[dd+j*d];
                                // EQ A10 (2nd term) IN AMA SGD APPENDIX GRD DELTA WITH MULTIPLICATIVE NOISE TERM 
                                grdDELTA[dd][t]     += 0.5*pow(    ( Rkl_minus_rij[t]   ) / VARij[t], 2)*grdVAR[dd][t];
                            }
                            // COMPUTE grd(Yk[i]): EQ A12 IN SGD APPENDIX... LEFT SIDE: grd(Yi[i]); RIGHT SIDE: grd(Yij[i])
                            double grdYij     = prod * ( grdDELTA[dd][t] - 0.5*grdVAR[dd][t]/VARij[t] );
                            grdY[dd][t][ctg] += grdYij;
                            grdZ[dd][t]      += grdYij;
                        }
                    } 
                }
            }
        }
        ///////////////////////////////////////////////////////////
        // COMPUTE POSTERIOR PROBABILITY AT THE CORRECT CATEGORY // 
        ///////////////////////////////////////////////////////////
        // INDEX OF CORRECT CATEGORY
        int ctg     = ctgInd[l]-1;     // ctgInd stores ctg from 1 to q. category must be indexed from 0 to q-1  
        // POSTERIOR PROBABILITY OF CORRECT CATEGORY
        pp[l]       = Y[ctg] / Z;   
        
        ///////////////////////////////////////////////////////////////
        // COMPUTE OPTIMAL ESTIMATE AND GRADIENT OF OPTIMAL ESTIMATE // 
        ///////////////////////////////////////////////////////////////
        double Xopt = 0.0f;
        for (int i = 0; i < nCtg; ++i){ 
            // POSTERIOR PROBABILITY
            ppAll[l+i*nStm] = Y[i] / Z;  
            // MEAN OF THE POSTERIOR (OPTIMAL ESTIMATE FOR MSE COST FUNCTION)
            Xopt += X[i]*ppAll[l+i*nStm];
            for (int dd = 0; dd < d; ++dd){
                for (int t = nFfix; t < q; ++t){ 
                    // EQ A18 IN AMA SGD APPENDIX: GRD[Xopt(kl)]
                    //                ( Xi   )*(   p( Xi | rkl ) )*(  grdY/Y - grdZ/Z )
                    grdXopt[dd][t] += ( X[i] )*( ppAll[l+i*nStm] )*( (grdY[dd][t][i])/Y[i] - (grdZ[dd][t])/Z); 
                }
            }
        }
        //////////////////////
        // COMPUTE GRADIENT // 
        //////////////////////
        for (int dd = 0; dd < d; ++dd){
            for (int t = nFfix; t < q; ++t){
                // EQ A17 IN AMA-SGD APPENDIX 3 (AVERAGED OVER ALL STIMS)
                grd[dd+t*d] += (1/double(nStm)) * 2 * ( Xopt - X[ctg] ) * grdXopt[dd][t] ; 
            }
        }
    }
}

#endif
