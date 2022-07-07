============
= Overview =
============
This repository provides a MATLAB implementation for learning the optimal receptive fields 
for specific tasks using Accuracy Maximization Analysis (AMA). 

Code is provided for

	+) AMA-SGD:  an implementation employing stochastic gradient descent (SGD)
		      	as described in the paper:
   		      	Johannes Burge & Priyank Jaini (2017)
	    	      	“Accuracy Maximization Analysis for Sensory-Perceptual Tasks: 
		      	Computational Improvements, Priors, and Coding Advantages for Scaled Additive Noise”
			PLoS Computational Biology, 13(2): e1005281. doi:10.1371/journal.pcbi.1005281

	+) AMA-Gauss: an implementation employing the assumption that the class-conditional distributions are 
		      	 Gaussian distributed as described in the paper:
		      	 Priyank Jaini & Johannes Burge (2017)
		      	 “Linking normative models of natural tasks with descriptive models of neural response." 
			 Journal of Vision, 17(12):16, 1-26 

	+) AMA:       the original method as described in the paper:
		      	 Wilson S. Geisler, Jiri Najemnik, Almon D. Ing (2009)
		      	 “Optimal stimulus encoders for natural tasks”
		      	 Journal of Vision, 9(13):17, 1-16

		      	 and first code release as described in the paper:
		      	 Johannes Burge & Wilson S. Geisler (2011)
		      	 “Optimal defocus estimation in individual natural images”
	              Proceedings of the National Academy of Sciences, 108(40): 16849-16854

	Please cite the appropriate work if you use this data in your research.

The repository also contains training sets of natural stimuli for two visual tasks:

	(+)  		binocular disparity estimation as described in the paper:
			 	Johannes Burge & Wilson S. Geisler (2014)
				“Optimal disparity estimation in natural stereo-images”
				Journal of Vision, 14(2):1, 1-18

	(+) 			retinal speed estimation as described in the paper:
			 	Johannes Burge & Wilson S. Geisler (2015)
				“Optimal speed estimation in natural image movies predicts human performance”
				Nature Communications, 6:7900, doi:10.1038/ncomms8900				

	Please cite the appropriate work if you use this data in your research.
	
=======================
= System Requirements =
=======================
The code has been tested on  MATLAB 2011a and later versions.  

Troubleshooting (Matlab R2011a, R2011b, R2012a and Xcode 4.2, Xcode 4.3). Any combination of these versions of Mathworks & Apple software may prevent Matlab from finding the C++ header files. To fix the issue, go to http://www.mathworks.com/support/solutions/en/data/1-FR6LXJ/ and follow the instructions

============
= Datasets =
============
		1) Binocular Disparity Estimation: AMAdataDisparity.mat
		2) Retinal Speed Estimation: 	   AMAdataSpeed.mat

==================
= Using the code =
==================
Step 1: Save the code repository and add it to your matlab path

Step 2: 	Compile the C files. Type the following commands at the prompt: 
		
		>> mex AMAengine.cpp
		>> mex AMAengineGradientMAP.cpp
		>> mex AMAengineGradientMSE.cpp
		
		NOTE! If, with one or all of these attempts to compile, the following error occurs:
		"error: cannot convert ‘const mwSize*’ {aka ‘const long unsigned int*’} to ‘const int*’ 
		in initialization" please instead attempt to compile by typing at the prompt:
		
		>> mex -compatibleArrayDims AMAengine.cpp
		>> mex -compatibleArrayDims AMAengineGradientMAP.cpp
		>> mex -compatibleArrayDims AMAengineGradientMSE.cpp
		
Step 3: Load the desired dataset using the file loadAMAdata.m 

	   	% LOAD DISPARITY DATA
	        [s, ctgInd, X] = loadAMAdata(‘Disparity’)	

		% LOAD SPEED DATA
		[s, ctgInd, X] = loadAMAdata(‘Speed’) 


Step 4: 	Run amaR01.m. 

		% EXAMPLE CALLS:
                    
		% LEARN 2 FILTERS SEQUENTIALLY W. FULL MODEL 
		[f E minTimeSec] = amaR01('FLL','MAP',2,0,1,[],s,ctgInd,X,5.7,0.5,0.23,5);   
                    
    		% LEARN 2 FILTERS SIMULATANEOUSLY W. FULL MODEL
    		[f E minTimeSec] = amaR01('FLL','MAP',2,0,2,[],s,ctgInd,X,5.7,0.5,0.23,5); 
   
    		% LEARN 2 FILTERS SIMULATANEOUSLY W. STOCHASTIC GRADIENT DESCENT             
    		[f E minTimeSec] = amaR01('SGD','MAP',2,0,2,[],s,ctgInd,X,5.7,0.5,0.23,5,570,15,.1,.001,.01);

		% LEARN 2 FILTERS SIMULTANEOUSLY W. GAUSSIAN MODEL
		[f E minTimeSec] = amaR01('GSS','MAP',2,0,2,[],s,ctgInd,X,5.7,0.5,0.23,5); 

		Note: To avoid non-global local minima in the optimization process, best practice is to 
       	 	initialize the filter learning process several times with different random seeds

Step 5: Plot the filters using the function plotFilters.m  

		% EXAMPLE CALLS:  

		% PLOT DISPARITY FILTERS
    		plotFilters('Disparity',f,[]);

    		% PLOT SPEED FILTERS
    		plotFilters('Speed',f,[]);

=========================
= Function Descriptions =
=========================
All Matlab functions in this repository contain detailed descriptions of the input and output parameters
