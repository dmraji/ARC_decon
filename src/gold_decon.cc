#include <string>
#include <iostream>
#include <cmath>
#include <vector>

#include "gold_decon.hh"

using namespace std;

/*	Vars to read-in;

		"source" is 1D array of spectrum read in by detector

		"respMatrix" is 2D array of recorded spectra per per each isotope

		"sizex" is number of channels per spectra in the recorded matrix

		"sizey" is number of spectra in recorded matrix

		"numberIterations" defines the number of trials that the program will use in
			order to find convergence

		"numberRepetitions" defines ...

		"boost" is effectively a boolean that tells whether to try for extra ...

*/

// initialization of vars using member initialization list
gold_decon::gold_decon(float **rM,
					 						 float *s,
					 						 int sx,
					 						 int sy,
					 						 int nI,
					 						 int nR,
					 						 double b) :
											 respMatrix(rM),
									     source(s),
									     ssizex(sx),
									     ssizey(sy),
									     numberIterations(nI),
									     numberRepetitions(nR),
									     boost(b)
					 						 {

		// Check size of matrix to read, check number of iterations
   	if (ssizex <= 0 || ssizey <= 0)
      	return "Wrong Parameters";
   	if (ssizex < ssizey)
      	return "Sizex must be greater than sizey)";
   	if (numberIterations <= 0)
      	return "Number of iterations must be positive";

		/*	initialize working_space, is a double array with certain size

				example: for sizeX = 6, sizeY = 3, we get working_space was size 78
		*/
   	double *working_space = new double[ssizex * ssizey + 2 * ssizey * ssizey + 4 * ssizex];

		/*	read respMatrix

				iterating with var "j" and terminal condition "lhx"

		*/

		// Cycle through all columns of respMatrix
   	for (j = 0; j < ssizey && lhx != -1; j++)
		{

				/*	initializing vars "area" and "lhx"

						lhx checks whether there are rows full of zeroes

						area keeps track of the sum of the values in respMatrix

				*/
      	area = 0;
      	lhx = -1;

				// Cycle through all rows of respMatrix
      	for (i = 0; i < ssizex; i++) {
         	lda = respMatrix[j][i];

					// changing the lhx value to account for the fact that the value is nonzero
         	if (lda != 0) {
            	lhx = i + 1;
         	}

					// Set the working_space index to the value from the respMatrix
         	working_space[j * ssizex + i] = lda;

					// Add the value from respMatrix to the total area
         	area = area + lda;
      	}

      	if (lhx != -1) {
         	for (i = 0; i < ssizex; i++)
					{
							// for nonzero value in at least one value in row of respMatrix, divide entire row by total row sum
            	working_space[j * ssizex + i] /= area; // working_space = working_space / area
					}
      	}
   	}

		// Fail with output message if zero row found
   	if (lhx == -1) {
      	delete [] working_space;

				// note: was formerly "zero column"
      	return ("ZERO COLUMN IN RESPONSE MATRIX; Spectra full of zero values");
   	}

		/*	Read source vector

				Cycle through all channels; this value needs to be the same for the
					respMatrix spectra and the source vector;

		*/
   	for (i = 0; i < ssizex; i++)

			/*	Take last two lengths of ssize in working_space and assign source
					values to first ssizex length;

					Example: for ssizex = 6, ssizey = 3 (working_space length = 78),
									 positions 67 thru 72 are assigned source vector values

			*/
    	working_space[ssizex * ssizey + 2 * ssizey * ssizey + 2 * ssizex + i] =
        	source[i];

			/*	create matrix at*a + at*y



			*/
   		for (i = 0; i < ssizey; i++) {
    	for (j = 0; j < ssizey; j++) {
       	lda = 0;
       	for (k = 0; k < ssizex; k++) {
          	ldb = working_space[ssizex * i + k];
          	ldc = working_space[ssizex * j + k];
          	lda = lda + ldb * ldc;
       	}
       	working_space[ssizex * ssizey + ssizey * i + j] = lda;
    	}
    	lda = 0;
    	for (k = 0; k < ssizex; k++) {
       	ldb = working_space[ssizex * i + k];
       	ldc =
           	working_space[ssizex * ssizey + 2 * ssizey * ssizey + 2 * ssizex +
                         	k];
       					lda = lda + ldb * ldc;
    		}
    	working_space[ssizex * ssizey + 2 * ssizey * ssizey + 3 * ssizex + i] =
        	lda;
   	}

			/*move vector at*y*/
   	for (i = 0; i < ssizey; i++)
    	working_space[ssizex * ssizey + 2 * ssizey * ssizey + 2 * ssizex + i] =
        	working_space[ssizex * ssizey + 2 * ssizey * ssizey + 3 * ssizex + i];

					/*create matrix at*a*at*a + vector at*a*at*y */
   		for (i = 0; i < ssizey; i++) {
    	for (j = 0; j < ssizey; j++) {
       	lda = 0;
       	for (k = 0; k < ssizey; k++) {
          	ldb = working_space[ssizex * ssizey + ssizey * i + k];
          	ldc = working_space[ssizex * ssizey + ssizey * j + k];
          	lda = lda + ldb * ldc;
       	}
       	working_space[ssizex * ssizey + ssizey * ssizey + ssizey * i + j] =
           	lda;
    	}
    	lda = 0;
    	for (k = 0; k < ssizey; k++) {
       	ldb = working_space[ssizex * ssizey + ssizey * i + k];
       	ldc =
           	working_space[ssizex * ssizey + 2 * ssizey * ssizey + 2 * ssizex +
                         	k];
       					lda = lda + ldb * ldc;
    		}
    	working_space[ssizex * ssizey + 2 * ssizey * ssizey + 3 * ssizex + i] =
        	lda;
   	}

			/*move at*a*at*y*/
   	for (i = 0; i < ssizey; i++)
      	working_space[ssizex * ssizey + 2 * ssizey * ssizey + 2 * ssizex + i] =
          	working_space[ssizex * ssizey + 2 * ssizey * ssizey + 3 * ssizex + i];

						/*initialization in resulting vector */
   	for (i = 0; i < ssizey; i++)
      	working_space[ssizex * ssizey + 2 * ssizey * ssizey + i] = 1;

        	/***START OF ITERATIONS***/
   	for (repet = 0; repet < numberRepetitions; repet++) {
    	if (repet != 0) {
       	for (i = 0; i < ssizey; i++)
          	working_space[ssizex * ssizey + 2 * ssizey * ssizey + i] = pow(working_space[ssizex * ssizey + 2 * ssizey * ssizey + i], boost);
    	}
    	for (lindex = 0; lindex < numberIterations; lindex++) {

       	for (i = 0; i < ssizey; i++) {
          	lda = 0;
          	for (j = 0; j < ssizey; j++) {
             	ldb =
                 	working_space[ssizex * ssizey + ssizey * ssizey + ssizey * i + j];
             	ldc = working_space[ssizex * ssizey + 2 * ssizey * ssizey + j];
             	lda = lda + ldb * ldc;
          	}
          	ldb =
              	working_space[ssizex * ssizey + 2 * ssizey * ssizey + 2 * ssizex + i];
          	if (lda != 0) {
             	lda = ldb / lda;
          	}

          	else
             	lda = 0;
          	ldb = working_space[ssizex * ssizey + 2 * ssizey * ssizey + i];
          	lda = lda * ldb;
          	working_space[ssizex * ssizey + 2 * ssizey * ssizey + 3 * ssizex + i] = lda;
       	}
       	for (i = 0; i < ssizey; i++)
          	working_space[ssizex * ssizey + 2 * ssizey * ssizey + i] =
              	working_space[ssizex * ssizey + 2 * ssizey * ssizey + 3 * ssizex + i];
    		}
   	}

		/*write back resulting spectrum*/
   	for (i = 0; i < ssizex; i++) {
      	if (i < ssizey)
         	source[i] = working_space[ssizex * ssizey + 2 * ssizey * ssizey + i];

      	else
         	source[i] = 0;
   	}
   	delete[]working_space;
   	return 0;

}

gold_decon::~gold_decon() {}
