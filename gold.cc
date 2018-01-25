#include <TMath.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

const char* doUnfold(float *source,
								   	 float **respMatrix,
                		 int ssizex,
										 int ssizey,
                		 int numberIterations,
                		 int numberRepetitions,
										 double boost)
{
   int i, j, k, lindex, lhx = 0, repet;
   double lda, ldb, ldc, area;
   if (ssizex <= 0 || ssizey <= 0)
      return "Wrong Parameters";
   if (ssizex < ssizey)
      return "Sizex must be greater than sizey)";
   if (numberIterations <= 0)
      return "Number of iterations must be positive";
   double *working_space =
       new double[ssizex * ssizey + 2 * ssizey * ssizey + 4 * ssizex];

/*read response matrix*/
   for (j = 0; j < ssizey && lhx != -1; j++) {
      area = 0;
      lhx = -1;
      for (i = 0; i < ssizex; i++) {
         lda = respMatrix[j][i];
         if (lda != 0) {
            lhx = i + 1;
         }
         working_space[j * ssizex + i] = lda;
         area = area + lda;
      }
      if (lhx != -1) {
         for (i = 0; i < ssizex; i++)
            working_space[j * ssizex + i] /= area;
      }
   }
   if (lhx == -1) {
      delete [] working_space;
      return ("ZERO COLUMN IN RESPONSE MATRIX");
   }

/*read source vector*/
   for (i = 0; i < ssizex; i++)
      working_space[ssizex * ssizey + 2 * ssizey * ssizey + 2 * ssizex + i] =
          source[i];

/*create matrix at*a + at*y */
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
            working_space[ssizex * ssizey + 2 * ssizey * ssizey + i] = TMath::Power(working_space[ssizex * ssizey + 2 * ssizey * ssizey + i], boost);
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
