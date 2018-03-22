#include <string>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "spatial_decon.hh"

using namespace std;

int factorial(int n)
{
    int i, x = 1;
    for (i = 1; i <= n; i++)
    {
        x *= i;
    }
    return x;
}

spatial_decon::spatial_decon(float *resp_space,
                             int resp_space_len,
                             float **source_space,
                             int fine,
                             int iso_count,
                             int iter
                             )
{
  /* Notes on response space:

    Is a 1D array, laid out in progressing depth from detector;

    Is laid out in terms of distance from det as follows (in cm):
    15.297, 9.487, 4.243, 4.243, 9.487, 15.297,
    17.493, 12.768, 9.487, 9.487, 12.768, 17.493,
    21.213, 17.493, 15.297, 15.297, 17.493, 21.213,
    25.807, 22.847, 21.213, 21.213, 22.847, 25.807

    again, this is a 1D array; for example, resp_space[4] = counts 4.243 cm away from det;

  */

  // Response cm distance vector, ordered to index in conj. with counts
  float resp_dist = {15.297, 9.487,  4.243,  4.243,  9.487,  15.297,
                     17.493, 12.768, 9.487,  9.487,  12.768, 17.493,
                     21.213, 17.493, 15.297, 15.297, 17.493, 21.213,
                     25.807, 22.847, 21.213, 21.213, 22.847, 25.807
                     };

  // Number of responses possible per mesh
  combos = factorial(fine*fine) / (factorial(iso_count) * factorial(fine * fine - iso_count));

  float*** response_spatial_reference = new float**[combos];
  for(int i = 0; i < combos; i++)
  {
    response_spatial_reference[i] = new float*[fine];
    for(int j = 0; j < fine; i++)
    {
      response_spatial_reference[i][j] = new float[fine];
    }
  }

  // set up gradient deducer
  for(int d = 0; d < (fine - 1) * 2 + 1; d++)
  {
    deg_temp = 0 - (fine - 1) + d;
    degrees.push_back(deg_temp);
  }

  // remove zero value
  degrees.erase (degrees.begin() + (fine - 1));

  // Assemble 3D response matrix
  for(int k = 0; k < fine * fine; k++)
  {
    response_spatial_reference[k][floor(k / fine)][k % fine] = resp_space[2];
    for(int deg_ind = 0; deg_ind < degrees.size(); deg_ind++)
    {
      try
      {
        response_spatial_reference[k][floor(k / fine) + degrees[deg_ind]][k % fine] = resp_space[2 + (6 * abs(degrees(deg_ind)))];
      }
      try
      {
        response_spatial_reference[k][floor(k / fine)][(k % fine) + degrees[deg_ind]] = resp_space[2 + (6 * abs(degrees(deg_ind)))];
      }
      try
      {
        response_spatial_reference[k][floor(k / fine) + degrees[deg_ind]][(k % fine) + degrees[deg_ind]] = resp_space[2 + (5 * abs(degrees(deg_ind)))];
      }
      try
      {
        response_spatial_reference[k][floor(k / fine) + degrees[deg_ind - 1]][(k % fine) + degrees[deg_ind]] = resp_space[2 + (5 * abs(degrees(deg_ind)))];
      }
    }


  }

  // Begin iterations
  for(int l = 0; l < iter; l++)
  {

  }

  for(i = 0; i < fine; ++i)
  {
    for(j = 0; j < fine; ++j)
    {
      delete [] response_spatial_reference[i][j];
    }
  }
  delete [] response_spatial_reference;

}

spatial_decon::~spatial_decon() {}
