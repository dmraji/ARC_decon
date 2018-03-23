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
                             int iter,
                             int future_sight
                             )
{

  // Setting up predictive pseudo-meshes
  if(future_sight == 1)
  {
    fine = fine + 2;
  }

  /* Notes on response space:

    Is a 1D array, laid out in progressing depth from detector;

    Is laid out in terms of distance from det as follows (in cm):


  */

  // Response cm distance vector, ordered to index in conj. with counts
  float resp_dist = {4.243,  9.487,  15.297, 21.213, 27.166, 33.136, 39.115, 45.100, 51.088,
                     9.487,  12.768, 17.493, 22.847, 28.460, 34.205, 40.025, 45.891, 51.788,
                     15.297, 17.493, 21.213, 25.807, 30.887, 36.249, 41.785, 47.434, 53.160,
                     21.213, 22.847, 25.807, 29.698, 34.205, 39.115, 44.294, 49.659, 55.154,
                     27.166, 28.460, 30.887, 34.205, 38.184, 42.638, 47.434, 52.479, 57.706,
                     33.136, 34.205, 36.249, 39.115, 42.638, 46.669, 51.088, 55.803, 60.745,
                     39.115, 40.025, 41.785, 44.294, 47.434, 51.088, 55.154, 59.548, 64.203,
                     45.100, 45.891, 47.434, 49.659, 52.479, 55.803, 59.548, 63.640, 68.015,
                     51.088, 51.788, 53.160, 55.154, 57.706, 60.745, 64.203, 68.015, 72.125
                     };

  /*  Number of layers in the response matrix

      We assume that, at max, each macro-mesh will see up to 3 isotopes

  */

  if(iso_count > 3)
  {
    isos = 3;
  }
  else
  {
    isos = iso_count;
  }
  iso_iter = isos;

  combos = 0;
  while(iso_iter > 0)
  {
    combos = combos + factorial(fine * fine) / (factorial(iso_iter) * factorial(fine * fine - iso_iter));
    isos_iter--;
  }

  float*** response_spatial_reference = new float**[combos];
  for(int i = 0; i < combos; i++)
  {
    response_spatial_reference[i] = new float*[fine];
    for(int j = 0; j < fine; i++)
    {
      response_spatial_reference[i][j] = new float[fine];
    }
  }

  // Set up gradient deducer
  for(int d = 1; d < fine - 1; d++)
  {
    degrees.push_back(d);
  }

  /*
      ~~~ BEGIN RESPONSE ASSEMBLY ~~~
  */

  // For one iso
  for(int k = 0; k < fine * fine; k++)
  {
    // "Seed" (source) mesh locations in macro-meshes
    response_spatial_reference[k][floor(k / fine)][k % fine] = resp_space[0];
    for(int deg_ind_x = 0; deg_ind_x < degrees.size(); deg_ind_x++)
    {
      for(int deg_ind_y = 0; deg_ind_y < degrees.size(); deg_ind_y++)
      {
        try
        {
          // Assembly of gradient
          response_spatial_reference[k][floor(k / fine) + degrees[deg_ind_x]][(k % fine) + degrees[deg_ind_y]] =
            resp_space[(deg_ind_x + 1) + (deg_ind_y * 9) - 1];
        }
      }
    }
  }

  // For more than one iso
  if(iso_count > 1)
  {
    // Storing the index of what the next layer of the response matrix will be
    layer_save = k;
    float over_temp[fine][fine];
    float temp_layer_a[fine][fine];
    float temp_layer_b[fine][fine];
    if(isos == 3)
    {
      float temp_layer_c[fine][fine];
    }

    // Iterating through the possibilities for multiple isotopes in space
    for(iso_ind = 2; iso_ind <= isos; iso_ind++)
    {
      for(int layer_ind_a = 0; layer_ind_a < fine * fine; layer_ind_a++)
      {
        for(int layer_ind_b = layer_ind_a + 1; layer_ind_b < fine * fine; layer_ind_b++)
        {
          if(iso_ind == 2)
          {
            for(int ty = 0; ty < fine; ty++)
            {
              for(int tx = 0; tx < fine; tx++)
              {
                temp_layer_a[ty][tx] = response_spatial_reference[layer_ind_a][ty][tx];
                temp_layer_b[ty][tx] = response_spatial_reference[layer_ind_b][ty][tx];

                over_temp[ty][tx] = temp_layer_a[ty][tx] + temp_layer_b[ty][tx];

                response_spatial_reference[k][ty][tx] = over_temp[ty][tx];

              }
            }

            // Next layer of 3D response matrix
            k = k + 1;
          }
          else
          {
            for(int layer_ind_c = layer_ind_b + 1; layer_ind_c < fine * fine; layer_ind_c++)
            {
              for(int ty = 0; ty < fine; ty++)
              {
                for(int tx = 0; tx < fine; tx++)
                {
                  temp_layer_a[ty][tx] = response_spatial_reference[layer_ind_a][ty][tx];
                  temp_layer_b[ty][tx] = response_spatial_reference[layer_ind_b][ty][tx];
                  temp_layer_c[ty][tx] = response_spatial_reference[layer_ind_c][ty][tx];

                  over_temp[ty][tx] = temp_layer_a[ty][tx] + temp_layer_b[ty][tx];

                  response_spatial_reference[k][ty][tx] = over_temp[ty][tx];

                }
              }

              // Next layer of 3D response matrix
              k = k + 1;
            }
          }
        }
      }
    }
  }

  // Deleting pseudo-meshes
  if(future_sight == 1)
  {
    // New 3D array to copy old "real" contents to
    float*** response_spatial_reference_expo = new float**[combos];
    for(i = 0; i < combos; i++)
    {
      response_spatial_reference_expo[i] = new float*[fine - 2];
      for(j = 0; j < fine - 2; i++)
      {
        response_spatial_reference_expo[i][j] = new float[fine - 2];
      }
    }

    for(int lay_i = 0; lay_i < combos; lay_i++)
    {
      for(i = 0; i < fine - 2; ++i)
      {
        for(j = 0; j < fine - 2; ++j)
        {
          response_spatial_reference_expo[lay_i][i][j] = response_spatial_reference[lay_i][i + 1][j + 1];
        }
      }
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
  /*
      ~~~ END OF RESPONSE ASSEMBLY ~~~
  */

  // Begin iterations
  for(int l = 0; l < iter; l++)
  {

  }

  if(future_sight == 0)
  {
    for(i = 0; i < fine; ++i)
    {
      for(j = 0; j < fine; ++j)
      {
        delete [] response_spatial_reference[i][j];
      }
    }
    delete [] response_spatial_reference;
  }

}

spatial_decon::~spatial_decon() {}
