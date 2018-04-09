#include <string>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <ctime>

#include "spatial_decon.hh"

using namespace std;

// int factorial(int n)
// {
//     int i, x = 1;
//     for (i = 1; i <= n; i++)
//     {
//         x *= i;
//         std::cout << x << '\n';
//     }
//     return x;
// }

spatial_decon::spatial_decon(float *resp_space,
                             int resp_space_len,
                             float **source_space,
                             int survsizey,
                             int survsizex,
                             int fine,
                             int iso_count,
                             int iter,
                             int future_sight
                             )
{

  std::cout << "Begin spatial decon." << '\n';

  // Validate input parameters
  if(fine != 3 && fine != 5 && fine != 7)
  {
    std::cout << "Invalid fineness specification. Terminating." << '\n';
    exit(1);
  }

  // Setting up predictive pseudo-meshes
  if(future_sight == 1)
  {
    fine = fine + 2;
  }

  conf_levels = {0.90, 0.925, 0.94, 0.96, 0.97};

  /* Notes on response space:

    Is a 1D array, laid out in progressing depth from detector;

    Is laid out in terms of distance from det as follows (in cm):


  */

  // Response cm distance vector, ordered to index in conj. with counts
  float resp_dist[81] = {4.243,  9.487,  15.297, 21.213, 27.166, 33.136, 39.115, 45.100, 51.088,
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

  // For sim stuff
  iso_count = 3;

  std::cout << "space 69" << '\n';

  if(iso_count > 3)
  {
    isos = 3;
  }
  else
  {
    isos = iso_count;
  }
  iso_iter = isos;

  // Combinations without large numbers
  combos = 0;

  while(iso_iter > 0)
  {
    it_prod = 1;
    combo_piece = 0;
    for(int it = 0; it < iso_iter; it++)
    {
      if(it == 0)
      {
        combo_piece = fine * fine;
      }
      else
      {
        combo_piece = combo_piece * (fine * fine - it);
        it_prod = it_prod * (it + 1);
      }

    }
    std::cout << "combo prog: " << combos << '\n';
    combos = combos + combo_piece / it_prod;
    iso_iter--;
  }



  // combos =+ factorial(fine * fine) / (factorial(iso_iter) * factorial(fine * fine - iso_iter));

  std::cout << "combos: " << combos << '\n';

  float*** response_spatial_reference = new float**[combos];
  for(int i = 0; i < combos; i++)
  {
    response_spatial_reference[i] = new float*[fine];
    for(int j = 0; j < fine; j++)
    {
      response_spatial_reference[i][j] = new float[fine];
    }
  }

  std::cout << "space 100" << '\n';

  // Set up gradient deducer
  for(int d = 0; d < 9; d++)
  {
    degrees.push_back(d);
  }

  /*
      ~~~ BEGIN RESPONSE ASSEMBLY ~~~
  */

  std::cout << "space 108" << '\n';

  int k_;

  // For one iso
  for(int k = 0; k < fine * fine; k++)
  {

    // "Seed" (source) mesh locations in macro-meshes
    response_spatial_reference[k][int (floor(k / fine))][int (k % fine)] = resp_space[0];

    for(int deg_ind_x = 0; deg_ind_x < degrees.size(); deg_ind_x++)
    {
      for(int deg_ind_y = 0; deg_ind_y < degrees.size(); deg_ind_y++)
      {
        // Assembly of gradient
        if(deg_ind_y == 0 && deg_ind_x == 0)
        {
          // std::cout << "seed" << '\n';
        }
        else
        {
          // std::cout << "k: " << k << '\n';
          if(floor(k / fine) + degrees[deg_ind_x] < fine && (k % fine) + degrees[deg_ind_y] < fine)
          {
            // std::cout << "xplus yplus" << '\n';
            response_spatial_reference[k][int (floor(k / fine) + degrees[deg_ind_x])][int ((k % fine) + degrees[deg_ind_y])] =
              resp_space[(deg_ind_x + 1) + (deg_ind_y * 9) - 1];
              // std::cout << (floor(k / fine) + degrees[deg_ind_x]) << ", " << ((k % fine) + degrees[deg_ind_y]) << '\n';
              // std::cout << resp_space[(deg_ind_x + 1) + (deg_ind_y * 9) - 1] << '\n';
          }
          if(deg_ind_x > 0)
          {
            if(floor(k / fine) - degrees[deg_ind_x] >= 0 && (k % fine) + degrees[deg_ind_y] < fine)
            {
              // std::cout << "xminus yplus" << '\n';
              response_spatial_reference[k][int (floor(k / fine) - degrees[deg_ind_x])][int ((k % fine) + degrees[deg_ind_y])] =
                resp_space[(deg_ind_x + 1) + (deg_ind_y * 9) - 1];
                // std::cout << (floor(k / fine) - degrees[deg_ind_x]) << ", " << ((k % fine) + degrees[deg_ind_y]) << '\n';
                // std::cout << resp_space[(deg_ind_x + 1) + (deg_ind_y * 9) - 1] << '\n';
            }
          }
          if(deg_ind_y > 0)
          {
            if(floor(k / fine) + degrees[deg_ind_x] < fine && (k % fine) - degrees[deg_ind_y] >= 0)
            {
              // std::cout << "xplus yminus" << '\n';
              response_spatial_reference[k][int (floor(k / fine) + degrees[deg_ind_x])][int ((k % fine) - degrees[deg_ind_y])] =
                resp_space[(deg_ind_x + 1) + (deg_ind_y * 9) - 1];
                // std::cout << (floor(k / fine) + degrees[deg_ind_x]) << ", " << ((k % fine) - degrees[deg_ind_y]) << '\n';
                // std::cout << resp_space[(deg_ind_x + 1) + (deg_ind_y * 9) - 1] << '\n';
            }
          }
          if(deg_ind_y > 0 && deg_ind_x > 0)
          {
            if(floor(k / fine) - degrees[deg_ind_x] >= 0 && (k % fine) - degrees[deg_ind_y] >= 0)
            {
              // std::cout << "xminus yminus" << '\n';
              response_spatial_reference[k][int (floor(k / fine) - degrees[deg_ind_x])][int ((k % fine) - degrees[deg_ind_y])] =
                resp_space[(deg_ind_x + 1) + (deg_ind_y * 9) - 1];
                // std::cout << (floor(k / fine) - degrees[deg_ind_x]) << ", " << ((k % fine) - degrees[deg_ind_y]) << '\n';
                // std::cout << resp_space[(deg_ind_x + 1) + (deg_ind_y * 9) - 1] << '\n';
            }
          }
        }
      }
    }

    // for (int i = 0; i < fine; i++)
    // {
    //    for (int j = 0; j < fine; j++)
    //    {
    //        cout << response_spatial_reference[k][i][j] << " ";
    //    }
    //    std::cout << '\n';
    // }
    //
    // std::cout << '\n';

    layer_save = k;
    k_ = k;

  }

  // for (int i = 0; i < fine; i++)
  // {
  //    for (int j = 0; j < fine; j++)
  //    {
  //        cout << response_spatial_reference[k_][i][j] << " ";
  //    }
  //    std::cout << '\n';
  // }

  std::cout << "space 210" << '\n';

  // For more than one iso
  if(iso_count > 1)
  {
    // Storing the index of what the next layer of the response matrix will be
    float over_temp[fine][fine];
    float temp_layer_a[fine][fine];
    float temp_layer_b[fine][fine];
    float temp_layer_c[fine][fine];

    // Iterating through the possibilities for multiple isotopes in space
    for(int iso_ind = 2; iso_ind <= isos; iso_ind++)
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
                // std::cout << "hi" << '\n';
                over_temp[ty][tx] = temp_layer_a[ty][tx] + temp_layer_b[ty][tx];
                // std::cout << "hihi" << '\n';
                // std::cout << k_ << '\n';
                response_spatial_reference[k_][ty][tx] = over_temp[ty][tx];
              }
            }

            // Next layer of 3D response matrix
            k_ = k_ + 1;
          }
          else
          {
            for(int layer_ind_c = layer_ind_b + 1; layer_ind_c < fine * fine; layer_ind_c++)
            {
              // std::cout << "layer a: " << layer_ind_a << '\n';
              // std::cout << "layer b: " << layer_ind_b << '\n';
              // std::cout << "layer c: " << layer_ind_c << '\n';
              // std::cout << '\n' << "k_: " << k_ << '\n';

              for(int ty = 0; ty < fine; ty++)
              {
                for(int tx = 0; tx < fine; tx++)
                {
                  temp_layer_a[ty][tx] = response_spatial_reference[layer_ind_a][ty][tx];
                  temp_layer_b[ty][tx] = response_spatial_reference[layer_ind_b][ty][tx];
                  temp_layer_c[ty][tx] = response_spatial_reference[layer_ind_c][ty][tx];



                  over_temp[ty][tx] = temp_layer_a[ty][tx] + temp_layer_b[ty][tx] + temp_layer_c[ty][tx];

                  response_spatial_reference[k_][ty][tx] = over_temp[ty][tx];

                }
              }

              // Next layer of 3D response matrix
              k_ = k_ + 1;
            }
          }
        }
      }
    }
  }

  std::cout << "space 238" << '\n';

  // New 3D array to copy old "real" contents to (for future_sight)
  float*** response_spatial_reference_expo = new float**[combos];
  for(int i = 0; i < combos; i++)
  {
    response_spatial_reference_expo[i] = new float*[fine - 2];
    for(int j = 0; j < fine - 2; j++)
    {
      response_spatial_reference_expo[i][j] = new float[fine - 2];
    }
  }

  std::cout << "space 251" << '\n';

  // Deleting pseudo-meshes
  if(future_sight == 1)
  {

    for(int lay_i = 0; lay_i < combos; lay_i++)
    {
      for(int i = 0; i < fine - 2; i++)
      {
        for(int j = 0; j < fine - 2; j++)
        {
          response_spatial_reference_expo[lay_i][i][j] = response_spatial_reference[lay_i][i + 1][j + 1];
        }
      }
    }

    std::cout << "space 291" << '\n';

    // Delete "finer" macro-mesh
    for(int i = 0; i < fine; ++i)
    {
      for(int j = 0; j < fine; ++j)
      {
        delete [] response_spatial_reference[i][j];
      }
    }
    delete [] response_spatial_reference;

  }

  /*
      ~~~ END OF RESPONSE ASSEMBLY ~~~
  */

  /*
      ~~~ BEGIN DECONVOLUTION ITERATIONS ~~~
  */

  // Init response matrix for iterations
  float*** response_spatial_reference_expo_iter = new float**[combos];
  for(int i = 0; i < combos; i++)
  {
    response_spatial_reference_expo_iter[i] = new float*[fine - 2];
    for(int j = 0; j < fine - 2; j++)
    {
      response_spatial_reference_expo_iter[i][j] = new float[fine - 2];
    }
  }

  float temp_arr[fine - 2][fine - 2];

  // Overall minimum differences are not reset in iteration
  //  These vectors are to be used for final result determination
  min_min_vec = {10000, 10000, 10000};
  compan_vec = {{-1, -1}, {-1, -1}, {-1, -1}};
  layer_seek = 0;
  iso_found = 0;

  // Each iteration resets the random error
  for(int l = 0; l < iter; l++)
  {

    // "Seed" response matrix for iteration
    for(int lay_i = 0; lay_i < combos; lay_i++)
    {
      for(int i = 0; i < fine - 2; i++)
      {
        for(int j = 0; j < fine - 2; j++)
        {
          response_spatial_reference_expo_iter[lay_i][i][j] = response_spatial_reference_expo[lay_i][i][j];
        }
      }
    }

    // Normalization by max value per layer of macro-mesh
    for(int lay_i = 0; lay_i < combos; lay_i++)
    {
      macro_max = 0;

      for(int i = 0; i < fine - 2; i++)
      {
        for(int j = 0; j < fine - 2; j++)
        {
          // Randomize by error in std dev
          rand_source = (float)rand() / (float)(RAND_MAX / 2) - 1;

          err = sqrt(response_spatial_reference_expo[lay_i][i][j]) * rand_source;


          // Uncomment this block to check that error randomization is working properly

          // int max_err = 0;
          //
          // if(abs(err) > max_err)
          // {
          //   max_err = err;
          //   std::cout << max_err << '\n';
          // }
          // else
          // {
          //   std::cout << "not big: " << err << '\n';
          //   std::cout << "resp: " << response_spatial_reference_expo[lay_i][i][j] << '\n';
          // }

          // if(response_spatial_reference_expo[lay_i][i][j] != response_spatial_reference_expo[lay_i][i][j] + err)
          // {
          //   std::cout << "err:                       " << err << '\n';
          //   std::cout << "iter: " << l << '\n';
          // }

          response_spatial_reference_expo_iter[lay_i][i][j] = response_spatial_reference_expo[lay_i][i][j] + err;

          if(response_spatial_reference_expo_iter[lay_i][i][j] > macro_max)
          {
            macro_max = response_spatial_reference_expo_iter[lay_i][i][j];
          }
        }
      }

      // Normalization by max value loop
      for(int i = 0; i < fine - 2; i++)
      {
        for(int j = 0; j < fine - 2; j++)
        {
          response_spatial_reference_expo_iter[lay_i][i][j] = response_spatial_reference_expo_iter[lay_i][i][j] / macro_max;
        }
      }
    }

    survy_save = 0;
    survx_save = 0;

    // Iterate through survey space to take differences between response and source
    for(int isos_ind = 1; isos_ind <= isos; isos_ind++)
    {
      // Re(set) temp and tracker vars
      min_dif = 1000;

      for(int clr_indy = 0; clr_indy < fine - 2; clr_indy++)
      {
        for(int clr_indx = 0; clr_indx < fine - 2; clr_indx++)
        {
          temp_arr[clr_indy][clr_indx] = 0;
        }
      }

      for(int survy = survy_save; survy < survsizey; survy++)
      {
        for(int survx = survx_save; survx < survsizex; survx++)
        {
          // Grab data from source the same size as macro-mesh
          real_ct = 0;
          for(int grady = 0; grady < fine - 2; grady++)
          {
            for(int gradx = 0; gradx < fine - 2; gradx++)
            {
              if( (survy - ((fine - 2) - 1) / 2 + grady > 0) && (survx - ((fine - 2) - 1) / 2 + gradx > 0) &&
                (survy - ((fine - 2) - 1) / 2 + grady < survsizey) && (survx - ((fine - 2) - 1) / 2 + gradx < survsizex) )
              {
                std::cout << survy - ((fine - 2) - 1) / 2 + grady << '\n';
                temp_arr[grady][gradx] = source_space[survy - ((fine - 2) - 1) / 2 + grady][survx - ((fine - 2) - 1) / 2 + gradx];
                real_ct = real_ct + 1;
              }

            }
          }

          // std::cout << "temp_arr constructed." << '\n';

          // std::cout << "363" << '\n';

          // for (int i = 0; i < fine - 2; i++)
          // {
          //     for (int j = 0; j < fine - 2; j++)
          //     {
          //         cout << temp_arr[i][j] << " ";
          //     }
          //     std::cout << '\n';
          // }

          // std::cout << "373" << '\n';

          // Ensure that survey location has enough data to comapre
          if(real_ct > round( ((fine - 2) * (fine - 2)) / 2 ))
          {
            std::cout << "live one" << '\n';

            // Normalize survey macro-mesh to highest value
            max_ele = 0;
            for(int grady = 0; grady < fine - 2; grady++)
            {
              for(int gradx = 0; gradx < fine - 2; gradx++)
              {
                // Find max data value
                if(temp_arr[grady][gradx] > max_ele)
                {
                  max_ele = temp_arr[grady][gradx];
                }
              }
            }

            // Divide survey macro-mesh by its max data value
            for(int grady = 0; grady < fine - 2; grady++)
            {
              for(int gradx = 0; gradx < fine - 2; gradx++)
              {
                temp_arr[grady][gradx] = temp_arr[grady][gradx] / max_ele;
              }
            }

            // Take difference of data
            for(int layer_thru = 0; layer_thru < combos; layer_thru++)
            {
              dif_sum = 0;
              resp_sum = 0;

              for(int grady = 0; grady < fine - 2; grady++)
              {
                for(int gradx = 0; gradx < fine - 2; gradx++)
                {
                  if(future_sight == 1)
                  {
                    dif_sum = dif_sum + abs(response_spatial_reference_expo_iter[layer_thru][grady][gradx] - temp_arr[grady][gradx]);
                    resp_sum = resp_sum + response_spatial_reference_expo_iter[layer_thru][grady][gradx];
                  }
                  else
                  {
                    dif_sum = dif_sum + abs(response_spatial_reference[layer_thru][grady][gradx] - temp_arr[grady][gradx]);
                    resp_sum = resp_sum + response_spatial_reference[layer_thru][grady][gradx];
                  }
                }
              }
              if(dif_sum < min_dif)
              {
                min_dif = dif_sum;
                // std::cout << min_dif << '\n';
                if(dif_sum < min_min_vec[isos_ind - 1])
                {
                  min_min_vec[isos_ind - 1] = dif_sum;
                  // std::cout << min_min_vec[isos_ind - 1] << '\n';
                  conf = 1 - dif_sum / resp_sum;
                  // std::cout << "conf: " << conf << '\n';
                  layer_seek = layer_thru;
                  if(conf > conf_levels[6 - (fine - 2)])
                  {
                    // std::cout << "iso #" << isos_ind << ", iter: " << l << '\n';
                    survx_save = survx + 1;
                    // std::cout << "survx_save: " << survx_save << '\n';
                    survy_save = survy + 1;
                    // std::cout << "survy_save: " << survy_save << '\n';
                    layer_thru = combos;
                    survx = survsizex;
                    survy = survsizey;
                    iso_found = 1;
                  }
                }
                // std::cout << "hit: " << dif_sum << '\n';
              }

            }

            max_ele_repo = 0;
            savey = 0;
            savex = 0;
            for(int grady = 0; grady < fine - 2; grady++)
            {
              for(int gradx = 0; gradx < fine - 2; gradx++)
              {
                if(future_sight == 1)
                {
                  // Find max data value
                  if(response_spatial_reference_expo_iter[layer_seek][grady][gradx] > max_ele_repo)
                  {
                    max_ele_repo = response_spatial_reference_expo_iter[layer_seek][grady][gradx];
                    savey = grady;
                    savex = gradx;
                  }
                }
                else
                {
                  if(response_spatial_reference[layer_seek][grady][gradx] > max_ele_repo)
                  {
                    max_ele_repo = response_spatial_reference[layer_seek][grady][gradx];
                    savey = grady;
                    savex = gradx;
                  }
                }
              }
            }

            // Write predicted locations to temp vectors
            y_end.push_back(survy - savey);
            x_end.push_back(survx - savex);

          }
          // std::cout << "x space: " << survx << '\n';
        }
        // std::cout << "y space: " << survy << '\n';
      }
      // std::cout << "iso ind: " << isos_ind << '\n';
      if(iso_found == 0)
      {
        std::cout << "No source located! Try lowering confidence required." << '\n';
        exit(1);
      }
    }
    std::cout << "iter: " << l << '\n';
  }

  end_ind.push_back(y_end);
  end_ind.push_back(x_end);

  // for (int i = 0; i < end_ind.size(); i++)
  // {
  //     for (int j = 0; j < end_ind[i].size(); j++)
  //     {
  //         cout << end_ind[i][j] << endl;
  //     }
  // }

  std::cout << "484" << '\n';

  if(future_sight == 0)
  {
    for(int i = 0; i < fine; ++i)
    {
      for(int j = 0; j < fine; ++j)
      {
        delete [] response_spatial_reference[i][j];
      }
    }
    delete [] response_spatial_reference;
  }

  std::cout << "End spatial decon." << '\n';

}

spatial_decon::~spatial_decon() {}
