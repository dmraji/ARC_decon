#include <string>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <algorithm>

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
                             int future_sight,
                             int sim,
                             int source_num_sim
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

    // Validate fineness selection considering source size
    while(survsizey * survsizex < (fine - 1) * (fine - 1))
    {
      if(fine >= 7)
      {
        std::cout << "Not enough source data to deconvolve. Lowering fineness..." << '\n';
        fine = fine - 2;
      }
      else
      {
        std::cout << "Not enough source data. Terminating." << '\n';
        exit(1);
      }
    }
  }

  conf_levels = {0.95, 0.97, 0.98};

  /* Notes on response space:

    Is a 1D array, laid out in progressing depth from detector;

    Is laid out in terms of distance from det as follows (in cm):


  */

  /* Response cm distance vector, ordered to index in conj. with counts
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

  */

  /*  Number of layers in the response matrix

      We assume that, at max, each macro-mesh will see up to 3 isotopes

  */

  // For sim stuff
  if(sim == 1)
  {
    iso_count = source_num_sim;
  }

  if(iso_count > 2)
  {
    isos = 2;
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

  // Set up gradient deducer
  for(int d = 0; d < 9; d++)
  {
    degrees.push_back(d);
  }

  /*
      ~~~ BEGIN RESPONSE ASSEMBLY ~~~
  */

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

  // For more than one iso
  if(isos > 1)
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
  for(int i = 0; i < iso_count + 10; i++)
  {
    compan_vec.push_back({-1, -1});
    min_min_vec.push_back(10000);
  }
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

    iso_pass = 0;
    pre_pass = 0;

    survy_save = 0;
    survx_save = 0;
    survskipx.clear();
    survskipy.clear();
    for(int i = 0; i < 9; i++)
    {
      survskipx.push_back(0);
      survskipy.push_back(0);
    }
    no_rep = 0;

    save_survy_vec.clear();
    save_survx_vec.clear();
    save_survy2_vec.clear();
    save_survx2_vec.clear();

    // Iterate through survey space to take differences between response and source
    for(int isos_ind = 1; isos_ind <= iso_count; isos_ind++)
    {

      if(iso_pass == 1)
      {
        pre_pass = 0;
        iso_pass = 0;
        goto label_i;
      }

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
                // std::cout << survy - ((fine - 2) - 1) / 2 + grady << '\n';
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

          // Ensure that survey location has enough data to comapre
          if(real_ct > round( ((fine - 2) * (fine - 2)) / 2 ))
          {
            // std::cout << "live one" << '\n';

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
            for(int layer_thru = 0; layer_thru < fine * fine; layer_thru++)
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
                  if(conf > conf_levels[2 - ((fine - 2) - 1) / 2 - 1])
                  {
                    // std::cout << "iso #" << isos_ind << ", iter: " << l << '\n';
                    survx_save = survx + 1;
                    // std::cout << "survx_save: " << survx_save << '\n';
                    survy_save = survy + 1;
                    // std::cout << "survy_save: " << survy_save << '\n';

                    // Record keeping block
                    max_ele_repo = 0;
                    max2_ele_repo = 0;
                    savey = -1;
                    savex = -1;
                    savey2 = -1;
                    savex2 = -1;

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

                            for(int cc = 0; cc < save_survy_vec.size(); cc++)
                            {
                              if(save_survy_vec[cc] == survy + savey - floor((fine - 2) / 2) &&
                                 save_survx_vec[cc] == survx + savex - floor((fine - 2) / 2))
                              {
                                pre_pass = 0;
                                savey = -1;
                                savex = -1;
                              }
                            }
                            for(int cc = 0; cc < save_survy2_vec.size(); cc++)
                            {
                              if(save_survy2_vec[cc] == survy + savey - floor((fine - 2) / 2) &&
                                 save_survx2_vec[cc] == survx + savex - floor((fine - 2) / 2))
                              {
                                pre_pass = 0;
                                savey = -1;
                                savex = -1;
                              }
                              else
                              {
                                pre_pass = 1;
                                save_survy_vec.push_back(survy + savey - floor((fine - 2) / 2));
                                save_survx_vec.push_back(survx + savex - floor((fine - 2) / 2));
                              }
                            }

                          }
                        }
                        // else
                        // {
                        //   if(response_spatial_reference[layer_seek][grady][gradx] > max_ele_repo)
                        //   {
                        //     max_ele_repo = response_spatial_reference[layer_seek][grady][gradx];
                        //     savey = grady;
                        //     savex = gradx;
                        //   }
                        // }
                      }
                    }

                    // If layer was built from multiple seeds
                    // if(layer_seek > fine * fine)
                    // {
                    //
                    //   for(int grady = 0; grady < fine - 2; grady++)
                    //   {
                    //     for(int gradx = 0; gradx < fine - 2; gradx++)
                    //     {
                    //       if(future_sight == 1)
                    //       {
                    //         // Find max data value
                    //         if(response_spatial_reference_expo_iter[layer_seek][grady][gradx] > max2_ele_repo &&
                    //            response_spatial_reference_expo_iter[layer_seek][grady][gradx] < max_ele_repo)
                    //         {
                    //           max2_ele_repo = response_spatial_reference_expo_iter[layer_seek][grady][gradx];
                    //           savey2 = grady;
                    //           savex2 = gradx;
                    //
                    //           // Verfiy that secondaries aren't repeated
                    //           for(int cc = 0; cc < save_survy_vec.size(); cc++)
                    //           {
                    //             if(save_survy_vec[cc] == survy + savey2 - floor((fine - 2) / 2) &&
                    //                save_survx_vec[cc] == survx + savex2 - floor((fine - 2) / 2))
                    //             {
                    //               savey2 = -1;
                    //               savex2 = -1;
                    //             }
                    //           }
                    //           for(int cc = 0; cc < save_survy2_vec.size(); cc++)
                    //           {
                    //             if(save_survy2_vec[cc] == survy + savey2 - floor((fine - 2) / 2) &&
                    //                save_survx2_vec[cc] == survx + savex2 - floor((fine - 2) / 2))
                    //             {
                    //               savey2 = -1;
                    //               savex2 = -1;
                    //             }
                    //             else
                    //             {
                    //               if(pre_pass == 1)
                    //               {
                    //                 iso_pass = 1;
                    //               }
                    //               save_survy2_vec.push_back(survy + savey2 - floor((fine - 2) / 2));
                    //               save_survx2_vec.push_back(survx + savex2 - floor((fine - 2) / 2));
                    //             }
                    //           }
                    //
                    //         }
                    //       }
                    //       else
                    //       {
                    //         if(response_spatial_reference[layer_seek][grady][gradx] > max2_ele_repo &&
                    //            response_spatial_reference[layer_seek][grady][gradx] < max_ele_repo)
                    //         {
                    //           max2_ele_repo = response_spatial_reference[layer_seek][grady][gradx];
                    //           savey2 = grady;
                    //           savex2 = gradx;
                    //
                    //         }
                    //       }
                    //     }
                    //   }
                    // }

                    // Write predicted locations to temp vectors
                    // for(int c = 0; c < compan_vec.size(); c++)
                    // {
                    //   if((compan_vec[c][0] == survy + savey - floor((fine - 2) / 2) &&
                    //       compan_vec[c][1] == survx + savex - floor((fine - 2) / 2)) ||
                    //      (abs(compan_vec[c][0] - survy + savey2 - floor((fine - 2) / 2)) +
                    //       abs(compan_vec[c][1] - survx + savex2 - floor((fine - 2) / 2)) < 3))
                    //   {
                    //     no_rep = 1;
                    //     if(savey2 > -1)
                    //     {
                    //       if((compan_vec[c][0] == survy + savey2 - floor((fine - 2) / 2) &&
                    //           compan_vec[c][1] == survx + savex2 - floor((fine - 2) / 2)) ||
                    //          (abs(compan_vec[c][0] - survy + savey2 - floor((fine - 2) / 2)) +
                    //           abs(compan_vec[c][1] - survx + savex2 - floor((fine - 2) / 2)) < 3))
                    //       {
                    //         no_rep = 1;
                    //       }
                    //     }
                    //   }
                    // }
                    if(no_rep == 0 && savey > -1)
                    {

                      // Skip already found source meshes
                      for(int com = 0; com < compan_vec.size(); com++)
                      {
                        passer = 0;
                        if((compan_vec[com][0] == survy + savey - floor((fine - 2) / 2)) && (compan_vec[com][1] == survx + savex - floor((fine - 2) / 2)))
                        {
                          goto label;
                        }
                      }

                      compan_vec[isos_ind - 1][0] = survy + savey - floor((fine - 2) / 2);
                      std::cout << "y res: " << survy + savey - floor((fine - 2) / 2) << '\n';

                      compan_vec[isos_ind - 1][1] = survx + savex - floor((fine - 2) / 2);
                      std::cout << "x res: " << survx + savex - floor((fine - 2) / 2) << '\n';

                      if(savey2 > -1)
                      {
                        // compan_vec[isos_ind][0] = survy + savey2 - floor((fine - 2) / 2);
                        // std::cout << "y res2: " << survy + savey2 - floor((fine - 2) / 2) << '\n';
                        //
                        // compan_vec[isos_ind][1] = survx + savex2 - floor((fine - 2) / 2);
                        // std::cout << "x res2: " << survx + savex2 - floor((fine - 2) / 2) << '\n';
                      }

                      for(int iii = 0; iii < 3; iii++)
                      {
                        for(int iiii = 0; iiii < 3; iiii++)
                        {
                          survskipy[iiii + (iii * 3)] = survy + savey - floor((fine - 2) / 2) - 1 + iii;
                          survskipx[iiii + (iii * 3)] = survx + savex - floor((fine - 2) / 2) - 1 + iii;
                        }
                      }

                      layer_thru = combos;
                      // survx = survsizex;
                      // survy = survsizey;
                      iso_found = 1;

                    }

                  }

                }
                // std::cout << "hit: " << dif_sum << '\n';
              }

            }

          }
          // std::cout << "x space: " << survx << '\n';
        }

        // Recieves skip from crater
        label:
          ;
        // std::cout << "y space: " << survy << '\n';
      }
      std::cout << "iso ind: " << isos_ind << '\n';
      if(iso_found == 0)
      {
        std::cout << "No source located! Try lowering confidence required." << '\n';
        exit(1);
      }

      // Recieves skip from double iso find
      label_i:
        ;

    }
    std::cout << "iter: " << l << '\n';
  }

  for(int e = 0; e < iso_count; e++)
  {
    end_ind.push_back(compan_vec[e]);
  }

  // for (int i = 0; i < end_ind.size(); i++)
  // {
  //     for (int j = 0; j < end_ind[i].size(); j++)
  //     {
  //         cout << end_ind[i][j] << endl;
  //     }
  // }


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
  else
  {
    for(int i = 0; i < fine - 2; ++i)
    {
      for(int j = 0; j < fine - 2; ++j)
      {
        delete [] response_spatial_reference_expo_iter[i][j];
      }
    }
    delete [] response_spatial_reference_expo_iter;
  }

  std::cout << "End spatial decon." << '\n';

}

spatial_decon::~spatial_decon() {}
