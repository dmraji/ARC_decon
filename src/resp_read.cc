// Reading in known isotope spectra
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "resp_read.hh"

using namespace std;

resp_read::resp_read()
{

  std::cout << "Begin response read." << '\n';

  // Spectral response
  spectra_names.push_back("3_1_bkg.txt");
  spectra_names.push_back("3_1_co57.txt");
  spectra_names.push_back("3_1_pu239.txt");
  spectra_names.push_back("3_2_ba133.txt");
  spectra_names.push_back("3_2_cf252.txt");
  spectra_names.push_back("3_2_cs137.txt");
  spectra_names.push_back("3_2_th232.txt");
  spectra_names.push_back("3_28_bkg_12hr.txt");
  spectra_names.push_back("3_28_am241_2hr.txt");
  spectra_names.push_back("3_28_ntrlU_2hr.txt");
  spectra_names.push_back("3_28_u235_2hr.txt");
  spectra_names.push_back("3_29_co6_2hr.txt");
  spectra_names.push_back("3_29_ra226_2hr.txt");

  // spectra_names.push_back("2_21_hour_BG.txt");
  // spectra_names.push_back("2_21_am241.txt"); // 100% intrin
  // spectra_names.push_back("2_21_co60.txt"); // 45% intrin
  // spectra_names.push_back("2_21_natU.txt"); // 40% intrin
  // spectra_names.push_back("2_21_ra226.txt"); // 60% intrin
  // spectra_names.push_back("2_21_cs137.txt"); // 58% intrin

  // eff.push_back(0);
  // eff.push_back(1);
  // eff.push_back(0.45);
  // eff.push_back(0.58);
  // eff.push_back(0.40);
  // eff.push_back(0.60);

  ifstream inFile;

  // Spectral response assemby
  for(int k = 0; k < 13; k++)
  {

    inFile.open(spectra_names[k]);
    if (!inFile)
    {
      cout << "Unable to open response file. Terminating." << endl;
      std::cout << k << '\n';
      exit(1); // terminate with error
    }

    temp.clear();
    spect_sum = 0;

    while (inFile >> x)
    {
      if(k == 0)
      {
        bg.push_back(floor(x / 6));
        // std::cout << x << '\n';
      }
      else if(k == 7)
      {
        bg2.push_back(floor(x / 6));
      }
      else
      {
        temp.push_back(x);
      }
    }

    if(k != 0 && k != 7)
    {
      for(int l = 0; l < 4096; l++)
      {
        // temp[l] = temp[l] / eff[k];
        // std::cout << temp[l] << '\n';
        if(k < 7)
        {
          temp[l] = temp[l] - bg[l];
        }
        else
        {
          temp[l] = temp[l] - bg2[l];
        }

        if(temp[l] < 0)
        {
          temp[l] = 0;
        }
        // std::cout << temp[l] << '\n';
        spect_sum = spect_sum + temp[l];
      }
    }

    // std::cout << spect_sum << '\n';

    // Normalization of spectra
    if(k != 0 && k != 7)
    {
      for(int cs = 0; cs < 4096; cs++)
      {
        temp[cs] = temp[cs] / spect_sum;
      }

      resp_mat_read.push_back(temp);
    }

    inFile.close();

  }

  std::cout << "Pass spectra read." << '\n';

  // Space response

  int xi[9] = {1, 4, 7, 10, 13, 16,  19,  22,  25};
  int yi[9] = {7, 4, 1, -2, -5, -8, -11, -14, -17};

  for(int yn = 0; yn < 9; yn++)
  {
    for(int xn = 0; xn < 9; xn++)
    {
      space_n = to_string(xi[xn]) + "-" + to_string(yi[yn]) + "_space_3-09.txt";
      inFile.open(space_n);
      if(!inFile)
      {
        space_n = to_string(xi[xn]) + "-" + to_string(yi[yn]) + "_space_3-14.txt";
        inFile.open(space_n);
        if(!inFile)
        {
          space_n = to_string(xi[xn]) + "-" + to_string(yi[yn]) + "_space_3-22.txt";
          inFile.open(space_n);
          if(!inFile)
          {

            // Provision to fill in missing data by "flipping" equidistant data
            if(yi[yn] < 1)
            {
              space_n = to_string(xi[yn]) + "-" + to_string(yi[xn]) + "_space_3-09.txt";
              inFile.open(space_n);
              if(!inFile)
              {
                space_n = to_string(xi[yn]) + "-" + to_string(yi[xn]) + "_space_3-14.txt";
                inFile.open(space_n);
                if(!inFile)
                {
                  space_n = to_string(xi[yn]) + "-" + to_string(yi[xn]) + "_space_3-22.txt";

                  inFile.open(space_n);
                  if(!inFile)
                  {
                    cout << "Unable to open file" << endl;
                    exit(1); // terminate with error
                  }
                  else
                  {
                    space_names.push_back(space_n);
                    inFile.close();
                  }
                }
                else
                {
                  space_names.push_back(space_n);
                  inFile.close();
                }
              }
              else
              {
                space_names.push_back(space_n);
                inFile.close();
              }
            }
            else
            {
              cout << "Unable to open file" << endl;
              exit(1); // terminate with error
            }
          }
          else
          {
            space_names.push_back(space_n);
            inFile.close();
          }
        }
        else
        {
          space_names.push_back(space_n);
          inFile.close();
        }
      }
      else
      {
        space_names.push_back(space_n);
        inFile.close();
      }
    }
  }

  // space_names.push_back("background_space_3-09.txt");
  // space_names.push_back("backgroundPost_space_3-09.txt");

  // Spatial response assembly

  count_time = 300;
  start_k = 0;

  for(int k = start_k; k < 81; k++)
  {

    inFile.open(space_names[k]);

    spect_sum = 0;

    while (inFile >> x)
    {
      spect_sum = spect_sum + x;
    }

    // std::cout << spect_sum << '\n';

    // Vector of count rates at mesh points from detector
    resp_space_mat_read.push_back(spect_sum / count_time);

    inFile.close();

  }

  std::cout << resp_space_mat_read.size() << '\n';

  cout << "end resp read." << endl;

  // return resp_mat_read;
}

resp_read::~resp_read() {}
