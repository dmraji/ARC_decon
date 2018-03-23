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

  // Spectral response
  // spectra_names.push_back("2_21_hour_BG.txt");
  spectra_names.push_back("2_21_am241.txt");
  spectra_names.push_back("2_21_co60.txt");
  spectra_names.push_back("2_21_cs137.txt");
  spectra_names.push_back("2_21_natU.txt");
  spectra_names.push_back("2_21_ra226.txt");

  ifstream inFile;

  // Spectral response assemby
  for(int k = 0; k < 5; k++)
  {

    inFile.open(spectra_names[k]);
    if (!inFile)
    {
      cout << "Unable to open file";
      exit(1); // terminate with error
    }

    temp.clear();
    spect_sum = 0;

    while (inFile >> x)
    {
      temp.push_back(x);
      spect_sum = spect_sum + x;
    }

    std::cout << spect_sum << '\n';

    // Normalization of spectra
    for(int cs = 0; cs < 4096; cs++)
    {
      temp[cs] = temp[cs] / spect_sum;
    }

    resp_mat_read.push_back(temp);

    inFile.close();

  }

  // Space response

  int xi = {1, 4, 7, 10, 13, 16,  19,  22,  25};
  int yi = {7, 4, 1, -2, -5, -8, -11, -14, -17};

  for(int yn = 0; yn < 9; yn++)
  {
    for(int xn = 0; xn < 9; xn++)
    {
      space_n = xi[xn] << "-" << yi[yn] << "_space_3-09.txt";
      inFile.open(space_n);
      if(!inFile)
      {
        space_n = xi[xn] << "-" << yi[yn] << "_space_3-14.txt";
        inFile.open(space_n);
        if(!inFile)
        {
          space_n = xi[xn] << "-" << yi[yn] << "_space_3-22.txt";
          inFile.open(space_n);
          if(!inFile)
          {

            // Provision to fill in missing data by "flipping" equidistant data
            if(yn < 1 && xn < 19)
            {
              space_n = xi[yn] << "-" << yi[xn] << "_space_3-09.txt";
              inFile.open(space_n);
              if(!inFile)
              {
                space_n = xi[yn] << "-" << yi[xn] << "_space_3-14.txt";
                inFile.open(space_n);
                if(!inFile)
                {
                  space_n = xi[yn] << "-" << yi[xn] << "_space_3-22.txt";
                  inFile.open(space_n);
                  if(!inFile)
                  {
                    cout << "Unable to open file";
                    exit(1); // terminate with error
                  }
                  else{
                    space_names.push_back(space_n);
                    inFile.close();
                  }
                }
                else{
                  space_names.push_back(space_n);
                  inFile.close();
                }
              }
              else{
                space_names.push_back(space_n);
                inFile.close();
              }
            }
            else
            {
              cout << "Unable to open file";
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

  for(k = start_k; k < 81; k++)
  {

    inFile.open(space_names[k]);

    spect_sum = 0;

    while (inFile >> x)
    {
      spect_sum = spect_sum + x;
    }

    std::cout << spect_sum << '\n';

    // Vector of count rates at mesh points from detector
    resp_space_mat_read.push_back(spect_sum / count_time);

    inFile.close();

  }

  cout << "end resp read." << endl;

  // return resp_mat_read;
}

resp_read::~resp_read() {}
