// Reading in known isotope spectra
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>

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

  // Space response
  space_names.push_back("1-1_space_3-09.txt");
  space_names.push_back("1-4_space_3-09.txt");
  space_names.push_back("1-7_space_3-09.txt");
  space_names.push_back("1-10_space_3-09.txt");
  space_names.push_back("1-13_space_3-09.txt");
  space_names.push_back("1-16_space_3-09.txt");
  space_names.push_back("4-1_space_3-09.txt");
  space_names.push_back("4-4_space_3-09.txt");
  space_names.push_back("4-7_space_3-09.txt");
  space_names.push_back("4-10_space_3-09.txt");
  space_names.push_back("4-13_space_3-09.txt");
  space_names.push_back("4-16_space_3-09.txt");
  space_names.push_back("7-1_space_3-09.txt");
  space_names.push_back("7-4_space_3-09.txt");
  space_names.push_back("7-7_space_3-09.txt");
  space_names.push_back("7-10_space_3-09.txt");
  space_names.push_back("7-13_space_3-09.txt");
  space_names.push_back("7-16_space_3-09.txt");
  space_names.push_back("10-1_space_3-09.txt");
  space_names.push_back("10-4_space_3-09.txt");
  space_names.push_back("10-7_space_3-09.txt");
  space_names.push_back("10-10_space_3-09.txt");
  space_names.push_back("10-13_space_3-09.txt");
  space_names.push_back("10-16_space_3-09.txt");

  // space_names.push_back("background_space_3-09.txt");
  // space_names.push_back("backgroundPost_space_3-09.txt");

  ifstream inFile;

  // Spectral response assemby
  for(int k=0; k < 5; k++)
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

  // Spatial response assembly

  count_time = 300;

  for(k=0; k < 24; k++)
  {

    inFile.open(space_names[k]);
    if (!inFile)
    {
      cout << "Unable to open file";
      exit(1); // terminate with error
    }

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
