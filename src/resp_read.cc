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

  // spectra_names.push_back("2_21_hour_BG.txt");
  spectra_names.push_back("2_21_am241.txt");
  spectra_names.push_back("2_21_co60.txt");
  spectra_names.push_back("2_21_cs137.txt");
  spectra_names.push_back("2_21_natU.txt");
  spectra_names.push_back("2_21_ra226.txt");

  ifstream inFile;

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

    for(int cs = 0; cs < 4096; cs++)
    {
      temp[cs] = temp[cs] / spect_sum;
    }

    resp_mat_read.push_back(temp);

    inFile.close();

  }

  cout << "end resp read." << endl;

  // return resp_mat_read;
}

resp_read::~resp_read() {}
