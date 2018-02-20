// Reading in known isotope spectra
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>

#include "resp_read.hh"

using namespace std;

float resp_read::resp_read()
{

  spectra_names.push_back("5m_BG_2-15-18.txt");
  spectra_names.push_back("5m_Am241_10micCi.txt");
  spectra_names.push_back("5m_Co60.txt");
  spectra_names.push_back("5m_Cs137_2-15-18.txt");
  spectra_names.push_back("5m_depU_2-15-18.txt");

  ifstream inFile;

  for(k=0; k < 5; k++)
  {

    inFile.open(spectra_names[k]);
    if (!inFile)
    {
      cout << "Unable to open file";
      exit(1); // terminate with error
    }

    while (inFile >> x)
    {
      temp.push_back(x);
    }

    resp_mat_read.push_back(temp);

    inFile.close();

  }

  return resp_mat_read;
}

resp_read::~resp_read() {}
