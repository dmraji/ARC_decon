#include <string>
#include <iostream>
#include <vector>

#include "histo.hh"

using namespace std;

histo::histo(float *source_matrix,
             int data_len)
{

  cout << "begin histo." << endl;

  source_decon[4096] = {0};

  // cout << data_len << endl;

  // for(int wksp_print=0; wksp_print < data_len; wksp_print++)
  // {
  //   std::cout << source_matrix[wksp_print] << '\n';
  // }

  for(int p=0; p<data_len-1; p++)
  {
    temp_var = int(source_matrix[p] + 0.5);

    source_decon[temp_var] = source_decon[temp_var] + 1;

    std::cout << source_decon[temp_var] << '\n';
  }

  cout << "end histo." << endl;

  // return source_decon;
}

histo::~histo() {}
