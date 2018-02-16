#include <string>
#include <iostream>
#include <vector>

#include "data_read.hh"
#include "gold_decon.hh"

using namespace std;

vector< vector<float> > data_mat;

int main(int argc, char** argv) {

  cout << "hi.\n";

  // Obtain vector of vectors via data read-in
  data_read data;

  // Take the energy values from the data and build a response matrix to put
  //   through Gold's energy deconvolution algorithm
  vector<float> resp_mat;
  for (int dat_ind = 0; dat_ind < data_mat.size(); dat_ind++)
  {
    resp_mat.push_back(data_mat[dat_ind][8]);
  }

  for (int i = 0; i < resp_mat.size(); i++)
  {
      cout << resp_mat[i] << endl;
  }

    gold_decon(*source,
               resp_mat,
               sizex,
               sizey,
               numberIterations,
               numberRepetitions,
               boost
             );

  return 0;
}
