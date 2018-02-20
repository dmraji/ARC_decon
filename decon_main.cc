#include <string>
#include <iostream>
#include <iomanip>
#include <ifstream>
#include <vector>

#include "data_read.hh"
#include "gold_decon.hh"

using namespace std;

vector< vector<float> > data_mat;
vector< vector<float> > resp_mat_read;

int main(int argc, char** argv) {

  cout << "hi.\n";

  // Build response matrix from experimental spectra
  int chs = 4096;
  int num_spectra = 2;

  // Call to function to read response files
  resp_read resp;

  float** response_matrix = new float*[chs];
  for(int i = 0; i < chs; i++)
  {
    response_matrix[i] = new int[num_spectra];
  }

  for(int j = 0; j < num_spectra; j++)
  {
    for(int i = 0; i < chs; i++)
    {
      response_matrix[i][j] = resp_mat_read[i][j];
    }
  }

  // Obtain vector of vectors via data read-in
  data_read data;

  // Take the energy values from the data and build a response matrix to put
  //   through Gold's energy deconvolution algorithm
  int x;
  ifstream inFile;

  inFile.open("source_sim.txt");
  if (!inFile) {
      cout << "Unable to open file";
      exit(1); // terminate with error
  }



  inFile.close();

  for (int dat_ind = 0; dat_ind < data_mat.size(); dat_ind++)
  {
    source_mat.push_back(data_mat[dat_ind][8]);
  }

  for (int i = 0; i < source_mat.size(); i++)
  {
      cout << source_mat[i] << endl;
  }

    gold_decon callit("VECOFVEC_RESP_VAR",
                      source_mat,
                      "INT_NUM_SPECTRA_IN_REC_MAT",
                      "INT_NUM_CHANNELS_IN_REC_MAT",
                      "INT_NUM_ITER",
                      "INT_NUM_REPET",
                      "DOUBLE_BOOST"
                      );

  for(int i = 0; i < chs; ++i) {
    delete [] response_matrix[i];
  }
  delete [] response_matrix;

  return 0;
}
