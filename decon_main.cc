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
float source_decon[4096];

int main(int argc, char** argv) {

  cout << "hi.\n";

  // Build response matrix from experimental spectra
  int chs = 4096;
  int num_spectra = 5;

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

  float* source_mat = new float[];

  for (int dat_ind = 0; dat_ind < data_mat.size(); dat_ind++)
  {
    source_mat[dat_ind] = data_mat[dat_ind][8];
  }

  histo makeHist(source_mat);

  int num_iter = 10000;
  int num_rep = 10;
  double boost = 10;

  gold_decon callIt(response_matrix,
                    spectrum_vec,
                    chs,
                    num_spectra,
                    num_iter,
                    num_rep,
                    boost
                    );

  spatial_decon spaceOut(
    // unfinished
  );

  for(int i = 0; i < chs; ++i) {
    delete [] response_matrix[i];
  }

  delete [] response_matrix;

  return 0;
}
