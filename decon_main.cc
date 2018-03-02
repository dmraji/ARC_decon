#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

#include "data_read.hh"
#include "gold_decon.hh"
#include "resp_read.hh"
#include "histo.hh"
#include "sim_source.hh"

using namespace std;

vector< vector<float> > data_mat;
int data_mat_len;

vector<float> data_sim;
vector< vector<float> > resp_mat_read;
float source_decon[4096];

int main(int argc, char** argv) {

  // cout << "hi.\n";

  // Build response matrix from experimental spectra
  int chs = 4096;
  int num_spectra = 6;

  // Call to function to read response files
  resp_read resp;

  float** response_matrix = new float*[chs];
  for(int i = 0; i < chs; i++)
  {
    response_matrix[i] = new float[num_spectra];
  }

  cout << "ln 39 main." << endl;

  for(int j = 0; j < num_spectra; j++)
  {
    for(int i = 0; i < chs; i++)
    {
      response_matrix[i][j] = resp_mat_read[i][j];
    }
  }

  cout << "ln 49 main." << endl;

  // Obtain vector of vectors via data read-in
  data_read data;

  // change this to stop using sim source
  int sim = 1;
  if(sim == 1)
  {
    sim_source sim;
  }

  cout << "ln 61 main." << endl;

  float* source_mat = new float[data_mat_len];

  if(sim == 0)
  {
    for (int dat_ind = 0; dat_ind < data_mat_len; dat_ind++)
    {
      source_mat[dat_ind] = data_mat[dat_ind][8];
    }
  }
  else
  {
    for (int dat_ind = 0; dat_ind < data_mat_len; dat_ind++)
    {
      source_mat[dat_ind] = data_sim[dat_ind];
    }
  }

  histo makeHist(source_mat,
                 data_mat_len);

  int num_iter = 10000;
  int num_rep = 10;
  double boost = 10;

  gold_decon callIt(response_matrix,
                    source_decon,
                    chs,
                    num_spectra,
                    num_iter,
                    num_rep,
                    boost
                    );

  for(int printer = 0; printer < 4096; printer++)
  {
    cout << source_decon[printer] << endl;
  }

  /*
  spatial_decon spaceOut(
    // unfinished
  );
  */

  for(int i = 0; i < chs; ++i) {
    delete [] response_matrix[i];
  }

  delete [] response_matrix;

  return 0;
}
