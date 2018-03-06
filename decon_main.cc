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
  int num_spectra = 5;

  // Call to function to read response files
  resp_read resp;

  // initialize 2D dynamic array for response
  float** response_matrix = new float*[num_spectra];
  for(int i = 0; i < num_spectra; i++)
  {
    response_matrix[i] = new float[chs];
  }

  cout << "ln 39 main." << endl;

  for(int j = 0; j < chs; j++)
  {
    for(int i = 0; i < num_spectra; i++)
    {
      response_matrix[i][j] = resp_mat_read[i][j];


    }
  }

  cout << "ln 49 main." << endl;

  // change this to stop using sim source
  int sim = 1;
  if(sim == 1)
  {
    sim_source sim;
  }
  else
  {
    // Obtain vector of vectors via data read-in
    data_read data;
  }

  cout << "ln 61 main." << endl;

  float* source_mat = new float[data_mat_len];

  if(sim == 0)
  {
    for (int dat_ind = 0; dat_ind < data_mat_len; dat_ind++)
    {
      source_mat[dat_ind] = data_mat[dat_ind][8];
    }
    histo makeHist(source_mat,
                   data_mat_len);
  }
  else
  {
    data_mat_len = 4096;
    for (int dat_ind = 0; dat_ind < data_mat_len; dat_ind++)
    {
      source_decon[dat_ind] = data_sim[dat_ind];
    }
  }

  int num_iter = 10000;
  int num_rep = 10;
  double boost = 2;

  // !!
  //
  // NORMALIZE RESP BY ITENTSITY (TOTAL COUNTS) AND SOURCE BY TIME (2*3600)
  //
  // !!
  gold_decon callIt(response_matrix,
                    source_decon,
                    chs,
                    num_spectra,
                    num_iter,
                    num_rep,
                    boost
                    );

  // for(int printer = 0; printer < 4096; printer++)
  // {
  //   cout << source_decon[printer] << endl;
  // }

  // spatial_decon spaceOut(// array of number of mesh points with counts as values
  //                        );


  std::cout << "main ln 110" << '\n';

  for(int i = 0; i < num_spectra; ++i) {
    delete [] response_matrix[i];
  }

  delete [] response_matrix;

  return 0;
}
