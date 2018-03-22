#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <ctime>
#include <cstdio>

#include "data_read.hh"
#include "gold_decon.hh"
#include "resp_read.hh"
#include "histo.hh"
#include "sim_source.hh"
#include "spatial_decon.hh"

using namespace std;

vector< vector<float> > data_mat;
int data_mat_len;

vector<float> data_sim;
float space_sim[24];

vector< vector<float> > resp_mat_read;
vector<float> resp_space_mat_read;
float source_decon[4096];

int iso_count;

int main(int argc, char** argv) {

  std::clock_t start;
  double duration;

  start = std::clock();

  // cout << "hi.\n";

  // Build response matrix from experimental spectra
  int chs = 4096;
  int num_spectra = 5;

  // Call to function to read response files
  resp_read resp;

  // initialize 2D dynamic array for responses
  float** response_matrix = new float*[num_spectra];
  for(int i = 0; i < num_spectra; i++)
  {
    response_matrix[i] = new float[chs];
  }

  for(int j = 0; j < chs; j++)
  {
    for(i = 0; i < num_spectra; i++)
    {
      response_matrix[i][j] = resp_mat_read[i][j];

    }
  }

  int space_depth = 4;
  int space_breadth = 6;
  int resp_space_len = space_depth * space_breadth;

  float* response_space_matrix = new float[resp_space_len];

  for(i = 0; i < resp_space_len; i++)
  {
    response_space_matrix[i] = resp_space_mat_read[i]
  }

  cout << "ln 73 main." << endl;

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

  cout << "ln 87 main." << endl;

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
    data_mat_len = chs;
    for (int dat_ind = 0; dat_ind < data_mat_len; dat_ind++)
    {
      source_decon[dat_ind] = data_sim[dat_ind];
    }
  }

  for(int cs = 0; cs < chs; cs++)
  {
    source_decon[cs] = source_decon[cs] / (2.0 * 3600);
  }

  int num_iter = 10000;
  int num_rep = 5;
  double boost = 10;

  // !!
  //
  // NORMALIZED RESP BY ITENTSITY (TOTAL COUNTS) AND SOURCE BY TIME (2*3600)
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

  // Specify fineness/coarseness of spatial deconvolution (3 most coarse, 9 most fine)
  int fine = 3;
  int space_iter = 50;

  spatial_decon spaceOut(response_space_matrix,
                         resp_space_len,
                         source_space,
                         fine,
                         iso_count,
                         space_iter
                         );

  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

  std::cout << "time: " << duration << '\n';

  for(int i = 0; i < num_spectra; ++i)
  {
    delete [] response_matrix[i];
  }

  delete [] response_matrix;

  delete [] response_space_matrix;

  return 0;
}
