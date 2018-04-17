#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <cstdio>

#include <postgresql/libpq-fe.h>

#include "data_read.hh"
#include "gold_decon.hh"
#include "resp_read.hh"
#include "histo.hh"
#include "sim_source.hh"
#include "spatial_decon.hh"

using namespace std;

// Data var init
vector< vector<float> > data_mat;
int data_mat_len;

// Sim var init
vector<float> data_sim;
vector< vector<float> > space_sim;
vector< vector<int> > mesh_tracker;
int source_num;

// Resp var init
vector< vector<float> > resp_mat_read;
vector<float> resp_space_mat_read;

// Spatial data init
int space_datax, space_datay;
vector< vector<int> > end_ind;

// Source var init
float source_decon[4096];

int iso_count;

int main(int argc, char** argv) {

  // Allow true random numbers
  srand(time(NULL));

  std::clock_t start;
  double duration;

  start = std::clock();

  // cout << "hi.\n";

  // Build response matrix from experimental spectra
  int chs = 4096;
  int num_spectra = 11;

  // Call to function to read response files
  resp_read respIt;

  // initialize 2D dynamic array for responses
  float** response_matrix = new float*[num_spectra];
  for(int i = 0; i < num_spectra; i++)
  {
    response_matrix[i] = new float[chs];
  }

  std::cout << "ln 53 main" << '\n';

  for(int j = 0; j < chs; j++)
  {
    for(int i = 0; i < num_spectra; i++)
    {
      response_matrix[i][j] = resp_mat_read[i][j];


    }
  }

  std::cout << "ln 62 main" << '\n';

  // UPDATE THIS
  int space_depth = 9;
  int space_breadth = 9;
  int resp_space_len = space_depth * space_breadth;

  float* response_space_matrix = new float[resp_space_len];

  for(int i = 0; i < resp_space_len; i++)
  {
    response_space_matrix[i] = resp_space_mat_read[i];
    // cout << response_space_matrix[i] << " ";
  }

  cout << "ln 89 main." << endl;

  int sourcex;
  int sourcey;

  int space_lenx;
  int space_leny;

  std::cout << "ln 100 main." << '\n';

  // change this to stop using sim source
  int sim = 1;
  if(sim == 1)
  {
    source_num = 2;

    sim_source simIt(response_space_matrix,
                     resp_space_len,
                     source_num
                     );

    // for (int i = 0; i < space_sim.size(); i++)
    // {
    //    for (int j = 0; j < space_sim[i].size(); j++)
    //    {
    //        cout << space_sim[i][j] << " ";
    //    }
    //    std::cout << '\n';
    // }

    space_lenx = space_sim.size();
    space_leny = space_sim[0].size();

  }
  else
  {
    // Obtain vector of vectors via data read-in
    data_read readIt;

    space_lenx = space_datax;
    space_leny = space_datay;
  }

  // Define var to be called with decon class
  float** source_space = new float*[space_lenx];
  for(int i = 0; i < space_lenx; i++)
  {
    source_space[i] = new float[space_leny];
  }

  if(sim == 1)
  {
    for(int spacey = 0; spacey < space_sim.size(); spacey++)
    {
      for(int spacex = 0; spacex < space_sim[spacey].size(); spacex++)
      {
        source_space[spacey][spacex] = space_sim[spacey][spacex];
      }
    }
  }

  cout << "ln 97 main." << endl;

  float* source_mat = new float[data_mat_len];

  if(sim == 0)
  {
    for (int dat_ind = 0; dat_ind < data_mat_len; dat_ind++)
    {
      source_mat[dat_ind] = data_mat[dat_ind][8];
    }
    histo histIt(source_mat,
                   data_mat_len);
  }
  else
  {
    data_mat_len = chs;
    for (int dat_ind = 0; dat_ind < data_mat_len; dat_ind++)
    {
      source_decon[dat_ind] = data_sim[dat_ind];
    }

    // ADD LINES FOR source_space SIM READ
  }

  for(int cs = 0; cs < chs; cs++)
  {
    source_decon[cs] = source_decon[cs] / (1.0 * 1);
  }

  int num_iter = 10000;
  int num_rep = 5;
  double boost = 10;

  // Energy deconvolution algorithm call

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

  // Spatial deconvolution algorithm call

  // Specify fineness/coarseness of spatial deconvolution (3 most coarse, 5 intermediate, 7 most fine)
  //  ONLY USE 3, 5 OR 7
  int fine = 7;
  int space_iter = 10;

  // Turn on for farther "view", better spatial algorithm; off for speed;
  int future_sight = 1;

  spatial_decon spaceIt(response_space_matrix,
                         resp_space_len,
                         source_space,
                         space_lenx,
                         space_leny,
                         fine,
                         iso_count,
                         space_iter,
                         future_sight,
                         sim,
                         source_num
                         );

  for(int i = 0; i < end_ind.size(); i++)
  {
    std::cout << "predicted space: " << end_ind[i][0] << ", " << end_ind[i][1] << '\n';
  }

  if(sim == 1)
  {
    for(int i = 0; i < source_num; i++)
    {
      std::cout << "real space: " << mesh_tracker[i][0] << ", " << mesh_tracker[i][1] << '\n';
    }
  }

  // Results & Data Export
  PGconn   *conn;
  PGresult *res;

  // Look at Michael's shim code for connection over socket
  conn = PQconnectdb("dbname=det_db user=postgres password=postgres hostaddr=127.0.0.1 port=5432");

  if (PQstatus(conn) == CONNECTION_BAD) {
         puts("We were unable to connect to the database");
         exit(0);
  }

  res = PQexec(conn, "delete from det_sim_space");

  // res = PQexec(conn,
  //        "create table det_sim_space (counts integer NOT NULL, neutron smallint NOT NULL, gamma smallint NOT NULL, xbin integer NOT NULL, ybin integer NOT NULL)");
  for(int i = 0; i < sim_size_y; i++)
  {
    for(int j = 0; j < sim_size_x; j++)
    {
      float rand_src = ((double) rand() / (RAND_MAX));
      std::cout << rand_src << '\n';
      string big_str = "insert into det_sim_space (ind, counts, neutron, gamma, xbin, ybin) values ("+to_string(j + (sim_size_x * (i)))+", "+to_string(space_sim[i][j])+", "+to_string(round(space_sim[i][j] * rand_src))+", "+
          to_string(round(space_sim[i][j] * (1 - rand_src)))+", "+to_string(j)+", "+to_string(i)+")";
      res = PQexec(conn, big_str.c_str());
    }
  }


  res = PQexec(conn,
         "select * from det_sim_space order by ind");

  if (PQresultStatus(res) != PGRES_TUPLES_OK) {
         puts("We did not get any data!");
         exit(0);
  }

  int rec_count = PQntuples(res);
  // cout << rec_count << '\n';

  printf("We received %d records.\n", rec_count);
  puts("==========================");

  for (int row=0; row<rec_count; row++) {
      for (int col=0; col<6; col++) {
          printf("%s\t", PQgetvalue(res, row, col));
      }
      puts("");
  }

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
