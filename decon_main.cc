#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <sstream>
#include <iterator>
#include <cmath>
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
vector< vector<float> > data_mat, data_space, data_space_gamma, data_space_starttimes, data_space_stoptimes;
vector <float> times, en, en_gamma, xbin, pos_xbin, ybin, pos_ybin;
vector <int> neut;
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
std::vector<float> conf_vec;

// Source var init
float source_decon[4096];

int iso_count;

int main(int argc, char** argv) {

  label:
    std::cout << "Re-init." << '\n';

  std::cout << "Decon start." << '\n';

  // Allow true random numbers
  srand(time(NULL));

  std::clock_t start;
  double duration;

  start = std::clock();

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

    for(int l = 0; l < data_mat_len; l++)
    {
      times.push_back(data_mat[3][l]);
      en.push_back(data_mat[4][l]);
      xbin.push_back(data_mat[5][l]);
      ybin.push_back(data_mat[6][l]);
      neut.push_back(data_mat[7][l]);

    }

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
  else
  {

    // Parsing data read

    int min_xbin = 0;
    int max_xbin = 0;
    int min_ybin = 0;
    int max_ybin = 0;
    float time_s = 0;
    int xbin_cur = 0;
    int ybin_cur = 0;

    string delim;
    std::vector < string > date_time;
    std::vector < string > hr_min_sec;
    float comp_time;

    // Adjust for negative space values
    pos_xbin.resize(xbin.size());
    pos_ybin.resize(ybin.size());

    for(int s = 0; s < xbin.size(); s++)
    {
      if(xbin[s] < min_xbin)
      {
        min_xbin = xbin[s];
      }
      if(xbin[s] > max_xbin)
      {
        max_xbin = xbin[s];
      }
      if(ybin[s] < min_ybin)
      {
        min_ybin = ybin[s];
      }
      if(ybin[s] > max_ybin)
      {
        max_ybin = ybin[s];
      }
    }
    for(int s = 0; s < xbin.size(); s++)
    {
      pos_xbin[s] = xbin[s] - min_xbin;
      pos_ybin[s] = ybin[s] - min_ybin;
    }

    max_xbin = max_xbin - min_xbin;
    max_ybin = max_ybin - min_ybin;

    data_space.resize(max_xbin);
    data_space_starttimes.resize(max_xbin);
    data_space_stoptimes.resize(max_xbin);
    for(int o = 0; o < max_xbin; o++)
    {
      data_space[o].resize(max_ybin);
      data_space_starttimes[o].resize(max_ybin);
      data_space_stoptimes[o].resize(max_ybin);
    }

    for(int s = 0; s < pos_xbin.size(); s++)
    {

      // Splitting timestamp strings
      istringstream f(to_string(times[s]));
      string str;
      while(getline(f, str, ' '))
      {
        date_time.push_back(str);
      }
      istringstream g(date_time[1]);
      while(getline(g, str, ':'))
      {
        hr_min_sec.push_back(str);
      }

      comp_time = (3600 * stoi(hr_min_sec[0])) + (60 * stoi(hr_min_sec[1])) + stoi(hr_min_sec[2]);

      if(data_space_starttimes[pos_xbin[s]][pos_ybin[s]] == 0)
      {

        data_space_starttimes[pos_xbin[s]][pos_ybin[s]] = comp_time;
      }
      data_space_stoptimes[pos_xbin[s]][pos_ybin[s]] = comp_time;

      data_space[pos_xbin[s]][pos_ybin[s]] = data_space[pos_xbin[s]][pos_ybin[s]] + 1;
      if(neut[s] == 0)
      {
        data_space_gamma[pos_xbin[s]][pos_ybin[s]] = data_space_gamma[pos_xbin[s]][pos_ybin[s]] + 1;
      }

    }

    // Transform to count rates
    for(int s = 0; s < pos_xbin.size(); s++)
    {
      data_space[pos_xbin[s]][pos_ybin[s]] = data_space[pos_xbin[s]][pos_ybin[s]] /
        data_space_stoptimes[pos_xbin[s]][pos_ybin[s]] - data_space_starttimes[pos_xbin[s]][pos_ybin[s]];
    }

  }

  float* source_mat = new float[data_mat_len];

  if(sim == 0)
  {
    for(int en_ct = 0; en_ct < en.size(); en_ct++)
    {
      if(neut[en_ct] == 0)
      {
        en_gamma.push_back(en[en_ct]);
      }
    }
    for(int dat_ind = 0; dat_ind < en_gamma.size(); dat_ind++)
    {
      source_mat[dat_ind] = en_gamma[dat_ind];
    }
    histo histIt(source_mat,
                 en_gamma.size()
                 );
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

  std::vector<string> iso_names;
  iso_names.push_back("Co-57");
  iso_names.push_back("Pu-239");
  iso_names.push_back("Ba-133");
  iso_names.push_back("Cf-252");
  iso_names.push_back("Cs-137");
  iso_names.push_back("Th-232");
  iso_names.push_back("Am-241");
  iso_names.push_back("Natural Uranium");
  iso_names.push_back("U-235");
  iso_names.push_back("Co-60");
  iso_names.push_back("Ra-226");

  gold_decon callIt(response_matrix,
                    source_decon,
                    chs,
                    num_spectra,
                    num_iter,
                    num_rep,
                    boost,
                    iso_names
                    );

  // for(int printer = 0; printer < 4096; printer++)
  // {
  //   cout << source_decon[printer] << endl;
  // }

  std::vector<float> spect_res;
  std::vector<string> iso_res;

  for(int i = 0; i < num_spectra; i++)
  {
    if(source_decon[i] > 0)
    {
      spect_res.push_back(source_decon[i]);
      iso_res.push_back(iso_names[i]);
    }
  }

  // Spatial deconvolution algorithm call

  // Specify fineness/coarseness of spatial deconvolution (3 most coarse, 5 intermediate, 7 most fine)
  //  ONLY USE 3, 5 OR 7
  int fine = 7;
  int space_iter = 10;

  // Turn on for farther "view", better spatial algorithm; off for speed; (warning: off not tested)
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

  /*

      ~~~~~~~~~~~~~~~~~~~~~

      RESULTS & DATA EXPORT

      ~~~~~~~~~~~~~~~~~~~~~

  */

  PGconn   *conn;
  PGresult *res;

  conn = PQconnectdb("dbname=det_db user=postgres password=postgres hostaddr=127.0.0.1 port=5432");

  if(PQstatus(conn) == CONNECTION_BAD)
  {
    puts("We were unable to connect to the database");
    exit(0);
  }

  // res = PQexec(conn, "delete from heatmap");

  // res = PQexec(conn,
  //        "create table det_sim_space (ind SERIAL, counts integer NOT NULL, neutron smallint NOT NULL, gamma smallint NOT NULL, xbin integer NOT NULL, ybin integer NOT NULL)");

  // Heatmap Results
  if(sim == 0)
  {

    res = PQexec(conn, "create table heatmap (ind SERIAL, counts integer NOT NULL, neutron smallint NOT NULL, gamma smallint NOT NULL, xbin integer NOT NULL, ybin integer NOT NULL)");

    for(int i = 0; i < data_space.size(); i++)
    {
      for(int j = 0; j < data_space[0].size(); j++)
      {

        string big_str = "insert into heatmap (ind, counts, neutron, gamma, xbin, ybin) values ("
                           +to_string(j + (data_space[0].size() * (i)))+", "
                           +to_string(data_space[i][j])+", "+to_string(round(data_space[i][j] - data_space_gamma[i][j]))
                           +", "+to_string(data_space_gamma[i][j])+", "+to_string(j)+", "+to_string(i)+")";

        res = PQexec(conn, big_str.c_str());
      }
    }
  }

  //
  // Spectral Decon Results
  //
  res = PQexec(conn, "create table spect_res (ind SERIAL, iso character varying(32) NOT NULL, intensity integer NOT NULL)");

  for(int i = 0; i < spect_res.size(); i++)
  {

      string spect_str = "insert into spect_res (ind, iso, intensity) values ("
                         +to_string(i)+", "+iso_res[i]+", "+to_string(spect_res[i])+")";

      res = PQexec(conn, spect_str.c_str());
  }

  //
  // Spatial Decon Results
  //
  res = PQexec(conn, "create table spatial_res (ind SERIAL, xbin integer NOT NULL, ybin integer NOT NULL, confidence real NOT NULL)");

  for(int i = 0; i < end_ind.size(); i++)
  {

      string spatial_str = "insert into spatial_res (ind, xbin, ybin, confidence) values ("
                         +to_string(i)+", "+to_string(end_ind[i][0])+", "+to_string(end_ind[i][1])
                         +", "+to_string(conf_vec[i])+")";

      res = PQexec(conn, spatial_str.c_str());
  }

  res = PQexec(conn, "select * from heatmap order by ind");

  if(PQresultStatus(res) != PGRES_TUPLES_OK)
  {
    puts("We did not get any data!");
    exit(0);
  }

  int rec_count = PQntuples(res);
  // cout << rec_count << '\n';

  printf("We received %d records.\n", rec_count);
  puts("==========================");

  for(int row = 0; row < rec_count; row++)
  {
    for(int col = 0; col < 6; col++)
    {
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

  // goto label;

  return 0;
}
