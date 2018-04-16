#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include <stdlib.h>
#include <stdio.h>

#include "sim_source.hh"

sim_source::sim_source(float *resp_space,
                       int resp_space_len,
                       int source_count
                       )
{

  cout << "begin sim source." << endl;

  // Assembly of spectral sim (reading in)
  ifstream inFile;

  inFile.open("source_sim.txt");
  if (!inFile)
  {
    cout << "Unable to open file" << endl;
    exit(1); // terminate with error
  }

  while (inFile >> x)
  {
    data_sim.push_back(x);
  }

  inFile.close();

  // Assembly of spatial sim
  sim_size_y = 80;
  sim_size_x = 100;

  float surv_sim[sim_size_y][sim_size_x] = {};

  // for (int i = 0; i < sim_size_y; i++)
  // {
  //    for (int j = 0; j < sim_size_x; j++)
  //    {
  //        cout << surv_sim[i][j] << " ";
  //    }
  //    std::cout << '\n';
  // }

  // Set up gradient deducer
  for(int d = 0; d < 9; d++)
  {
    degrees.push_back(d);
  }

  // source_count = rand() % 5 + 1;
  // std::cout << source_count << '\n';

  for(int source_ctr = 0; source_ctr < source_count; source_ctr++)
  {
    randy = rand() % sim_size_y;
    randx = rand() % sim_size_x;

    mesh_temp.push_back(randy);
    mesh_temp.push_back(randx);

    mesh_tracker.push_back(mesh_temp);

    mesh_temp.clear();

    surv_sim[randy][randx] = surv_sim[randy][randx] + resp_space[0];
    // std::cout << "resp: " << resp_space[0] << '\n';

    // Assembly of gradient
    for(int deg_ind_y = 0; deg_ind_y < degrees.size(); deg_ind_y++)
    {
      for(int deg_ind_x = 0; deg_ind_x < degrees.size(); deg_ind_x++)
      {
        if(deg_ind_y == 0 && deg_ind_x == 0)
        {
          // std::cout << "seed" << '\n';
        }
        else
        {
          if(randy + degrees[deg_ind_y] < sim_size_y && randx + degrees[deg_ind_x] < sim_size_x)
          {
            // std::cout << "yplus xplus" << '\n';
            surv_sim[randy + degrees[deg_ind_y]][randx + degrees[deg_ind_x]] =
              surv_sim[randy + degrees[deg_ind_y]][randx + degrees[deg_ind_x]] +
              resp_space[(deg_ind_x + 1) + (deg_ind_y * 9) - 1];
            // std::cout << randy + degrees[deg_ind_y] << ", " << randx + degrees[deg_ind_x] << '\n';
            // std::cout << "resp: " << resp_space[(deg_ind_x + 1) + (deg_ind_y * 9) - 1] << '\n';
            // if(surv_sim[randy + degrees[deg_ind_y]][randx + degrees[deg_ind_x]] < 0)
            // {
            //   std::cout << "y deg: " << degrees[deg_ind_y] << " x deg: " << degrees[deg_ind_x] << '\n';
            //   std::cout << "resp: " << resp_space[(deg_ind_x + 1) + (deg_ind_y * 9) - 1] << '\n';
            // }
          }
          if(deg_ind_y > 0)
          {
            if(randy - degrees[deg_ind_y] >= 0 && randx + degrees[deg_ind_x] < sim_size_x)
            {
              // std::cout << "yminus xplus" << '\n';
              surv_sim[randy - degrees[deg_ind_y]][randx + degrees[deg_ind_x]] =
                surv_sim[randy - degrees[deg_ind_y]][randx + degrees[deg_ind_x]] +
                resp_space[(deg_ind_x + 1) + (deg_ind_y * 9) - 1];
              // std::cout << randy - degrees[deg_ind_y] << ", " << randx + degrees[deg_ind_x] << '\n';
              // std::cout << "resp: " << resp_space[(deg_ind_x + 1) + (deg_ind_y * 9) - 1] << '\n';
              // if(surv_sim[randy - degrees[deg_ind_y]][randx + degrees[deg_ind_x]] < 100)
              // {
              //   std::cout << "y deg: " << degrees[deg_ind_y] << " x deg: " << degrees[deg_ind_x] << '\n';
              //   std::cout << "resp: " << resp_space[(deg_ind_x + 1) + (deg_ind_y * 9) - 1] << '\n';
              // }
            }

          }
          if(deg_ind_x > 0)
          {
            if(randy + degrees[deg_ind_y] < sim_size_y && randx - degrees[deg_ind_x] >= 0)
            {
              // std::cout << "yplus xminus" << '\n';
              surv_sim[randy + degrees[deg_ind_y]][randx - degrees[deg_ind_x]] =
                surv_sim[randy + degrees[deg_ind_y]][randx - degrees[deg_ind_x]] +
                resp_space[(deg_ind_x + 1) + (deg_ind_y * 9) - 1];
              // std::cout << randy + degrees[deg_ind_y] << ", " << randx - degrees[deg_ind_x] << '\n';
              // std::cout << "resp: " << resp_space[(deg_ind_x + 1) + (deg_ind_y * 9) - 1] << '\n';
              // if(surv_sim[randy + degrees[deg_ind_y]][randx - degrees[deg_ind_x]] < 0)
              // {
              //   std::cout << "y deg: " << degrees[deg_ind_y] << " x deg: " << degrees[deg_ind_x] << '\n';
              //   std::cout << "resp: " << resp_space[(deg_ind_x + 1) + (deg_ind_y * 9) - 1] << '\n';
              // }
            }

          }
          if(deg_ind_y > 0 && deg_ind_x > 0)
          {
            if(randy - degrees[deg_ind_y] >= 0 && randx - degrees[deg_ind_x] >= 0)
            {
              // std::cout << "yminus xminus" << '\n';
              surv_sim[randy - degrees[deg_ind_y]][randx - degrees[deg_ind_x]] =
                surv_sim[randy - degrees[deg_ind_y]][randx - degrees[deg_ind_x]] +
                resp_space[(deg_ind_x + 1) + (deg_ind_y * 9) - 1];
              // std::cout << randy - degrees[deg_ind_y] << ", " << randx - degrees[deg_ind_x] << '\n';
              // std::cout << "resp: " << resp_space[(deg_ind_x + 1) + (deg_ind_y * 9) - 1] << '\n';
              // if(surv_sim[randy - degrees[deg_ind_y]][randx - degrees[deg_ind_x]] < 0)
              // {
              //   std::cout << "y deg: " << degrees[deg_ind_y] << " x deg: " << degrees[deg_ind_x] << '\n';
              //   std::cout << "resp: " << resp_space[(deg_ind_x + 1) + (deg_ind_y * 9) - 1] << '\n';
              // }
            }

          }
        }
      }
    }
  }

  // Normalization by max element
  max = 0;

  for(int y_ind = 0; y_ind < sim_size_y; y_ind++)
  {
    for(int x_ind = 0; x_ind < sim_size_x; x_ind++)
    {
      if(surv_sim[y_ind][x_ind] > max)
      {
        max = surv_sim[y_ind][x_ind];
      }
    }
  }

  for(int y_ind = 0; y_ind < sim_size_y; y_ind++)
  {
    temp.clear();
    for(int x_ind = 0; x_ind < sim_size_x; x_ind++)
    {
      surv_sim[y_ind][x_ind] = surv_sim[y_ind][x_ind] / max;
      temp.push_back(surv_sim[y_ind][x_ind]);
    }
    space_sim.push_back(temp);
  }

  // for (int i = 0; i < space_sim.size(); i++)
  // {
  //    for (int j = 0; j < space_sim[i].size(); j++)
  //    {
  //        cout << space_sim[i][j] << " ";
  //    }
  //    std::cout << '\n';
  // }

  cout << "end sim source." << endl;

}

sim_source::~sim_source() {}

/*
    ~~~~~~~~~~~~~~~~~~~

    FOR TESTING AS MAIN

    ~~~~~~~~~~~~~~~~~~~
*/

// #include <string>
// #include <iostream>
// #include <fstream>
// #include <vector>
// #include <cmath>
//
// #include <stdlib.h>
// #include <stdio.h>
// #include <postgresql/libpq-fe.h>
//
// #include "sim_source.hh"
// #include "data_read.hh"
//
// sim_source::sim_source(float *resp_space,
//                        int resp_space_len
//                        )
// {
//
//   cout << "begin sim source." << endl;
//
//   // Assembly of spectral sim (reading in)
//   ifstream inFile;
//
//   inFile.open("source_sim.txt");
//   if (!inFile)
//   {
//     cout << "Unable to open file" << endl;
//     exit(1); // terminate with error
//   }
//
//   while (inFile >> x)
//   {
//     data_sim.push_back(x);
//   }
//
//   inFile.close();
//
//   // Assembly of spatial sim
//   sim_size_y = 10;
//   sim_size_x = 10;
//
//   float surv_sim[sim_size_y][sim_size_x] = {{0}};
//
//   // Set up gradient deducer
//   for(int d = 0; d < 8; d++)
//   {
//     degrees.push_back(d);
//   }
//
//   // source_count = rand() % 5 + 1;
//   source_count = 1;
//   // std::cout << source_count << '\n';
//
//   for(int source_ctr = 0; source_ctr < source_count; source_ctr++)
//   {
//     randy = rand() % sim_size_y;
//     randx = rand() % sim_size_x;
//     surv_sim[randy][randx] = surv_sim[randy][randx] + resp_space[0];
//     // std::cout << "resp: " << resp_space[0] << '\n';
//
//     // Assembly of gradient
//     for(int deg_ind_y = 0; deg_ind_y < degrees.size(); deg_ind_y++)
//     {
//       for(int deg_ind_x = 0; deg_ind_x < degrees.size(); deg_ind_x++)
//       {
//         if(deg_ind_y == 0 && deg_ind_x == 0)
//         {
//           // std::cout << "seed" << '\n';
//         }
//         else
//         {
//           if(randy + degrees[deg_ind_y] < sim_size_y && randx + degrees[deg_ind_x] < sim_size_x)
//           {
//             // std::cout << "yplus xplus" << '\n';
//             surv_sim[randy + degrees[deg_ind_y]][randx + degrees[deg_ind_x]] =
//               surv_sim[randy + degrees[deg_ind_y]][randx + degrees[deg_ind_x]] +
//               resp_space[(deg_ind_x + 1) + (deg_ind_y * 9) - 1];
//             // std::cout << randy + degrees[deg_ind_y] << ", " << randx + degrees[deg_ind_x] << '\n';
//             // std::cout << "resp: " << resp_space[(deg_ind_x + 1) + (deg_ind_y * 9) - 1] << '\n';
//             // if(surv_sim[randy + degrees[deg_ind_y]][randx + degrees[deg_ind_x]] < 0)
//             // {
//             //   std::cout << "y deg: " << degrees[deg_ind_y] << " x deg: " << degrees[deg_ind_x] << '\n';
//             //   std::cout << "resp: " << resp_space[(deg_ind_x + 1) + (deg_ind_y * 9) - 1] << '\n';
//             // }
//           }
//           if(deg_ind_y > 0)
//           {
//             if(randy - degrees[deg_ind_y] >= 0 && randx + degrees[deg_ind_x] < sim_size_x)
//             {
//               // std::cout << "yminus xplus" << '\n';
//               surv_sim[randy - degrees[deg_ind_y]][randx + degrees[deg_ind_x]] =
//                 surv_sim[randy - degrees[deg_ind_y]][randx + degrees[deg_ind_x]] +
//                 resp_space[(deg_ind_x + 1) + (deg_ind_y * 9) - 1];
//               // std::cout << randy - degrees[deg_ind_y] << ", " << randx + degrees[deg_ind_x] << '\n';
//               // std::cout << "resp: " << resp_space[(deg_ind_x + 1) + (deg_ind_y * 9) - 1] << '\n';
//               // if(surv_sim[randy - degrees[deg_ind_y]][randx + degrees[deg_ind_x]] < 100)
//               // {
//               //   std::cout << "y deg: " << degrees[deg_ind_y] << " x deg: " << degrees[deg_ind_x] << '\n';
//               //   std::cout << "resp: " << resp_space[(deg_ind_x + 1) + (deg_ind_y * 9) - 1] << '\n';
//               // }
//             }
//
//           }
//           if(deg_ind_x > 0)
//           {
//             if(randy + degrees[deg_ind_y] < sim_size_y && randx - degrees[deg_ind_x] >= 0)
//             {
//               // std::cout << "yplus xminus" << '\n';
//               surv_sim[randy + degrees[deg_ind_y]][randx - degrees[deg_ind_x]] =
//                 surv_sim[randy + degrees[deg_ind_y]][randx - degrees[deg_ind_x]] +
//                 resp_space[(deg_ind_x + 1) + (deg_ind_y * 9) - 1];
//               // std::cout << randy + degrees[deg_ind_y] << ", " << randx - degrees[deg_ind_x] << '\n';
//               // std::cout << "resp: " << resp_space[(deg_ind_x + 1) + (deg_ind_y * 9) - 1] << '\n';
//               // if(surv_sim[randy + degrees[deg_ind_y]][randx - degrees[deg_ind_x]] < 0)
//               // {
//               //   std::cout << "y deg: " << degrees[deg_ind_y] << " x deg: " << degrees[deg_ind_x] << '\n';
//               //   std::cout << "resp: " << resp_space[(deg_ind_x + 1) + (deg_ind_y * 9) - 1] << '\n';
//               // }
//             }
//
//           }
//           if(deg_ind_y > 0 && deg_ind_x > 0)
//           {
//             if(randy - degrees[deg_ind_y] >= 0 && randx - degrees[deg_ind_x] >= 0)
//             {
//               // std::cout << "yminus xminus" << '\n';
//               surv_sim[randy - degrees[deg_ind_y]][randx - degrees[deg_ind_x]] =
//                 surv_sim[randy - degrees[deg_ind_y]][randx - degrees[deg_ind_x]] +
//                 resp_space[(deg_ind_x + 1) + (deg_ind_y * 9) - 1];
//               // std::cout << randy - degrees[deg_ind_y] << ", " << randx - degrees[deg_ind_x] << '\n';
//               // std::cout << "resp: " << resp_space[(deg_ind_x + 1) + (deg_ind_y * 9) - 1] << '\n';
//               // if(surv_sim[randy - degrees[deg_ind_y]][randx - degrees[deg_ind_x]] < 0)
//               // {
//               //   std::cout << "y deg: " << degrees[deg_ind_y] << " x deg: " << degrees[deg_ind_x] << '\n';
//               //   std::cout << "resp: " << resp_space[(deg_ind_x + 1) + (deg_ind_y * 9) - 1] << '\n';
//               // }
//             }
//
//           }
//         }
//       }
//     }
//   }
//
//   // Normalization by max element
//   max = 0;
//
//   for(int y_ind = 0; y_ind < sim_size_y; y_ind++)
//   {
//     for(int x_ind = 0; x_ind < sim_size_x; x_ind++)
//     {
//       if(surv_sim[y_ind][x_ind] > max)
//       {
//         max = surv_sim[y_ind][x_ind];
//       }
//       if(surv_sim[y_ind][x_ind] > 25000 || surv_sim[y_ind][x_ind] < 10)
//       {
//         surv_sim[y_ind][x_ind] = 0;
//       }
//     }
//   }
//
//   for(int y_ind = 0; y_ind < sim_size_y; y_ind++)
//   {
//     temp.clear();
//     for(int x_ind = 0; x_ind < sim_size_x; x_ind++)
//     {
//       // surv_sim[y_ind][x_ind] = surv_sim[y_ind][x_ind] / max;
//       temp.push_back(surv_sim[y_ind][x_ind]);
//     }
//     space_sim.push_back(temp);
//   }
//
//   for (int i = 0; i < space_sim.size(); i++)
//   {
//      for (int j = 0; j < space_sim[i].size(); j++)
//      {
//          cout << space_sim[i][j] << " ";
//      }
//      std::cout << '\n';
//   }
//
//   cout << "end sim source." << endl;
//
//   //
//   //
//   // For kyle and dem boyz
//   //
//   //
//
//   cout << "begin data read." << endl;
//
//   PGconn   *conn;
//   PGresult *res;
//
//   // Look at Michael's shim code for connection over socket
//   conn = PQconnectdb("dbname=sim_space_db user=postgres password=postgres hostaddr=127.0.0.1 port=5432");
//
//   if (PQstatus(conn) == CONNECTION_BAD) {
//          puts("We were unable to connect to the database");
//          exit(0);
//   }
//
//   res = PQexec(conn, "delete from det_sim_space");
//
//   // res = PQexec(conn,
//   //        "create table det_sim_space (counts integer NOT NULL, neutron smallint NOT NULL, gamma smallint NOT NULL, xbin integer NOT NULL, ybin integer NOT NULL)");
//   for(int i = 0; i < sim_size_y; i++)
//   {
//     for(int j = 0; j < sim_size_x; j++)
//     {
//       float rand_src = ((double) rand() / (RAND_MAX));
//       std::cout << rand_src << '\n';
//       string big_str = "insert into det_sim_space (ind, counts, neutron, gamma, xbin, ybin) values ("+to_string(j + (sim_size_x * (i)))+", "+to_string(space_sim[i][j])+", "+to_string(round(space_sim[i][j] * rand_src))+", "+
//           to_string(round(space_sim[i][j] * (1 - rand_src)))+", "+to_string(j)+", "+to_string(i)+")";
//       res = PQexec(conn, big_str.c_str());
//     }
//   }
//
//
//   res = PQexec(conn,
//          "select * from det_sim_space order by ind");
//
//   if (PQresultStatus(res) != PGRES_TUPLES_OK) {
//          puts("We did not get any data!");
//          exit(0);
//   }
//
//   int rec_count = PQntuples(res);
//   // cout << rec_count << '\n';
//
//   printf("We received %d records.\n", rec_count);
//   puts("==========================");
//
//   for (int row=0; row<rec_count; row++) {
//       for (int col=0; col<6; col++) {
//           printf("%s\t", PQgetvalue(res, row, col));
//       }
//       puts("");
//   }
//
//   for (int r=0; r<rec_count; r++)
//   {
//       std::cout << "values= ";
//
//       for (int c=0; c<5; c++)
//       {
//           char* value = PQgetvalue(res, r, c);
//           float fval = atof(value);
//
//           std::cout << fval << " ";
//
//           // Push the int into the temporary vector<float>
//           temp.push_back(fval);
//       }
//       std::cout << std::endl;
//       data_mat.push_back(temp);
//   }
//
//   // for (int i = 0; i < data_mat.size(); i++)
//   // {
//   //     for (int j = 0; j < data_mat[i].size(); j++)
//   //     {
//   //         cout << data_mat[i][j] << endl;
//   //     }
//   // }
//
//   data_mat_len = data_mat[8].size();
//
//   puts("==========================");
//
//   PQclear(res);
//
//   PQfinish(conn);
//
//   cout << "end data read." << endl;
//   exit(1);
//
// }
//
// sim_source::~sim_source() {}
