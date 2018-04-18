// Reading in radiation detector output files
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <stdlib.h>
#include <stdio.h>
#include <postgresql/libpq-fe.h>

#include "data_read.hh"

using namespace std;

data_read::data_read()
{

  cout << "begin data read." << endl;

  PGconn   *conn;
  PGresult *res;

  // Look at Michael's shim code for connection over socket
  conn = PQconnectdb("dbname=det_db user=postgres password=postgres hostaddr = 127.0.0.1 port = 5432");

  if (PQstatus(conn) == CONNECTION_BAD) {
         puts("We were unable to connect to the database");
         exit(0);
  }

  res = PQexec(conn,
         "update det_out set psd=\'0.0005\' where ind=3");

  res = PQexec(conn,
         "select * from det_out order by ind");

  if (PQresultStatus(res) != PGRES_TUPLES_OK) {
         puts("We did not get any data!");
         exit(0);
  }

  rec_count = PQntuples(res);
  // cout << rec_count << '\n';

  printf("We received %d records.\n", rec_count);
  puts("==========================");

  for (row=0; row<rec_count; row++) {
      for (col=0; col<8; col++) {
          printf("%s\t", PQgetvalue(res, row, col));
      }
      puts("");
  }

  for (int r=0; r<rec_count; r++)
  {
      std::cout << "values= ";

      for (int c=0; c<8; c++)
      {
          char* value = PQgetvalue(res, r, c);
          float fval = atof(value);

          std::cout << fval << " ";

          // Push the int into the temporary vector<float>
          temp.push_back(fval);
      }
      std::cout << std::endl;
      data_mat.push_back(temp);
  }

  for (int i = 0; i < data_mat.size(); i++)
  {
      for (int j = 0; j < data_mat[i].size(); j++)
      {
          cout << data_mat[i][j] << endl;
      }
  }

  data_mat_len = data_mat[0].size();

  puts("==========================");

  PQclear(res);

  PQfinish(conn);

  cout << "end data read." << endl;

}

// object destructor

data_read::~data_read() {}
