// Reading in radiation detector output files
#include <iostream>
#include <string>
#include <vector>

#include <stdlib.h>
#include <stdio.h>
#include <postgresql/libpq-fe.h>

using namespace std;

data_read::data_read() {

  PGconn   *conn;
  PGresult *res;
  int      rec_count;
  int      row;
  int      col;

  conn = PQconnectdb("dbname=det_test_db user=postgres password=postgres hostaddr = 127.0.0.1 port = 5432");

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
      for (col=0; col<10; col++) {
          printf("%s\t", PQgetvalue(res, row, col));
      }
      puts("");
  }

  vector<vector<float> > read_mat(10, vector<float>(rec_count, 1));
  for (int r=0; r<rec_count; r++)
  {
      std::cout << "values= ";
      vector<float> temp;
      for (int c=1; c<10; c++)
      {
          char* value = PQgetvalue(res, r, c);
          float fval = atof(value);

          std::cout << fval << " ";
          temp.push_back(fval);     // Push the int into the temporary vector<int>
      }
      std::cout << std::endl;
      read_mat.push_back(temp);
  }

  puts("==========================");

  PQclear(res);

  PQfinish(conn);

  return 0;

}

// object destructor

data_read::~data_read() {}
