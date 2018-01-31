// Reading in radiation detector output files
#include <pqxx/pqxx>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>

using namespace std;
using namespace pqxx;

data_read::data_read(int argc, char** argv) {

  ifstream in_file;

  in_file.open("test.txt");

  // Check that file actually exists
  if (!in_file) {
    cout << "File open failed.\n";
    exit(1);
  }

  std::vector<int>counts_by_ch;
  int num;

  while(in_file >> num) {
    counts_by_ch.push_back(num);
  }

  for (std::vector<int>::const_iterator i = counts_by_ch.begin(); i != counts_by_ch.end(); ++i) {
    std::cout << *i << ' ';
  }

  in_file.close();

  cout << "\nRead successful.\n";
  return 0;
}

// object destructor

data_read::~data_read() {}
