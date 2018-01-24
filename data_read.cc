// Reading in radiation detector output files

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

using namespace std;

int main() {

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

  cout << "\ndone\n";
  return 0;
}
