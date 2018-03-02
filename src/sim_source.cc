#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "sim_source.hh"

sim_source::sim_source()
{

  cout << "begin sim source." << endl;

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

  cout << "end sim source." << endl;


}

sim_source::~sim_source() {}
