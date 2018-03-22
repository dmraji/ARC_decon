#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "sim_source.hh"

sim_source::sim_source()
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

  cout << "end sim source." << endl;

  space_sim = [
               5650,  2030,  1010,  201,  161,  332,
               2101,  1409,  610,   148,  321,  590,
               956,   532,   147,   215,  599,  946,
               368,   189,   130,   104,  298,  513,
               134,   96,    72,    65,   112,  304
               ];

}

sim_source::~sim_source() {}
