#include <string>
#include <iostream>
#include <vector>
#include <array>

#include "gold_decon.hh"

using namespace std;

histo::histo(float *source_matrix)
{

  source_decon[] = {0};

  for(int p=0; p<source_matrix.size()-1; p++)
  {
    source_decon[source_matrix[p]] = source_decon[source_matrix[p]] + 1;
  }

  return source_decon;
}

~histo::histo() {}
