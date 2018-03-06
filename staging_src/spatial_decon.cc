#include <string>
#include <iostream>
#include <cmath>
#include <vector>

#include "spatial_decon.hh"

using namespace std;

spatial_decon::spatial_decon(// array of number of mesh points with counts as values
                             )
{
  // This var is here for show, will be read in from main in reality

  // Describes 25 mesh points
  int* mesh_counts = new int[
   0,  1,  2,  3,  4,
   5,  6,  7,  8,  9,
   10, 11, 12, 13, 14,
   15, 16, 17, 18, 19,
   20, 21, 22, 23, 24
  ]

  /*  Sample sub-matrices from source matrix

      Size of sub-mat depends on sampling coarse or fine, minimum size of 2x2

      Search against 3D response matrix, normalized via count rate

  */

  for 

}

spatial_decon::~spatial_decon() {}
