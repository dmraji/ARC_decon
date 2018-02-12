#include <string>
#include <iostream>
#include <vector>

#include "gold_decon.hh"

using namespace std;

int main() {

  cout << "hi.\n";

  // Obtain vector of vectors via data read-in
  data_mat = data_read();

  // Take the energy values from the data and build a response matrix to put
  //   through Gold's energy deconvolution algorithm
  vector<float> resp_mat;
  for (int dat_ind = 0; dat_ind < data_mat.size(); dat_ind++)
  {
    resp_mat.push_back(data_mat[dat_ind][8]);
  }

  energy_decon_data_mat = gold_decon(new float *source,
                                     resp_mat,
                                     int sizex,
                                     int sizey,
                                     int numberIterations,
                                     int numberRepetitions,
                                     double boost
                                     );

  return 0;
}
