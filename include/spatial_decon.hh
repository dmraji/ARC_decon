#ifndef spatial_decon_hh
#define spatial_decon_hh

using namespace std;

class spatial_decon {

  public:

    spatial_decon(float*,
                  int,
                  float**,
                  int,
                  int,
                  int,
                  int
                  );
    ~spatial_decon();

  private:

    int resp_index, combos, deg_temp, layer_save, finer, isos, iso_iter;
    float att_cps;
    std::vector<int> degrees;

};

#endif
