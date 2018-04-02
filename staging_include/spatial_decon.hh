#ifndef spatial_decon_hh
#define spatial_decon_hh

using namespace std;

extern std::vector< std::vector<int> > end_ind;

class spatial_decon {

  public:

    spatial_decon(float*,
                  int,
                  float**,
                  int,
                  int,
                  int,
                  int,
                  int,
                  int
                  );
    ~spatial_decon();

  private:

    int resp_index, combos, deg_temp, layer_save, finer, isos, iso_iter, err, layer_seek, savey, savex;
    float att_cps, macro_max, rand_source, min_dif, max_ele, dif_sum;
    std::vector<int> degrees;

    std::vector<int> y_end;
    std::vector<int> x_end;

};

#endif
