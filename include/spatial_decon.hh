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
                  int,
                  int,
                  int
                  );
    ~spatial_decon();

  private:

    int resp_index, combos, combo_piece, it_prod, deg_temp, layer_save,
        finer, isos, iso_iter, err, layer_seek, survy_save, survx_save,
        savey, savex, real_ct, iso_found, passer, no_rep;
    float att_cps, macro_max, rand_source, min_dif, max_ele, dif_sum,
          resp_sum, conf, max_ele_repo;

    std::vector<int> degrees;
    std::vector<int> y_end, x_end, survskipx, survskipy;
    std::vector<float> min_min_vec, conf_levels;
    std::vector< std::vector<int> > compan_vec;

};

#endif
