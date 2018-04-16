#ifndef resp_read_hh
#define resp_read_hh

using namespace std;

extern std::vector<std::vector <float> > resp_mat_read;
extern std::vector <float> resp_space_mat_read;

class resp_read {

  public:

    resp_read();
    ~resp_read();

  private:

    std::vector <string> spectra_names;
    std::vector <string> space_names;

    std::vector <float> temp;
    std::vector<float> eff;
    std::vector<float> bg;
    std::vector<float> bg2;

    int x, count_time, start_k;
    float spect_sum;
    string space_n;

};

#endif
