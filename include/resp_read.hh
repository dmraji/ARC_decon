#ifndef resp_read_hh
#define resp_read_hh

using namespace std;

class resp_read {

  public:

    resp_read();
    ~resp_read();

    extern std::vector<std::vector <float> > resp_mat_read;

  private:

    std::vector <string> spectra_names;

    std::vector <float> temp;

    int x;

};

#endif
