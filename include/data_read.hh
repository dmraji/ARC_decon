#ifndef data_read_hh
#define data_read_hh

using namespace std;

class data_read {

  public:

    data_read();
    ~data_read();

    extern std::vector<std::vector <float> > data_mat;

  private:

    vector<float> temp;

};

#endif
