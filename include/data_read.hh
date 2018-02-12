#ifndef data_read_hh
#define data_read_hh

using namespace std;

class data_read() {

  public:
    data_read();
    ~data_read();

    vector<vector<float> > read_mat(0, vector<float>(rec_count, 1));

};

#endif
