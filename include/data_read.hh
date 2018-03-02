#ifndef data_read_hh
#define data_read_hh

using namespace std;

extern std::vector<std::vector <float> > data_mat;
extern int data_mat_len;

class data_read
{

  public:

    data_read();
    ~data_read();

  private:

    vector<float> temp;

    int rec_count;
    int row;
    int col;

    float x;

};

#endif
