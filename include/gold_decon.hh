#ifndef gold_decon_hh
#define gold_decon_hh

using namespace std;

extern int iso_count;

class gold_decon
{

  public:
    gold_decon(float**,
               float*,
               int,
               int,
               int,
               int,
               double
               );
    ~gold_decon();

  private:

    int i, j, k, lindex, lhx = 0, repet;
    int wksp_len;
   	double lda, ldb, ldc, area;
    std::vector<double> lda_vec;
    std::vector<string> iso_names;

};

#endif
