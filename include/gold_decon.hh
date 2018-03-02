#ifndef gold_decon_hh
#define gold_decon_hh

using namespace std;

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
   	double lda, ldb, ldc, area;

};

#endif
