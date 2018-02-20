#ifndef gold_decon_hh
#define gold_decon_hh

using namespace std;

class gold_decon {

  public:
    gold_decon(std::vector< std::vector <float> > rM,
               std::vector <float> s,
               int sx,
               int sy,
               int nI,
               int nR,
               double b);
    ~gold_decon();

    int sizex, sizey;
    const int numberIterations, numberRepetitions;
    const double boost;

    int i, j, k, lindex, lhx = 0, repet;
   	double lda, ldb, ldc, area;

  // Member initialization list
  private:
    std::vector< std::vector <float> > respMatrix(rM);
    std::vector <float> source(s);
    const int ssizex(sx);
    const int ssizey(sy);
    const int numberIterations(nI);
    const int numberRepetitions(nR);
    const double boost(b);

};

#endif
