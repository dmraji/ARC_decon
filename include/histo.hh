#ifndef histo_hh
#define histo_hh

using namespace std;

extern float source_decon[4096];

class histo
{

  public:
    histo(float*,
          int);
    ~histo();

  private:
    int temp_var;

};

#endif
