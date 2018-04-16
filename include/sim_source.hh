#ifndef sim_source_hh
#define sim_source_hh

using namespace std;

extern std::vector<float> data_sim;
extern std::vector< std::vector<float> > space_sim;
extern std::vector< std::vector<int> > mesh_tracker;

class sim_source
{

  public:
    sim_source(float*,
               int,
               int
               );
    ~sim_source();

  private:

    int sim_size_y, sim_size_x, randy, randx, source_count;
    float x, max;

    std::vector<int> degrees, mesh_temp;
    std::vector<float> temp;

};

#endif
