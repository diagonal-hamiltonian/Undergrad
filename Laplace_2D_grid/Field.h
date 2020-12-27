#This class acts as a 2D grid with elements, phi(i,j), which can be accessed and replaced through the brackets operator (.. , ..)

#include <vector>
class Field
{
  private:
    int             qx_;    // This is nx+1, the number of grid
                    // points in the x direction when
                    // the range is sliced into nx parts
    int             qy_;
    std::vector<double>        data_;

    int index_(int x, int y) const { return x + qx_*y; }
    
  public:
    Field(int nx, int ny) : qx_(nx+1), qy_(ny+1), data_(qx_*qy_)
    {
      int n=qx_*qy_;
      for (int i=0;i<n;i++) data_[i] = 0.0;
    }

    double  operator() (int x, int y) const { return data_[index_(x,y)]; }
    double& operator() (int x, int y)       { return data_[index_(x,y)]; }
    
};

