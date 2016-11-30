#include <igl/readOBJ.h>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <igl/readSTL.h>
#include <igl/viewer/Viewer.h>
#include "tutorial_shared_path.h"
#include <igl/writeOFF.h> 
#include <igl/writePLY.h>

using namespace Eigen;  
using namespace std;
using namespace igl;

struct Mesh
{
  Eigen::MatrixXd V,U,N,UV; // really, N, UV are very useless @
  Eigen::MatrixXi T,F;
} ahead,behind,scene;


//void convertFromObjToPly(std::string);

int main(int argc, char *argv[])
{
  //if(!readSTL("horseBehind.stl",behind.V,behind.F,behind.N))
  if(!readSTL("horseBehind.stl",behind.V,behind.F,behind.N))
  {
    cout<<"failed to load horse Behind stl "<<endl;
  } 
  if(!readSTL("horseAhead.stl",behind.V,behind.F,behind.N))
  {
    cout<<"failed to load horse Ahead stl "<<endl;
  }
  /***********************************************************/ 
  // SETUP LibIgl Viewer 
  /***********************************************************/ 
  igl::viewer::Viewer viewer;
  viewer.data.set_mesh(behind.V, behind.F);
  viewer.launch();
}
