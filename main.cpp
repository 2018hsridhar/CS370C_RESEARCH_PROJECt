#include <igl/readSTL.h>
#include <igl/readOFF.h>
#include <igl/readOBJ.h>

#include <igl/writeSTL.h>
#include <igl/viewer/Viewer.h>
#include "tutorial_shared_path.h"
#include <igl/rotation_matrix_from_directions.h>
#include <igl/writeOFF.h>
#include <igl/per_vertex_normals.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/doublearea.h>
#include <igl/exterior_edges.h> 
#include <igl/edges.h> 
#include <igl/is_boundary_edge.h>
#include <igl/on_boundary.h> 
#include <igl/is_border_vertex.h>
#include <set>

// my own libraries ( note :: test libraries seperately, themselves, THEN include them ) 
#include "meanCurvatureFlow.h" 

using namespace Eigen; 
using namespace std;

Eigen::MatrixXd V;
Eigen::MatrixXi F;

int main(int argc, char *argv[])
{

  // LOAD mesh data ( OFF format )
  igl::readOFF(TUTORIAL_SHARED_PATH "/camelhead.off",V,F);
  //bool hasBndry = meshHasBoundary(V,F);
  //std::cout << hasBndry << std::endl;

  // PLOT initial mesh 
  igl::viewer::Viewer viewer;
  
  //viewer.callback_key_down = &key_down;
  viewer.data.set_mesh(V, F);
  viewer.launch();

  return 0;
}
