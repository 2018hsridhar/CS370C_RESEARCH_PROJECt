// IN mainc.pp file for TUTORIAL 102
#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include "tutorial_shared_path.h"
#include <igl/rotation_matrix_from_directions.h>
#include <igl/writeOFF.h>
#include <igl/per_vertex_normals.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h> // remember !! use barycentric ! 


using namespace Eigen; 
using namespace std;

Eigen::MatrixXd V_one;
Eigen::MatrixXi F_one;

Eigen::MatrixXd V_mcf;
Eigen::MatrixXi F_mcf;

// note :: this is ONLY FOR POINTS ... if they were vectors ... then I need a different method !
Eigen::MatrixXd convertToHomogenousForm(const Ref<const MatrixXd>& mat )
{
  Eigen::MatrixXd homogenoizedMatrix;
  homogenoizedMatrix = mat.transpose().colwise().homogeneous().transpose(); 
  return homogenoizedMatrix;
}

Eigen::MatrixXd normalizeHomogenousMatrix (const Ref<const MatrixXd>& mat )
{
  Eigen::MatrixXd normalizedMatrix;
  normalizedMatrix = mat.transpose().colwise().hnormalized().transpose(); 
  return normalizedMatrix;
}

// function is called when keyboard buttons are pressed down. useful for alternating amongst a set of differing views
bool key_down( igl::viewer::Viewer& viewer, unsigned char key, int modifier)
{
  std::cout << "Key : " << key << (unsigned int) key << std::endl;
  if ( key == '1' )
  {
    // clear data before drawing mesh
    viewer.data.clear();
    viewer.data.set_mesh(V_one,F_one);
    viewer.core.align_camera_center(V_one,F_one); 
  }
  else if ( key == '2' ) 
  {
    // clear data before drawing mesh
    viewer.data.clear();
    viewer.data.set_mesh(V_mcf,F_mcf);
    viewer.core.align_camera_center(V_mcf,F_mcf);
  } 
  return false;
}

int main(int argc, char *argv[])
{


  // *********************************************************** 
  // PRINT OUT CRITICAL INFORMATION TO USER / DEVELOPER 
  // Load mesh data , in OFF format
  // Q1. what asssumptions can I make about my mesh input data?? I'm pretty sure I cannot assume that they are the same size !
  // *********************************************************** 

  std::cout << R"(
1 switch to initial view
2 Switch to mean curvature based view ( convert to rotated, as best as possible ) 
    )";

  igl::readOFF(TUTORIAL_SHARED_PATH "/bunny.off", V_one, F_one);


  /***********************************************************/ 
  // #TODO :: FILL IN THIS ENTRY LATER 
  /***********************************************************/ 

  int numVerticesMeshOne = V_one.rows();
  Eigen::SparseMatrix<float> stiffnessMatrix (numVerticesMeshOne,numVerticesMeshOne); 
  Eigen::SparseMatrix<float> massMatrix (numVerticesMeshOne,numVerticesMeshOne); 

  igl::cotmatrix(V_one,F_one, sitffnessMatrix );  
  MassMatrixType mcfType = MASSMATRIX_TYPE_BARYCENTRIC;
  igl::massmatrix(V_one,F_one, mcfType , massMatrix);  

	// split columns
	// solve the system
	// put_back_values 
	// update mass-stiffness matrices    

  // Run for timesteps \delta , for max_iter number of iterations
  int k;
  int maxIters = 100;
  double delta = 0.1; 
  V_mcf = V_one;
  F_mcf = F_one;
  for ( k = 1; k <= maxIters; k += 1 ) 
  {
      for ( m = 0; m < 3; m++ )
      {
          MatrixXd vertex_dimComponent = V_mcf.col(m);        
//          solve(); #TODO :: this should be a method , when I think about it ! 

      }
  }

  /***********************************************************/ 
  // Plot the initial mesh 
  /***********************************************************/ 
  igl::viewer::Viewer viewer;
  viewer.callback_key_down = &key_down;
  viewer.data.set_mesh(V_one, F_one);
  viewer.launch();
}

  /*// note :: do store a couple of test matrices ... they seem to be useful for l8r purposes :-) 
  Eigen::Matrix3f m = Matrix3f::Random();
  std::cout << m << std::endl;

  Eigen::MatrixXf homogenous = convertToHomogenousForm(m);
  std::cout << "\n\n\n" << std::endl;
  std::cout << homogenous << std::endl;

  Eigen::MatrixXf deHomogenous = normalizeHomogenousMatrix(homogenous);
  std::cout << "\n\n\n" << std::endl;
  std::cout << deHomogenous << std::endl;
*/

















