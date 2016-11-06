/*
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

Eigen::MatrixXd V_two;
Eigen::MatrixXi F_two;

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
    viewer.data.set_mesh(V_two,F_two);
    viewer.core.align_camera_center(V_two,F_two);
  } 
  else if ( key == '3' ) 
  {
    
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
2 Switch to rotated view 1
3 Switch to mean curvature based view ( convert to rotated, as best as possible ) 
    )";

  igl::readOFF(TUTORIAL_SHARED_PATH "/bunny.off", V_one, F_one);
  igl::readOFF(TUTORIAL_SHARED_PATH "pointCloud2.off", V_two, F_two); // needs a better name ( idenitity mesh, goal mesh ) 
*/

  

  /***********************************************************/ 
  // TODO : set up as a preallocation step 
  // CALCULATE the rotated form of the mesh 
  /***********************************************************/ 
/*
  Eigen::Vector3d e1 ( 0, 0, -1 ); 
  Eigen::Vector3d e2 ( 1, 0, 0 ); 
  Eigen::Matrix3d rot_matrix_zAxis = igl::rotation_matrix_from_directions (e1, e2 );
  V_two = V_one * rot_matrix_zAxis;
  F_two = F_one; 
*/

  /***********************************************************/ 
  // STEP (1) :: SELECT control points pi \in P ( i = 1..N), , compute surface normals n_pi, set initial transformation matrix T_0
  /***********************************************************/ 

  /* TODO : IS it here that I need to perform barycentric based random sampling ? */
/*
   Eigen::MatrixXd control_points = V_one;
   int numControlPoints = control_points.rows();

   // I am very unsure? Plus there is the issue of getting normals too , from barycentric ( not sure about that too ! ) 

  Eigen::MatrixXd transforMat_k = Eigen::MatrixXd::Identity(numRowsDataSet1, 3);
  igl::per_vertex_normals(V_one, F_one, meshOneVertexNormals);
*/

/*
  int k;
  int maxIters = 100;
  // so mean curvature is determined nby # of basis functions ... and I believe we need "n" basis functions for "n" verticse. NONETHELESS ... I SHOULD TAKE A LOOK @ THIS AGAIN ! 
  int numVerticesMeshOne = V_one.rows();
  Eigen::SparseMatrix<float> stiffnessMatrix(numVerticesMeshOne,numVerticesMeshOne); 
  Eigen::SparseMatrix<float> massMatrix(numVerticesMeshOne,numVerticesMeshOne); 

  igl::cotmatrix(V_one,F_one, sitffnessMatrix );  
  //MassMatrixType mcfType = MASSMATRIX_TYPE_BARYCENTRIC;
  igl::cotmatrix(V_one,F_one, MASSMATRIX_TYPE_BARYCENTRIC, sitffnessMatrix );  

 /// shit ... how to get the set of initial coefficients  vecotr(x(0))??  any libigl utility to solve this, and hat basis too ? 
  // Run for timesteps \delta , for max_iter number of iterations
  for ( k = 1; k <= maxIters; k += 1 ) 
  {
      
       







  }

  /***********************************************************/ 
  // Plot the initial mesh 
  /***********************************************************/ 
  /*igl::viewer::Viewer viewer;
  //viewer.callback_key_down = &key_down;
  viewer.data.set_mesh(V_one, F_one);
  viewer.launch();
}
*/
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

















