/*
// IN mainc.pp file for TUTORIAL 102
#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include "tutorial_shared_path.h"
#include <igl/rotation_matrix_from_directions.h>
#include <igl/writeOFF.h>
#include <igl/per_vertex_normals.h>

Eigen::MatrixXd V_one;
Eigen::MatrixXi F_one;

Eigen::MatrixXd V_two;
Eigen::MatrixXi F_two;

Eigen::MatrixXd V_Approx;
Eigen::MatrixXi F_Approx;

Eigen::MatrixXd initMeshVertexNormals;

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


  // PRINT OUT CRITICAL INFORMATION TO USER / DEVELOPER 
  std::cout << R"(
1 switch to identity view
2 Switch to rotated view 1
3 Switch to ICP view 
    )";

  // Load a mesh in OFF format
  // Q1. what asssumptions can I make about my mesh input data?? I'm pretty sure I cannot assume that they are the same size !
  igl::readOFF(TUTORIAL_SHARED_PATH "/bunny.off", V_one, F_one);
  igl::readOFF(TUTORIAL_SHARED_PATH "pointCloud2.off", V_two, F_two); // needs a better name ( idenitity mesh, goal mesh ) 

  int numRowsDataSet1 = V_one.rows();
  int numColsDataSet1 = 3;

  // needs to be set to a PREALLOCATION STEP 
  Eigen::Vector3d e1 ( 0, 0, -1 ); 
  Eigen::Vector3d e2 ( 1, 0, 0 ); 
  Eigen::Matrix3d rot_matrix_zAxis = igl::rotation_matrix_from_directions (e1, e2 );
  V_two = V_one * rot_matrix_zAxis;
  F_two = F_one; 

 

  // step (1) :: select control points pi \in P ( i = 1..N), , compute surface normals n_pi, set initial transformation matrix T_0

  // note :: is it here that I need to perform barycentric based random sampling?? not sure?? there is the issue of getting normals from this too !!

  Eigen::MatrixXd T_0 = Eigen::MatrixXd::Identity(numRowsDataSet1, 3);
  igl::per_vertex_normals(V_one, F_one, initMeshVertexNormals);

  int numIters;
  int maxIters = 1000000;
  for ( numIters = 0; numIters <= maxIters; numIters += 1 ) // note :: isn't there a nice, C++ way of writing a loop ??
  {
  






  }

  // Plot the ROTATED mesh 
  //igl::viewer::Viewer viewer;
  //viewer.callback_key_down = &key_down;
  //viewer.data.set_mesh(V, F);
  //viewer.launch();


}
*/
