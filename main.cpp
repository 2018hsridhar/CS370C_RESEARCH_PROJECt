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

 // igl::readOFF(TUTORIAL_SHARED_PATH "/bunny.off", V_one, F_one); //---> the naive approach is slow for the bunny
  igl::readOFF(TUTORIAL_SHARED_PATH "/cow.off", V_one, F_one); 


  /***********************************************************/ 
  // #TODO :: FILL IN THIS ENTRY LATER 
  /***********************************************************/ 

  Eigen::SparseMatrix<double> massMatrix_iterK ;
  Eigen::SparseMatrix<double> stiffnessMatrix_iterK ; 

//  igl::cotmatrix(V_one,F_one, stiffnessMatrix_iterK );  //  technically, the stiffness matrix is calculated only once !
  igl::MassMatrixType mcfType = igl::MASSMATRIX_TYPE_BARYCENTRIC;
  igl::massmatrix(V_one,F_one, mcfType, massMatrix_iterK);  

	// split columns
	// solve the system
	// put_back_values 
	// update mass-stiffness matrices    
    // Run for timesteps \delta , for max_iter number of iterations

  int k;
  int maxIters = 15; // for some reason ... I get a weird error @ say, 5 iterations ...vs 2 iterations ... WHAT!! ... there should be no error with iterations ! not sure why mass + stiffness matrix is getting updated weirdly
  double delta = 0.001; 
  V_mcf = V_one;
  F_mcf = F_one;
  for ( k = 1; k <= maxIters; k += 1 ) 
  {
      std::cout << "here 1\n";
      // newVertices.setZero(); // see alec jacobson's post for why this was used  .. no need :: just use auto !

      Eigen::SparseMatrix<double> A = ( massMatrix_iterK - ( delta * stiffnessMatrix_iterK ));
      Eigen::MatrixXd B = ( massMatrix_iterK * V_mcf); 
      // Q1 :: should I be solving for an EXACT sol ( direct method ) or APPROX sol ( iterative method ) ?? 
      Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > solver; 
	  // choice was made, since , ACCORDING to EIGEn, this was the msot basic sparse-matrix solver 
 	  // PLUS cotan matrix is self-adjoint ( i believe. .. need to check ). other matrix properties do not fit here !
      std::cout << "here 2\n";
      solver.compute(A); 
      if(solver.info() != Eigen::Success) {
          std::cout << "Decomposition of A failed." << std::endl;
      } 
      auto updatedMeshVertices = solver.solve(B);
      if ( solver.info() != Eigen::Success )  {
          std::cout << "Solving B failed." << std::endl;
      }
      V_mcf = updatedMeshVertices.eval(); // what is the difference between "solve" and "eval" ???
     
      // update vertices, laplacian matrix, and mass matrix
      std::cout << "here 3\n"; // there is an issue here! ( line immediately below ... I wonder why ! ) 
      //igl::cotmatrix(V_mcf,F_mcf, stiffnessMatrix_iterK );   // error is here? I wonder why ??
      igl::MassMatrixType mcfType = igl::MASSMATRIX_TYPE_BARYCENTRIC;
      igl::massmatrix(V_mcf,F_mcf, mcfType, massMatrix_iterK);  
      std::cout << "here 4\n";
  }

  /***********************************************************/ 
  // Plot the initial mesh 
  /***********************************************************/ 
  igl::viewer::Viewer viewer;
  viewer.callback_key_down = &key_down;
  viewer.data.set_mesh(V_one, F_one);
  viewer.launch();
}
