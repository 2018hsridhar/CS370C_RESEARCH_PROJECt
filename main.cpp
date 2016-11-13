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

Eigen::SparseMatrix<double> massMatrix_iterK ;
Eigen::SparseMatrix<double> stiffnessMatrix_iterK ; 
double delta = 0.001; 
int k = 0;

void applyOneTimeStepOfMeanCurvatureFlow();
bool key_down( igl::viewer::Viewer& , unsigned char , int );

// function is called when keyboard buttons are pressed down. useful for alternating amongst a set of differing views
bool key_down( igl::viewer::Viewer& viewer, unsigned char key, int modifier)
{
  std::cout << "Key : " << key << (unsigned int) key << std::endl;
  if ( key == '1' ) {
    // clear data before drawing mesh
    viewer.data.clear();
    viewer.data.set_mesh(V_one,F_one);
    viewer.core.align_camera_center(V_one,F_one); 
  }
  else if ( key == '2' ) {
    // clear data before drawing mesh
    viewer.data.clear();
	applyOneTimeStepOfMeanCurvatureFlow();
    viewer.data.set_mesh(V_mcf,F_mcf);
    viewer.core.align_camera_center(V_mcf,F_mcf);
  } 
  else if ( key == '3' ) {
    // clear data before drawing mesh
    viewer.data.clear();
    V_mcf = V_one;
    F_mcf = F_one;
    k = 0;
    viewer.data.set_mesh(V_mcf,F_mcf);
    viewer.core.align_camera_center(V_mcf,F_mcf);
  } 
  return false;
}

void applyOneTimeStepOfMeanCurvatureFlow()
{
	Eigen::SparseMatrix<double> A = ( massMatrix_iterK - ( delta * stiffnessMatrix_iterK ));
	Eigen::MatrixXd B = ( massMatrix_iterK * V_mcf); 
	std::cout << " [1] Caluclated (M-delta*L) and (M*v) matrices \n";
	// Q1 :: should I be solving for an EXACT sol ( direct method ) or APPROX sol ( iterative method ) ?? NOT SURE !! SEE NOTES for when to tell there is or is not an exact solution ! 
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > solver; 
	// choice was made, since , ACCORDING to EIGEN, this was the msot basic sparse-matrix solver 
	// PLUS cotan matrix is self-adjoint ( i believe. .. need to check ). other matrix properties do not fit here !
	solver.compute(A); 
	if(solver.info() != Eigen::Success) {
		std::cout << "Decomposition of A failed." << std::endl;
	} 
	auto updatedMeshVertices = solver.solve(B);
	if ( solver.info() != Eigen::Success )  {
		std::cout << "Solving B failed." << std::endl;
	}
	V_mcf = updatedMeshVertices.eval(); // what is the difference between "solve" and "eval" ???
	std::cout << " [2] Passed solver tests. \n";

	// update vertices, laplacian matrix, and mass matrix
	igl::cotmatrix(V_mcf,F_mcf, stiffnessMatrix_iterK );   
	/* #TODO :: find out what is causing the error here
	* TECHNICALLy, the stiffness matrix has to be calcualted ONLY once ! Although this too, should work !
	* IN ADDITION, the mass matrix is also experiencing issues ... but this is only l8r !
	*/
	igl::MassMatrixType mcfType = igl::MASSMATRIX_TYPE_BARYCENTRIC;
	igl::massmatrix(V_mcf,F_mcf, mcfType, massMatrix_iterK);  
	std::cout << "[3] Succesfully updated mass and stiffness matrices \n";
    k += 1;
	std::cout << " [4] Finished iteration "<< (k-1) << "." <<  std::endl;
}

int main(int argc, char *argv[])
{


  // *********************************************************** 
  // PRINT OUT CRITICAL INFORMATION TO USER / DEVELOPER 
  // LOAD mesh data ( OFF format )
  // Q1. what asssumptions can I make about my mesh input data?? I'm pretty sure I cannot assume that they are the same size !
  // *********************************************************** 

  std::cout << R"(
1 Switch to initial view
2 Run mean curvature based view 
3 Reset mean curvature based view  ( can rerun again ) 
    )";

//  igl::readOFF(TUTORIAL_SHARED_PATH "/bunny.off", V_one, F_one); 
  igl::readOFF(TUTORIAL_SHARED_PATH "/cow.off", V_one, F_one); 
  V_mcf = V_one;
  F_mcf = F_one;

  /***********************************************************/ 
  // #TODO :: FILL IN THIS ENTRY LATER 
  /***********************************************************/ 

  igl::cotmatrix(V_one,F_one, stiffnessMatrix_iterK );  
  igl::MassMatrixType mcfType = igl::MASSMATRIX_TYPE_BARYCENTRIC;
  igl::massmatrix(V_one,F_one, mcfType, massMatrix_iterK);  

  /***********************************************************/ 
  // PLOT initial current mesh 
  /***********************************************************/ 
  igl::viewer::Viewer viewer;
  viewer.callback_key_down = &key_down;
  viewer.data.set_mesh(V_one, F_one);
  viewer.launch();
  return 0;
}
