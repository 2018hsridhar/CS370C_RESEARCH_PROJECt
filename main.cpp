// SADLY, volume works ONLY for tetrahedral meshes :-( . so too , is face_areas.h !
// NOT sure if I should use vector_area.h or double_area.h?? ... double area ! vector_area only inhputs FACES ( gives you a generic matrix really ) ... no vertex info thoguh !
#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/viewer/Viewer.h>
#include "tutorial_shared_path.h"
#include <igl/rotation_matrix_from_directions.h>
#include <igl/writeOFF.h>
#include <igl/per_vertex_normals.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h> 
#include <igl/doublearea.h>

using namespace Eigen; 
using namespace std;

Eigen::MatrixXd V_one;
Eigen::MatrixXi F_one;

Eigen::MatrixXd V_mcf;

Eigen::SparseMatrix<double> massMatrix_iterK ;
Eigen::SparseMatrix<double> stiffnessMatrix_iterK ; 
//double delta = 0.00001; clearluy, a smaller delta = more stability, but also more computaiton
double delta = 0.001;
int k = 0;

void applyOneTimeStepOfMeanCurvatureFlow();
bool key_down( igl::viewer::Viewer& , unsigned char , int );
bool convertObjToOff( std::string );

bool convertObjToOff( std::string fileOfInterest )
{
    igl::readOBJ(TUTORIAL_SHARED_PATH "/sphere.obj", V_mcf, F_one);
    igl::writeOFF("proper_sphere.off",V_mcf,F_one); 
}

bool key_down( igl::viewer::Viewer& viewer, unsigned char key, int modifier)
{
  std::cout << "Key : " << key << (unsigned int) key << std::endl;
  if ( key == '1' ) {
    viewer.data.clear();
    viewer.data.set_mesh(V_one,F_one);
    viewer.core.align_camera_center(V_one,F_one); 
  }
  else if ( key == '2' ) {
    viewer.data.clear();
	applyOneTimeStepOfMeanCurvatureFlow();
    viewer.data.set_mesh(V_mcf,F_one);
    viewer.core.align_camera_center(V_mcf,F_one);
  } 
  else if ( key == '3' ) {
    viewer.data.clear();
    V_mcf = V_one;
    k = 0;
    igl::cotmatrix(V_one,F_one, stiffnessMatrix_iterK );  
    igl::MassMatrixType mcfType = igl::MASSMATRIX_TYPE_BARYCENTRIC;
    igl::massmatrix(V_one,F_one, mcfType, massMatrix_iterK);  
    viewer.data.set_mesh(V_one,F_one);
    viewer.core.align_camera_center(V_mcf,F_one);
  } 
  return false;
}

void applyOneTimeStepOfMeanCurvatureFlow()
{
    k = 0;
		Eigen::SparseMatrix<double> A = ( massMatrix_iterK - ( delta * stiffnessMatrix_iterK ));
		Eigen::MatrixXd B = ( massMatrix_iterK * V_mcf); 
	//	std::cout << " [1] Caluclated (M-delta*L) and (M*v) matrices \n";
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
		auto newVertices = updatedMeshVertices.eval(); // what is the difference between "solve" and "eval" ???
		if ( solver.info() != Eigen::Success )  {
			std::cout << "Solving B (actually) failed." << std::endl;
		}
        
	    //	std::cout << " [2] Passed solver tests. \n";

		// update vertices, laplacian matrix, and mass matrix
		//igl::cotmatrix(V_mcf,F_one, stiffnessMatrix_iterK );   
		/* #TODO :: find out what is causing the error here
		* TECHNICALLY, the stiffness matrix has to be calcualted ONLY once ! Although this too, should work !
		* IN ADDITION, the mass matrix is also experiencing issues ... but this is only l8r !
		*/
        V_mcf = newVertices;
		igl::MassMatrixType mcfType = igl::MASSMATRIX_TYPE_BARYCENTRIC;
		igl::massmatrix(V_mcf,F_one, mcfType, massMatrix_iterK);  
		//std::cout << "[3] Succesfully updated mass and stiffness matrices \n";

		Eigen::VectorXd dbla;
		igl::doublearea(V_mcf,F_one,dbla);
		double newVolume = 0.5 * dbla.sum(); 
		V_mcf /= std::sqrt(newVolume);

		//std::cout << "[4] Rescaled by mesh volume / unit area \n"; 
		k += 1;
		//std::cout << " [5] Finished iteration "<< (k-1) << "." <<  std::endl;
}

int main(int argc, char *argv[])
{

  // PRINT OUT CRITICAL INFORMATION TO USER / DEVELOPER 
  std::cout << R"(
	1 Switch to initial view
	2 Run mean curvature based view 
	3 Reset mean curvature based view  ( can rerun again ) 
    )";

  // LOAD mesh data ( OFF format )
  //igl::readOFF(TUTORIAL_SHARED_PATH "/bunny.off", V_one, F_one); 
  igl::readOFF(TUTORIAL_SHARED_PATH "/cow.off", V_one, F_one); 
  //igl::readOFF(TUTORIAL_SHARED_PATH "/proper_sphere.off", V_one, F_one);  // it straight up does not work for this case !
  V_mcf = V_one;

  // CALCUALTE initial mass and stiffness matrices
  igl::cotmatrix(V_one,F_one, stiffnessMatrix_iterK );  
  igl::MassMatrixType mcfType = igl::MASSMATRIX_TYPE_BARYCENTRIC;
  igl::massmatrix(V_one,F_one, mcfType, massMatrix_iterK);  

  // PLOT initial mesh 
  igl::viewer::Viewer viewer;
  viewer.callback_key_down = &key_down;
  viewer.data.set_mesh(V_one, F_one);
  viewer.launch();

  return 0;
}







  // *********************************************************** 
  // Q1. what asssumptions can I make about my mesh input data?? I'm pretty sure I cannot assume that they are the same size !
  // *********************************************************** 
