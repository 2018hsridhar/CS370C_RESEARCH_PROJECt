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
#include <igl/is_border_vertex.h>
#include <set>

#include "meanCurvatureFlow.h"

/*
using namespace Eigen; 
using namespace std;

Eigen::MatrixXd V_one;
Eigen::MatrixXi F_one;

Eigen::MatrixXd V_mcf;

Eigen::SparseMatrix<double> massMatrix_iterK ;
Eigen::SparseMatrix<double> stiffnessMatrix_iterK ; 
double delta = 0.001;
int k = 0;
*/

/*
void applyOneTimeStepOfMcfBoundaryCase();
void applyOneTimeStepOfMcfWatertightCase();
bool key_down( igl::viewer::Viewer& , unsigned char , int );
bool meshHasBoundary();
void applyMeanCurvatureFlow();
*/

bool meshHasBoundary(Ref<MatrixXd>& V, Ref<MatrixXi>& F)
{
	unsigned int boundaryCount = 0;	
    std::vector<bool> boundaryVerticesStatus = igl::is_border_vertex(V, F);  
	for (auto i: boundaryVerticesStatus) {
        if ( i == true)
            boundaryCount++;
    }
    return ( boundaryCount > 0 );
}

/*
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
	//applyOneTimeStepOfMcfWatertightCase();
	applyOneTimeStepOfMcfBoundaryCase();
    viewer.data.set_mesh(V_mcf,F_one);
    viewer.core.align_camera_center(V_mcf,F_one);
    std::string output_of_mesh = "meanCurvaureFlowOutput[" + std::to_string(k) + "].stl";
    std::cout << "WRITING MCF data for iteration [" << std::to_string(k) << " ].\n"; 
    igl::writeOFF(output_of_mesh, V_mcf, F_one);
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

void applyOneTimeStepOfMcfBoundaryCase()
{
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

		// update vertices ( coefficent vector ) 
        //std::vector<int> boundaryVerticesIndexes = igl::boundary_loop(F_one); 
       // [0] solve for boundary vertices of object, and rewrite those @ end. ignore them in curent analysis !
		// NOTE :: we still keep same dim for Mass-Stiffness matrices; we'll just update coeffs l8r 

 /// let us quicklyu assert that "boundaryEdges-edges" actually makes sense !

		Eigen::MatrixXi boundaryEdges;
		Eigen::MatrixXi edges;
*/
/* ... DIAGNOSING ...
        igl::exterior_edges(F_one,boundaryEdges);
        igl::edges(F_one, edges);

        std::cout << "num boundary edges = " << boundaryEdges.rows() << std::endl;
        std::cout << "total number edges = " << edges.rows() << std::endl;
*/
/*
        int boundaryCount = 0;	
        std::vector<bool> boundaryVerticesStatus = igl::is_border_vertex(V_mcf, F_one);  
		for (auto i: boundaryVerticesStatus) {
            if ( i == true)
                boundaryCount++;
        }
        //std::cout << " num boundary vertices = " << boundaryCount << std::endl;
        //std::cout << " total num vertices = " << V_mcf.rows() << std::endl;

        Eigen::MatrixXd oldVertices = V_mcf;
        for(int i = 0; i < boundaryVerticesStatus.size(); i++) {
            if ( boundaryVerticesStatus[i] ) 
                V_mcf.row(i) = oldVertices.row(i);
            else
                V_mcf.row(i) = newVertices.row(i); 
        }
        //std::cout << [3] Succesfully updated vertices \n";

 		// update mass matrix
		igl::MassMatrixType mcfType = igl::MASSMATRIX_TYPE_BARYCENTRIC;
		igl::massmatrix(V_mcf,F_one, mcfType, massMatrix_iterK);  
		//std::cout << "[3] Succesfully updated mass matrix \n";

		k += 1;
		//std::cout << " [5] Finished iteration "<< (k-1) << "." <<  std::endl;
}

void applyOneTimeStepOfMcfWatertightCase()
{
		Eigen::SparseMatrix<double> A = ( massMatrix_iterK - ( delta * stiffnessMatrix_iterK ));
		Eigen::MatrixXd B = ( massMatrix_iterK * V_mcf); 
	//	std::cout << " [1] Caluclated (M-delta*L) and (M*v) matrices \n";
		Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > solver; 
		solver.compute(A); 
		if(solver.info() != Eigen::Success) {
			std::cout << "Decomposition of A failed." << std::endl;
		} 
		auto updatedMeshVertices = solver.solve(B);
		if ( solver.info() != Eigen::Success )  {
			std::cout << "Solving B failed." << std::endl;
		}
		auto newVertices = updatedMeshVertices.eval(); 
		if ( solver.info() != Eigen::Success )  {
			std::cout << "Solving B (actually) failed." << std::endl;
		}
        
	    //	std::cout << " [2] Passed solver tests. \n";

		// update vertices ( coefficent vector ) and mass matrix
		Eigen::VectorXd dbla;
		igl::doublearea(V_mcf,F_one,dbla);
		double oldVolume = 0.5 * dbla.sum(); 

        V_mcf = newVertices;
		igl::MassMatrixType mcfType = igl::MASSMATRIX_TYPE_BARYCENTRIC;
		igl::massmatrix(V_mcf,F_one, mcfType, massMatrix_iterK);  
		//std::cout << "[3] Succesfully updated mass and stiffness matrices \n";

		Eigen::VectorXd dbla_new;
		igl::doublearea(V_mcf,F_one,dbla_new);
		double newVolume = 0.5 * dbla_new.sum(); 

		V_mcf /= std::cbrt(newVolume); 

		//std::cout << "[4] Rescaled by mesh volume / unit area \n"; 
		k += 1;
		//std::cout << " [5] Finished iteration "<< (k-1) << "." <<  std::endl;
}

int main(int argc, char *argv[])
{

  std::cout << R"(
	1 Switch to initial view
	2 Run mean curvature based view 
	3 Reset mean curvature based view  ( can rerun again ) 
    )";

  // LOAD mesh data ( OFF format )
  igl::readOFF(TUTORIAL_SHARED_PATH "/camelhead.off",V_one,F_one);
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
void applyMeanCurvatureFlow()
{
  // CALCUALTE initial mass and stiffness matrices
  igl::cotmatrix(V_one,F_one, stiffnessMatrix_iterK );  
  igl::MassMatrixType mcfType = igl::MASSMATRIX_TYPE_BARYCENTRIC;
  igl::massmatrix(V_one,F_one, mcfType, massMatrix_iterK);  

  // #TODO :: like Kazhdan, I need to identify numerical instabilities ( terminate this algorithm early ) 
  bool haveBndryVertices= meshHasBoundary();
  if(haveBndryVertices) {
  	for ( int i = 0; i < 512; i++) {
		applyOneTimeStepOfMcfBoundaryCase();
  	}
  }
  else { 
  	for ( int i = 0; i < 512; i++) {
		applyOneTimeStepOfMcfWatertightCase();
  	}
  }

}

*/

// STL caauses issues ... OBJ and OFF formats do not ... not sure why though !
  //igl::readOFF(TUTORIAL_SHARED_PATH "/cow.off", V_one, F_one); 
  //igl::readOFF(TUTORIAL_SHARED_PATH "/proper_sphere.off", V_one, F_one);  // it straight up does not work for this case !
  //igl::readSTL("decimated-max.stl", V_one, F_one,N_one);  // it straight up does not work for this case !
  //igl::readOBJ(TUTORIAL_SHARED_PATH "/decimated-max.obj",V_one,F_one);
  //igl::readOFF("horseAhead.off", V_one, F_one);  // it straight up does not work for this case !
  //igl::readOFF(TUTORIAL_SHARED_PATH "/beetle.off",V_one,F_one); /// this finally works 
