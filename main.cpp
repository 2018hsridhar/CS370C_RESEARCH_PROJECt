#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include "tutorial_shared_path.h"
#include <igl/rotation_matrix_from_directions.h>
#include <igl/writeOFF.h> // #include <igl/per_vertex_normals.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/procrustes.h>

using namespace Eigen; 
using namespace std;

Eigen::MatrixXd V_one;
Eigen::MatrixXi F_one;

Eigen::MatrixXd V_two;
Eigen::MatrixXi F_two;

Eigen::MatrixXd V_approx;
Eigen::MatrixXi F_approx;

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
    // clear data before drawing mesh
    viewer.data.clear();
    viewer.data.set_mesh(V_approx,F_approx);
    viewer.core.align_camera_center(V_approx,F_approx);
  }
  return false;
}

int main(int argc, char *argv[])
{

  /***********************************************************/ 
  // PRINT out critical information to user / developer 
  // LOAD mesh data , in OFF format
  /***********************************************************/ 
  std::cout << R"(
1 switch to identity view
2 Switch to rotated view
3 Switch to ICP view 
    )";

  if(!igl::readOFF(TUTORIAL_SHARED_PATH "/bunny.off", V_one, F_one))
  {
    cout << " Failed to load mesh " << endl;
  }
  if(!igl::readOFF(TUTORIAL_SHARED_PATH "/dataset2.off", V_two, F_two)) // needs a better name ( idenitity mesh, goal mesh ) 
  {
    cout << " Failed to load mesh " << endl;
  }


  /***********************************************************/ 
  // PREALLOCATE matrices for computation ( p_i, q, transformation_iter_k)
  /***********************************************************/ 
  
  /* TODO : IS it here that I need to perform barycentric based random sampling ? */
  Eigen::MatrixXd p_i = convertToHomogenousForm(V_one); // #TODO check if need for sampling from here
  Eigen::MatrixXd q = convertToHomogenousForm(V_two); // #TODO check if need for sampling from here

  int numControlPoints = p_i.rows();
  Eigen::MatrixXd transformMat_iterK = Eigen::MatrixXd::Identity(4, 4); 
  Eigen::VectorXd hack = Eigen::Vector4d (0,0,0,0);

   /*
	std::cout << "Analysis of transformMat_iterK ( rotate ) " << "(" << transformMat_iterK.rows() << "," << transformMat_iterK.cols() << ")" << std::endl;
	std::cout << "Analysis of hack ( translate ) " << "(" << hack.rows() << "," << hack.cols() << ")" << std::endl;
	std::cout << transformMat_iterK << std::endl;
	std::cout << hack << std::endl;
   */

  Eigen::MatrixXd transformMatAppliedTo_p_i_iterK = (p_i * transformMat_iterK).rowwise() + hack.transpose(); // #TODO :: why error here? seems to be a dumb typing thing TBH


  //Eigen::MatrixXd transformMatAppliedTo_p_i_iterK = p_i; //#TODO :: get previous like working correctly? weird type error occuring
  // #TODO :: understand how construction MatrixXd XPrime = (X*R).rowwise() + t.transpose() makes sense

  /***********************************************************/ 
  // ITERATIVELY improve rigid ICP method , for solving the transformation matrix
  /***********************************************************/ 
  int k;
  int maxIters = 10; // so this blew up @ 100 iterations ... was something incorrect here. #TODO assert correctness !
  for ( k = 1; k <= maxIters; k += 1 ) 
  {
    // we need to create a new matrix of minimum q's !!
    // take a row of current transformation matrix applied to  p_i, find closest q in mesh-two
    // I need to make sure that I am passing in the correct option !
    //  #TODO :: understand why "Ele" is needed in closest-points code ! 
    // SADLY, this does not accomadate Homogenous coordinates

    /***********************************************************/ 
    // DISCOVER closest points in q, from transformMat_iter_k * p_i
    /***********************************************************/ 
    Eigen::VectorXi Ele = Eigen::VectorXi::LinSpaced(q.rows(),0,q.rows() - 1);
    Eigen::VectorXd smallestSquaredDists;
    Eigen::VectorXi smallestDistIndxs;
    Eigen::MatrixXd q_j_k; 				// this coresponds to closestPointsTo_p_i_From_q ;
    igl::point_mesh_squared_distance(normalizeHomogenousMatrix(transformMatAppliedTo_p_i_iterK),normalizeHomogenousMatrix(q),Ele,smallestSquaredDists,smallestDistIndxs,q_j_k);  
    Eigen::MatrixXd q_j_k_homog = convertToHomogenousForm(q_j_k);

    /***********************************************************/ 
    // APPLY Procrstues to solve for T, that minimizes norm (T*p_i - q_j_k )
    // UPDATE transformation matrix via T_k = T * T_{k-1} 
    /***********************************************************/ 
    Eigen::Matrix4d Rotate;
    Eigen::Vector4d Translate;
    double Scale;
    igl::procrustes(transformMatAppliedTo_p_i_iterK, q_j_k_homog, false,false,Scale,Rotate,Translate ); 

	/*
		std::cout << "Analysis of rotate    " << "(" << Rotate.rows() << "," << Rotate.cols() << ")" << std::endl;
		std::cout << "Analysis of translate " << "(" << Translate.rows() << "," << Translate.cols() << ")" << std::endl;
		std::cout << Rotate << std::endl;
		std::cout << Translate << std::endl;
	*/

    // #TODO :: assert the correctness of the update step  
    Eigen::Matrix4d newTransformMatrix = Rotate;
    newTransformMatrix.col(3) += Translate;  
	transformMat_iterK = newTransformMatrix * transformMat_iterK;  
  }

  /***********************************************************/ 
  // CALCULATE the mesh based on RIGID ICP transformation 
  /***********************************************************/ 
  F_approx = F_one;
  Eigen::MatrixXd rigidIcpMesh = (p_i * transformMat_iterK).rowwise() + hack.transpose(); 
  V_approx = normalizeHomogenousMatrix(rigidIcpMesh);

  /***********************************************************/ 
  // SETUP LibIgl Viewer 
  /***********************************************************/ 
  igl::viewer::Viewer viewer;
  viewer.callback_key_down = &key_down;
  viewer.data.set_mesh(V_one, F_one);
  viewer.launch();
}

/***********************************************************/ 
//APPROACH 2 ( using surface normals ... code this up later ) 
//Eigen::MatrixXd meshOneVertexNormals_prime = transformMat_k * meshOneVertexNormals;
// find intersection qi_k on surface Q, based on normal lines !
// STEP (1) :: SELECT control points pi \in P ( i = 1..N), , compute surface normals n_pi, set initial transformation matrix T_0
// I am very unsure? Plus there is the issue of getting normals too , from barycentric ( not sure about that too ! ) 
/***********************************************************/ 

