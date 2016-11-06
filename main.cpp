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

Eigen::MatrixXd V_Approx;
Eigen::MatrixXi F_Approx;

Eigen::MatrixXd meshOneVertexNormals;

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

int main(int argc, char *argv[])
{

  /***********************************************************/ 
  // PRINT OUT CRITICAL INFORMATION TO USER / DEVELOPER 
  // Load mesh data , in OFF format
  // Q1. what asssumptions can I make about my mesh input data?? I'm pretty sure I cannot assume that they are the same size !
  /***********************************************************/ 
  std::cout << R"(
1 switch to identity view
2 Switch to rotated view 1
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

  // *********************************************************** 
  // TODO : set up as a preallocation step  ( wait ... what is this ?? ) 
  // CALCULATE the rotated form of the mesh 
  // *********************************************************** 

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
   // I am very unsure? Plus there is the issue of getting normals too , from barycentric ( not sure about that too ! ) 

  Eigen::MatrixXd p_i = convertToHomogenousForm(V_one); // #TODO check if need for sampling from here
  Eigen::MatrixXd q = convertToHomogenousForm(V_two); // #TODO check if need for sampling from here
  //Eigen::MatrixXd p_i = V_one; // #TODO check if need for sampling from here
  //Eigen::MatrixXd q = V_two; // #TODO check if need for sampling from here

/* ... can be deleted later ..
  Eigen::MatrixXd Rotate;
  Eigen::VectorXd Translate;
  double Scale;
  std::cout << "Here1\n";
  igl::procrustes(p_i, q, false,false,Scale,Rotate,Translate ); 
  Eigen::MatrixXd applyProc = (p_i * Rotate).rowwise() + Translate.transpose();
  std::cout << "Applied procrustes" << std::endl;
  std::cout << "Analysis of matrix p_i" << "(" << p_i.rows() << "," << p_i.cols() << ")" << std::endl;
  std::cout << "Analysis of rotate    " << "(" << Rotate.rows() << "," << Rotate.cols() << ")" << std::endl;
  std::cout << "Analysis of translate " << "(" << Translate.rows() << "," << Translate.cols() << ")" << std::endl;
  std::cout << "Analysis of applyProc " << "(" << applyProc.rows() << "," << applyProc.cols() << ")" << std::endl;
  std::cout << "Analyze dimensions of all matrices" << std::endl;
*/

  int numControlPoints = p_i.rows();
  Eigen::MatrixXd transformMat_iterK = Eigen::MatrixXd::Identity(3, 3); 
  //Eigen::MatrixXd transformMatAppliedTo_p_i_iterK = (p_i * transformMat_iterK).rowwise(); // err here WHY???
  Eigen::MatrixXd transformMatAppliedTo_p_i_iterK = p_i; //#TODO :: get previous like working correctly? weird type error occuring
  // #TODO :: understand how construction MatrixXd XPrime = (X*R).rowwise() + t.transpose() makes sense

  // used in the computation of closest distance points 
  Eigen::VectorXd smallestSquaredDists;
  Eigen::VectorXi smallestDistIndxs;
  Eigen::MatrixXd q_j_k; 				// really, this coresponds to closestPointsTo_p_i_From_q ;

  int k;
  int maxIters = 10;
  for ( k = 1; k <= maxIters; k += 1 ) 
  {
    int j; 

    // we need to create a new matrix of minimum q's !!
    // take a row of current transformation matrix applied to  p_i, find closest q in mesh-two
    // I need to make sure that I am passing in the correct option !
    Eigen::VectorXi Ele = Eigen::VectorXi::LinSpaced(q.rows(),0,q.rows() - 1);
    //  #TODO :: understand why "Ele" is needed in closest-points code ! 
    // Sadly, this does not accomadate Homogenous coordinates
    igl::point_mesh_squared_distance(normalizeHomogenousMatrix(transformMatAppliedTo_p_i_iterK),normalizeHomogenousMatrix(q),Ele,smallestSquaredDists,smallestDistIndxs,q_j_k);  
    Eigen::MatrixXd q_j_k_homog = convertToHomogenousForm(q_j_k);

    // apply Procrstues to solve for T, that minimizes norm (T*p_i - q_j_k )
    Eigen::MatrixXd Rotate;
    Eigen::VectorXd Translate;
    double Scale;
    igl::procrustes(transformMatAppliedTo_p_i_iterK, q_j_k_homog, false,false,Scale,Rotate,Translate ); 
    // #TODO :: should "r" be allowed to be a reflection matrix? not sure what this exactly means ??
    Eigen::MatrixXd newTransformMatrix = Rotate + Translate.transpose();  // this should work ...I think
    // #TODO :: chck if this caluation is correct. ( line above me ) 
	transformMat_iterK = newTransformMatrix * transformMat_iterK;  
  }

	//approach 2 ( code this up l8r ) 
	//Eigen::MatrixXd meshOneVertexNormals_prime = transformMat_k * meshOneVertexNormals;
	// find intersection qi_k on surface Q, based on normal lines !

  // *********************************************************** 
  // Plot the initial mesh 
  // *********************************************************** 

  igl::viewer::Viewer viewer;
  viewer.data.set_mesh(V_one, F_one);
  viewer.launch();
}
  
  /*
  Eigen::Matrix3f m = Matrix3f::Random();
  std::cout << m << std::endl;

  Eigen::MatrixXf homogenous = convertToHomogenousForm(m);
  std::cout << "\n\n\n" << std::endl;
  std::cout << homogenous << std::endl;

  Eigen::MatrixXf deHomogenous = normalizeHomogenousMatrix(homogenous);
  std::cout << "\n\n\n" << std::endl;
  std::cout << deHomogenous << std::endl;
*/




