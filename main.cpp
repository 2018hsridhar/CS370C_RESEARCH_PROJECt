// what is this :: a greedy surface reconstruction algorithm, from 2 range images
// solve a closest points scheme, for both sets of vertices, and just add edges. 
#include <igl/readOBJ.h>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <igl/readSTL.h>
#include <igl/viewer/Viewer.h>
#include "tutorial_shared_path.h"
#include <igl/writeOFF.h>
#include <igl/writePLY.h>
#include <igl/cat.h> 
#include <igl/exterior_edges.h> 
#include <igl/is_boundary_edge.h> 
#include <igl/on_boundary.h> 
#include <igl/point_mesh_squared_distance.h>
#include <igl/adjacency_list.h>
#include <igl/slice.h>
#include <igl/slice_mask.h>

#include <igl/is_border_vertex.h>
#include <set> 

#include <algorithm>
#include <iterator> 

using namespace Eigen;  
using namespace std;
using namespace igl;

// a set of useful methods for this code
int countBoundaryVertices(vector<bool> boundaryVerticesStatus);
void fillBoundaryIndexes(vector<bool> verticesBoundaryStatus, int numBoundaryVertices,int* indexesArray);

struct Mesh
{
  Eigen::MatrixXd V; 
  Eigen::MatrixXi F;
  Eigen::MatrixXi E; 
} scan1,scan2,scans,scene,interpolatedSurface;

int main(int argc, char *argv[])
{
  //if(!readOFF(TUTORIAL_SHARED_PATH "/camelhead.off",scan1.V,scan1.F))
  if(!readOFF(TUTORIAL_SHARED_PATH "/planexy.off",scan1.V,scan1.F))
  {
    cout<<"failed to load partial scan one "<<endl;
  } 
  //if(!readOFF(TUTORIAL_SHARED_PATH "/camelhead2.off",scan2.V,scan2.F))
  if(!readOFF(TUTORIAL_SHARED_PATH "/planexy2.off",scan2.V,scan2.F))
  {
    cout<<"failed to load partial scan two "<<endl;
  }


  // solve for vertex adjacency lists of the two partial scans ... will be used for determining minimal edges
  vector<vector<double>> Adjacency_Scan1;
  igl::adjacency_list(scan1.F,Adjacency_Scan1);
   
  vector<vector<double>> Adjacency_Scan2;
  igl::adjacency_list(scan2.F,Adjacency_Scan2);

  // discover boundary vertices
  std::vector<bool> boundaryVerticesStatus_scan1 = igl::is_border_vertex(scan1.V, scan1.F);  
  std::vector<bool> boundaryVerticesStatus_scan2 = igl::is_border_vertex(scan2.V, scan2.F);  

  // count # of boundary vertices 
  int numBoundaryVerticesScan1 = countBoundaryVertices(boundaryVerticesStatus_scan1);
  int numBoundaryVerticesScan2 = countBoundaryVertices(boundaryVerticesStatus_scan2);
  int totalNumBoundaryVertices = numBoundaryVerticesScan1 + numBoundaryVerticesScan2;

  // sovle for boundary vertex indices
  int boundaryVerticesIdxs_scan1_array[numBoundaryVerticesScan1];  
  fillBoundaryIndexes(boundaryVerticesStatus_scan1,numBoundaryVerticesScan1,boundaryVerticesIdxs_scan1_array);

  int boundaryVerticesIdxs_scan2_array[numBoundaryVerticesScan2];  
  fillBoundaryIndexes(boundaryVerticesStatus_scan2,numBoundaryVerticesScan2,boundaryVerticesIdxs_scan2_array);

// since we know total # boundary vertices ... we know how many vertice and faces the interpolating surface will have! 
  igl::cat(1,scan1.V,scan2.V,interpolatedSurface.V);
   interpolatedSurface.F = Eigen::MatrixXi::Zero(totalNumBoundaryVertices,3);  // #TODO :: check if this must be changed !

  // construct the sets of boundary vertices
  // we can pre-alloc a ZERO() matrix of a given size, as we know the # of boundary vertices
  // then iteratively replace ith row with corresponding boundary vertex !

// I will leave this section as is, BUT it should be refactorzed!
  Eigen::MatrixXd boundaryVertices_scan1 = Eigen::MatrixXd::Zero(numBoundaryVerticesScan1,3);  
  int idx = 0;
  int i;
  for ( i = 0; i < numBoundaryVerticesScan1; i++)
  {
     int indexIntoVerticesOfScan1 = boundaryVerticesIdxs_scan1_array[i];
     if ( indexIntoVerticesOfScan1 != -1 ) 
      {
		boundaryVertices_scan1.row(idx) = (scan1.V).row(indexIntoVerticesOfScan1); 
        idx++;
      }
  }

  Eigen::MatrixXd boundaryVertices_scan2 = Eigen::MatrixXd::Zero(numBoundaryVerticesScan2,3); 
  int idx2 = 0;
  for ( i = 0; i < numBoundaryVerticesScan2; i++)
  {
     int indexIntoVerticesOfScan2 = boundaryVerticesIdxs_scan2_array[i];
     if ( indexIntoVerticesOfScan2 != -1 ) 
      {
		boundaryVertices_scan2.row(idx2) = (scan2.V).row(indexIntoVerticesOfScan2); 
        idx2++;
      }
  }

  //std::cout << " boundary vertices, scan 1 are " << std::endl;
  //std::cout << boundaryVertices_scan1 << std::endl;

/////////////////////////////////////////////////////////////////////////////
// CAN CONFIRM :: I seem to be getting the correct set of boundary vertices  !
/////////////////////////////////////////////////////////////////////////////

/*
  [1] solve for a seed edge :: choose a rand point in scan_1, find closest point in scan_2
  [2] keep alternating edge solving , and use the adjacency lists  
  [3] end once you have the original edge data !
 */

// [1] solve for a seed edge :: choose a rand point in scan_1, find closest point in scan_2
  int scan1SeedPointIndex = boundaryVerticesIdxs_scan1_array[0];
  Eigen::MatrixXd scan1SeedPoint = boundaryVertices_scan1.row(0);  
  Eigen::MatrixXd closestPointsFrom_Scan1_To_Scan2;
  Eigen::VectorXi Ele = Eigen::VectorXi::LinSpaced(boundaryVertices_scan2.rows(), 0, boundaryVertices_scan2.rows() - 1);
  Eigen::VectorXd smallestSquaredDists;
  Eigen::VectorXi smallestDistIndxs;
  igl::point_mesh_squared_distance(scan1SeedPoint,boundaryVertices_scan2,
                                    Ele,
									smallestSquaredDists,smallestDistIndxs,
									closestPointsFrom_Scan1_To_Scan2);

  int scan2ClosestPointToSeedIndex = smallestDistIndxs(0,0) + scan1.V.rows(); 
  Eigen::VectorXi seedEdge = Eigen::Vector2i( scan1SeedPointIndex, scan2ClosestPointToSeedIndex);

  //std::cout << seedEdge << std::endl;





//////////////////////////////////////////////////////////////////////////////
// OLD ALGORITHM :: THIS WAS NOT CORRECT IN THE FIRST PLACE //////////////////
//////////////////////////////////////////////////////////////////////////////

  /* For scan 1 :: 
   * [a] SOLVE for closest vertex in scan 2
   * [b] SOVLE for closest vertex in boundaryVertices_scan1 ( note :: has to be 1 put, else, feeding back to self again issue !)
   * CONSTRUCT a triangle, corresponding to the three indices found here !
   */
 
  // part (a)  
/*
  Eigen::MatrixXd closestPointsFrom_Scan1_To_Scan2;
  Eigen::VectorXi Ele = Eigen::VectorXi::LinSpaced(boundaryVertices_scan2.rows(), 0, boundaryVertices_scan2.rows() - 1);
  Eigen::VectorXd smallestSquaredDists;
  Eigen::VectorXi smallestDistIndxs;
  igl::point_mesh_squared_distance(boundaryVertices_scan1,boundaryVertices_scan2,
                                    Ele,
									smallestSquaredDists,smallestDistIndxs,
									closestPointsFrom_Scan1_To_Scan2);

  // part (b)  
  Eigen::MatrixXd closestBoundaryPointIn_Scan1;
  Eigen::VectorXi Ele_Scan1 = Eigen::VectorXi::LinSpaced(boundaryVertices_scan1.rows(), 0, boundaryVertices_scan1.rows() - 1);
  Eigen::VectorXd smallestSquaredDists_Scan1;
  Eigen::VectorXi smallestDistIndxs_Scan1;

  int j = 0;
  for ( int i = 0; i < numBoundaryVerticesScan1; i++)
  {
    int scan1_boundaryPoint = boundaryVerticesIdxs_scan1_array[i];

    Eigen::MatrixXd scan1CurrentPoint = boundaryVertices_scan1.row(i);  
    // since always closest to self; use dummyVertices, to set the current BOUNDARY VERTEX to very huge (x,y,z) vals! 
    Eigen::MatrixXd dummyVertices = boundaryVertices_scan1;
    dummyVertices.row(i) = Eigen::Vector3d(10000,100000,100000); 
  	igl::point_mesh_squared_distance(scan1CurrentPoint,dummyVertices,
                                    	Ele_Scan1,
										smallestSquaredDists_Scan1,smallestDistIndxs_Scan1,
										closestBoundaryPointIn_Scan1);

    // construct face data ... output to file
    int scan2_closestPoint = smallestDistIndxs(i,0) + scan1.V.rows(); 
    int scan1_closestPoint = boundaryVerticesIdxs_scan1_array[smallestDistIndxs_Scan1(0,0)]; // needs to be fixed!
    Eigen::VectorXi newFace = Eigen::Vector3i( scan1_boundaryPoint, scan1_closestPoint, scan2_closestPoint);
	(interpolatedSurface.F).row(j) = newFace.transpose();
    j++; 
  }
*/

	/////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////

  // @ this j, we are now working with the 2nd partial scan
  /* For scan 2 :: 
   * [a] SOLVE for closest vertex in scan 1 
   * [b] SOVLE for closest vertex in boundaryVertices_scan2 ( note :: has to be 1 put, else, feeding back to self again issue !)
   * CONSTRUCT a triangle, corresponding to the three indices found here !
   */

/* 
  // part (a)  
  Eigen::MatrixXd closestPointsFrom_Scan2_To_Scan1;
  Eigen::VectorXi Ele_Prime = Eigen::VectorXi::LinSpaced(boundaryVertices_scan1.rows(), 0, boundaryVertices_scan1.rows() - 1);
  Eigen::VectorXd smallestSquaredDists_Prime;
  Eigen::VectorXi smallestDistIndxs_Prime;
  igl::point_mesh_squared_distance(boundaryVertices_scan2,boundaryVertices_scan1,
                                    Ele_Prime,
									smallestSquaredDists_Prime,smallestDistIndxs_Prime,
									closestPointsFrom_Scan2_To_Scan1);

  // part (b)  
  Eigen::MatrixXd closestBoundaryPointIn_Scan2;
  Eigen::VectorXi Ele_Scan2 = Eigen::VectorXi::LinSpaced(boundaryVertices_scan2.rows(), 0, boundaryVertices_scan2.rows() - 1);
  Eigen::VectorXd smallestSquaredDists_Scan2;
  Eigen::VectorXi smallestDistIndxs_Scan2;

  int k = j;
  for ( int i = 0; i < numBoundaryVerticesScan2; i++)
  {
    int scan2_boundaryPoint = boundaryVerticesIdxs_scan2_array[i] + scan1.V.rows();

    Eigen::MatrixXd scan2CurrentPoint = boundaryVertices_scan2.row(i);  
    // since always closest to self; use dummyVertices, to set the current BOUNDARY VERTEX to very huge (x,y,z) vals! 
    Eigen::MatrixXd dummyVertices_prime = boundaryVertices_scan2;
    dummyVertices_prime.row(i) = Eigen::Vector3d(10000,100000,100000); 
  	igl::point_mesh_squared_distance(scan2CurrentPoint,dummyVertices_prime,
                                    	Ele_Scan2,
										smallestSquaredDists_Scan2,smallestDistIndxs_Scan2,
										closestBoundaryPointIn_Scan2);

    // construct face data ... output to file
    int scan1_closestPoint = smallestDistIndxs(i,0); 
    int scan2_closestPoint = boundaryVerticesIdxs_scan2_array[smallestDistIndxs_Scan2(0,0)] + scan1.V.rows(); // needs to be fixed!
    Eigen::VectorXi newFace_prime = Eigen::Vector3i( scan2_boundaryPoint, scan2_closestPoint, scan1_closestPoint);
	(interpolatedSurface.F).row(k) = newFace_prime.transpose();
    k++; 
  }
*/

  // CREATE ONE HUGE MESH containing the two partial scans and interpoalted surface
  igl::cat(1,scan1.V,scan2.V,scans.V);
  igl::cat(1,scans.V,interpolatedSurface.V,scene.V); 

  igl::cat(1,scan1.F, MatrixXi(scan2.F.array() + scan1.V.rows()), scans.F);
  igl::cat(1,scans.F, interpolatedSurface.F, scene.F);

  /***********************************************************/ 
  // SETUP LibIgl Viewer 
  /***********************************************************/ 
  igl::viewer::Viewer viewer;
  viewer.data.set_mesh(scene.V, scene.F); 
  viewer.launch();
 
}

int countBoundaryVertices(vector<bool> boundaryVerticesStatus)
{
  int numBoundaryVertices = 0;
  int size = boundaryVerticesStatus.size();
  for ( int i = 0; i < size; ++i)
      if(boundaryVerticesStatus[i] ) 
          numBoundaryVertices++;
  return numBoundaryVertices;
}

void fillBoundaryIndexes(vector<bool> verticesBoundaryStatus, int numBoundaryVertices,int* indexesArray)
{
  for ( int i = 0; i < numBoundaryVertices; i++)
      indexesArray[i] = -1;

  int i;
  int cur = 0;
  for ( i = 0; i < verticesBoundaryStatus.size(); ++i) {
      if(verticesBoundaryStatus[i] ) {
          indexesArray[cur] = i;
          cur++;
      }
  }
} 







// create Slice Stack ( this could have been approach, but Nathan Clement's code would need to be updated ! ) 
  // generate the enclosing bounding box for the two mesh manifolds/partial scan
  //     - this is automatically created from SliceStack constructor 
/*
  std::string objectName = "toSolveForInterpSurfaceMesh";
  std::string pathToSlices = "/u/hari2018/CS370C/LIBS/libigl/tutorial/CS370C_RESEARCH_PROJECt/build/partial_scans";

  SliceStack ss(pathToSlices.c_str(), objectName.c_str());


  Eigen::MatrixXd TV;
  Eigen::MatrixXd TF;
  Eigen::MatrixXi TT; // not sure ( tetrahedrals, I believe ... not the same as vertices | faces, i think ... unsure ) 
  Eigen::MatrixXi TO; // not sure ( original markers ... what do they represent ) ? 

  // and tetrahedralize the common surface

  SliceStack::tetrahedralizeSlice( 
      scan1.V, scan1.F,
      scan2.V, scan2.F,
      scan1.O, scan2.O,
      TV,TT,TF,TO);

   generate a set of "pseudo"-temperature values for each vertex, via [ heat-flow ] 
  int slice_no = 0; // huh?   what does this do, exactly? 
  Eigen::VectorXd Temperatures;
  SliceStack::computeLaplace(slice_no,TV,TF,TT,TO,Temperatures);
*/

// include files for code , from project Tiling 
/*
#include "tiling/SliceStack.h" // I believe it has too be this format. I cannot just to "SliceStack.h" directly, though 
#include "tiling/SliceStack.cpp" // not correct #TODO fix linking s.t. files can be included properly 
#include "tiling/curvatureFlow.h"
#include "tiling/glob_defs.h"
#include "tiling/offsetSurface.h"
#include "tiling/viewTetMesh.h"
#include "tiling/Helpers.h"
*/

  //std::vector<bool> v1 = igl::is_border_vertex(scan1.V, scan1.F);  
  //Eigen::VectorXi J (v1.data());
  //Eigen::VectorXi J = VectorXi::Map(v1.data(), v1.size());
  //const Eigen::Array<bool,Eigen::Dynamic,1> keep = J.array(); ... do this, but once u have everything else working !

/*
  std::cout << scan1.V << std::endl;
  Eigen::MatrixXd partialImage = igl::slice_mask(scan1.V,keep,1);
  std::cout << partialImage << std::endl;
  return 0;
*/




