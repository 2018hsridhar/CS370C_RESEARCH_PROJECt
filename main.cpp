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
#include <igl/cat.h> // used to help concatenate sets of vertices, edges, normals ( matrices in general ... MATLAB-func based )
#include <igl/exterior_edges.h> // used to solve for set of all edges  ( this nmight get boundary edges  ?? )
#include <igl/is_boundary_edge.h> 
#include <igl/on_boundary.h> 
#include <igl/point_mesh_squared_distance.h>
#include <igl/adjacency_list.h>

#include <igl/is_border_vertex.h>
#include <set> 

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

using namespace Eigen;  
using namespace std;
using namespace igl;

struct Mesh
{
  Eigen::MatrixXd V,U,N,UV; // really, N, UV are very useless @ ... to be deleted? not sure exactly atm !
  Eigen::MatrixXi T,F;
  Eigen::MatrixXi E; // (Ex2) dimensinoal matrix of all edges of a particular mesh Eigen::MatrixXd O; // (O) not sure, but it is needed for tetraHedralize Slice method in "tiling" . what data must go in? not sure !  } ahead,behind,scene;
} behind,ahead,scans,scene,interpolatedSurface;

// in my case, I expect the # boundary vertices = 16 ( turns out 25 ... I wonder why ) 
// and I expect the # boundary edges = 16 ! 
int main(int argc, char *argv[])
{
  //if(!readOFF(TUTORIAL_SHARED_PATH "/camelhead.off",behind.V,behind.F))
  if(!readOFF(TUTORIAL_SHARED_PATH "/planexy.off",behind.V,behind.F))
  {
    cout<<"failed to load parital scan one "<<endl;
  } 
  //if(!readOFF(TUTORIAL_SHARED_PATH "/camelhead2.off",ahead.V,ahead.F))
  if(!readOFF(TUTORIAL_SHARED_PATH "/planexy2.off",ahead.V,ahead.F))
  {
    cout<<"failed to load partial scan two "<<endl;
  }

  // solve for vertex adjacency lists of the two partial scans ... will be used for determining minimal edges
  vector<vector<double>> Adjacency_Scan1;
  igl::adjacency_list(behind.F,Adjacency_Scan1);
   
  vector<vector<double>> Adjacency_Scan2;
  igl::adjacency_list(ahead.F,Adjacency_Scan2);

  // discover boundary vertices
  // note :: convert to a set of indices of boundary_vertices :: will be easier to solve this problem then
  std::vector<bool> boundaryVerticesStatus_scan1 = igl::is_border_vertex(behind.V, behind.F);  
  std::vector<bool> boundaryVerticesStatus_scan2 = igl::is_border_vertex(ahead.V, ahead.F);  

  int numBoundaryVerticesScan1 = boundaryVerticesStatus_scan1.size();
  int numBoundaryVerticesScan2 = boundaryVerticesStatus_scan2.size();
  int totalNumBoundaryVertices = numBoundaryVerticesScan1 + numBoundaryVerticesScan2;

  std::cout << "# boundary vertices , scan 1 = " << numBoundaryVerticesScan1 << std::endl;
  std::cout << "# boundary vertices , scan 2 = " << numBoundaryVerticesScan2 << std::endl;

// this portion is not needed in the code ... its just useful for collecting statistics 
  Eigen::MatrixXi boundaryEdges_Scan1;
  Eigen::MatrixXi boundaryEdges_Scan2;
  igl::exterior_edges(behind.F,boundaryEdges_Scan1); 
  igl::exterior_edges(ahead.F,boundaryEdges_Scan2); 
  int numBoundaryEdgesScan1 = boundaryEdges_Scan1.rows();
  int numBoundaryEdgesScan2 = boundaryEdges_Scan2.rows(); 

  std::cout << "# boundary edges , scan 1 = " << numBoundaryEdgesScan1 << std::endl;
  std::cout << "# boundary edges , scan 2 = " << numBoundaryEdgesScan2 << std::endl; // so this makes sense too me! 

  // since we know total # boundary vertices ... we know how many vertice and faces the interpolating surface will have! 
  //interpolatedSurface.V = Eigen::MatrixXd::Identity(totalNumBoundaryVertices,3);
  //interpolatedSurface.F = Eigen::MatrixXi::Identity(totalNumBoundaryVertices,3); // is this the correct # of cols??
  interpolatedSurface.V = Eigen::MatrixXd::Zero(numBoundaryVerticesScan1,3);
  interpolatedSurface.F = Eigen::MatrixXi::Zero(numBoundaryVerticesScan1,3); 


  // convert <boolean> vector to <int> vector , as the indices of the boundary vertices will be needed to generate the 
  // boundary interpolating surface mesh
  // SCREW SETS for requiring iterator,s instead of being random-access , like ARRAYS

  int boundaryVerticesIdxs_scan1_array[numBoundaryVerticesScan1]; 
	// this dummy is used to find idxs quickly, ( access in sets is O(n), arrays is O(1))
  set<int> boundaryVerticesIdxs_scan1;
  int i;
  for ( i = 0; i < numBoundaryVerticesScan1; ++i)
  {
      if(boundaryVerticesStatus_scan1[i] ) {
          boundaryVerticesIdxs_scan1_array[i] = i;
          boundaryVerticesIdxs_scan1.insert(boundaryVerticesIdxs_scan1.end(),i);
      }
  }

  int boundaryVerticesIdxs_scan2_array[numBoundaryVerticesScan2];
  set<int> boundaryVerticesIdxs_scan2;
  int j;
  for ( j = 0; j < numBoundaryVerticesScan2; ++j)
  {
      if(boundaryVerticesStatus_scan2[j] ) {
          boundaryVerticesIdxs_scan2.insert(boundaryVerticesIdxs_scan2.end(),i);
          boundaryVerticesIdxs_scan2_array[i] = i;
      }
  } 



  // construct the sets of boundary vertices
  // we can pre-alloc a ZERO() matrix of a given size, as we know the # of boundary vertices
  // then iteratively replace ith row with corresponding boundary vertex !

  // 	SO THIS IS WHERE THINGS START TO BREAK !

  Eigen::MatrixXd boundaryVertices_scan1 = Eigen::MatrixXd::Zero(numBoundaryVerticesScan1,3); 
  int idx = 0;
  for ( i = 0; i < numBoundaryVerticesScan1; ++i)
  {
     // #TODO rain check thsi code
     if ( boundaryVerticesIdxs_scan1_array[i] == i )  
      {
		boundaryVertices_scan1.row(idx) = (behind.V).row(i); // this is where THE ISSUE is @ 
      }
      idx++;
  }
  
  Eigen::MatrixXd boundaryVertices_scan2 = Eigen::MatrixXd::Zero(numBoundaryVerticesScan2,3); 
  int idx2 = 0;
  for ( i = 0; i < numBoundaryVerticesScan2; ++i)
  {
      if ( boundaryVerticesIdxs_scan2_array[i] == i ) 
      {
		boundaryVertices_scan2.row(idx2) = (ahead.V).row(i);
      }
      idx2++;
  }

  /* For scan 1 :: 
   * [a] SOLVE for closest vertex in scan 2
   * [b] SOVLE for closest vertex in boundaryVertices_scan1 ( note :: has to be 1 put, else, feeding back to self again issue !)
   * CONSTRUCT a triangle, corresponding to the three indices found here !
   */
 
  // part (a)  
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


//////////// HUGE NOTE TO SELF /////////////////////
///////////// WHEN OUTPUTTING THIS FILE //////////////
///////////// Scan1's vertice first, then Scan's 2 //////////////
//////////// ACCOUNT FOR THIS OFFSET /////////////// 
// total number of face to add = total number of boundary vertices !  // 

  j = 0;
  for ( int i = 0; i < numBoundaryVerticesScan1; i++)
  {
/* PART 1 
    Eigen::MatrixXd scan1CurrentPoint = boundaryVertices_scan1.row(i);  
  	igl::point_mesh_squared_distance(scan1CurrentPoint,boundaryVertices_scan1,
                                    	Ele_Scan1,
										smallestSquaredDists_Scan1,smallestDistIndxs_Scan1,
										closestBoundaryPointIn_Scan1);
*/

/* PART 2
    // construct face data ... output to file
    int scan1_boundaryPoint = boundaryVerticesIdxs_scan1_array[i];
    int scan1_closestPoint = smallestDistIndxs_Scan1(0,0);
    int scan2_closestPoint = smallestDistIndxs(i,0);  // assert if correct

    Eigen::VectorXi newFace = Eigen::Vector3i( scan1_boundaryPoint, scan1_closestPoint, scan2_closestPoint);
	(interpolatedSurface.F).row(j) = newFace.transpose();
*/
    j++; 
  }
  // @ this j, we are now working with the 2nd partial scan

  // CREATE ONE HUGE MESH containing both partial scan pieces ( inspired by example 407 ) 
  // and the interpolated surface

  igl::cat(1,behind.V,ahead.V,scene.V);
  //igl::cat(1,scene.V,interpolatedSurface.V,scene.V); // technically, we are not adding new vertices to our system !, so why is this even here??

  igl::cat(1,behind.F, MatrixXi(ahead.F.array() + behind.V.rows()), scans.F);
  igl::cat(1,scans.F, MatrixXi(interpolatedSurface.F.array() + scene.V.rows()), scene.F);  // GETTING ISSUE HERE! ( works with scans.V.rows(), not scene.V.rows() ??? WTF !! ) 

  /***********************************************************/ 
  // SETUP LibIgl Viewer 
  /***********************************************************/ 
  igl::viewer::Viewer viewer;
  viewer.data.set_mesh(scene.V, scene.F); 
  viewer.launch();
 

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
      behind.V, behind.F,
      ahead.V, ahead.F,
      behind.O, ahead.O,
      TV,TT,TF,TO);

   generate a set of "pseudo"-temperature values for each vertex, via [ heat-flow ] 
  int slice_no = 0; // huh?   what does this do, exactly? 
  Eigen::VectorXd Temperatures;
  SliceStack::computeLaplace(slice_no,TV,TF,TT,TO,Temperatures);
*/
