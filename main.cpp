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

#include <algorithm>
#include <iterator> // ostream_iterator

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
  Eigen::MatrixXd V; 
  Eigen::MatrixXi F;
  Eigen::MatrixXi E; // (Ex2) dimensinoal matrix of all edges of a particular mesh Eigen::MatrixXd O; // (O) not sure, but it is needed for tetraHedralize Slice method in "tiling" . what data must go in? not sure !  } scan2,scan1,scene;
} scan1,scan2,scans,scene,interpolatedSurface;

// in my case, I expect the # boundary vertices = 16 = # boundary edges ( for the planar data case )
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
  // note :: convert to a set of indices of boundary_vertices :: will be easier to solve this problem then
  std::vector<bool> boundaryVerticesStatus_scan1 = igl::is_border_vertex(scan1.V, scan1.F);  
  std::vector<bool> boundaryVerticesStatus_scan2 = igl::is_border_vertex(scan2.V, scan2.F);  

/*
    MY SANITY CHECK for boolean vector of boundary vertex status  ... this makes sense to me 
	for (auto const& c : boundaryVerticesStatus_scan1)
       std::cout << c << std::endl;

	for (auto const& c : boundaryVerticesStatus_scan2)
       std::cout << c << std::endl;
    return 0;
*/

/*
// TECHNICALLY, this portion is not needed in the code ... its just useful for collecting statistics & sanity checking !
  Eigen::MatrixXi boundaryEdges_Scan1;
  Eigen::MatrixXi boundaryEdges_Scan2;
  igl::exterior_edges(scan1.F,boundaryEdges_Scan1); 
  igl::exterior_edges(scan2.F,boundaryEdges_Scan2); 
  int numBoundaryEdgesScan1 = boundaryEdges_Scan1.rows();
  int numBoundaryEdgesScan2 = boundaryEdges_Scan2.rows(); 

  std::cout << "# boundary edges , scan 1 = " << numBoundaryEdgesScan1 << std::endl;
  std::cout << "# boundary edges , scan 2 = " << numBoundaryEdgesScan2 << std::endl; // so this makes sense too me! 
*/

  // CONVERT <boolean> vector to <int> vector , as the boundary vertex indices 
  // WILL GENERATE the boundary interpolating surface mesh

  int boundaryVerticesIdxs_scan1_array[scan1.V.rows()];  // I can get away with this, but its not good practice !
  for ( int i = 0; i < scan1.V.rows(); i++)
      boundaryVerticesIdxs_scan1_array[i] = -1;

	// this dummy is used to find idxs quickly, ( access in sets is O(n), arrays is O(1))
  set<int> boundaryVerticesIdxs_scan1;
  int i;
  int numBoundaryVerticesScan1 = 0;
  for ( i = 0; i < scan1.V.rows(); ++i)
  {
      if(boundaryVerticesStatus_scan1[i] ) {
          boundaryVerticesIdxs_scan1_array[i] = i;
          boundaryVerticesIdxs_scan1.insert(boundaryVerticesIdxs_scan1.end(),i);
          numBoundaryVerticesScan1++;
      }
  }

  int boundaryVerticesIdxs_scan2_array[scan2.V.rows()];
  for ( int i = 0; i < scan1.V.rows(); i++)
      boundaryVerticesIdxs_scan2_array[i] = -1;

  set<int> boundaryVerticesIdxs_scan2;
  int j;
  int numBoundaryVerticesScan2 = 0;
  for ( j = 0; j < scan2.V.rows(); ++j)
  {
      if(boundaryVerticesStatus_scan2[j] ) {
          boundaryVerticesIdxs_scan2_array[j] = j;
          boundaryVerticesIdxs_scan2.insert(boundaryVerticesIdxs_scan2.end(),j);
          numBoundaryVerticesScan2++;
      }
  } 

 // GLARING BUG :: boundaryverticesStatus_scan1 is a bool vector ... u have to count though !
  int totalNumBoundaryVertices = numBoundaryVerticesScan1 + numBoundaryVerticesScan2;
  std::cout << "# boundary vertices , scan 1 = " << numBoundaryVerticesScan1 << std::endl;
  std::cout << "# boundary vertices , scan 2 = " << numBoundaryVerticesScan2 << std::endl;

  std::cout << "# vertices in boundary VertexIdx, set 1 , scan 1 = " << boundaryVerticesIdxs_scan1.size() << std::endl;
  std::cout << "# vertices in boundary VertexIdx, set 2 , scan 2 = " << boundaryVerticesIdxs_scan2.size() << std::endl;


   std::ostream_iterator< int > output( cout, " " );

   //cout << "set 1, scan 1, contains: ";
   //std::copy( boundaryVerticesIdxs_scan1.begin(), boundaryVerticesIdxs_scan1.end(), output ); // this is right ! missing (8th), which is expected ! 
   //std::copy( boundaryVerticesIdxs_scan2.begin(), boundaryVerticesIdxs_scan2.end(), output ); // this is right ! missing (8th), which is expected ! 
//return 0;


// since we know total # boundary vertices ... we know how many vertice and faces the interpolating surface will have! 
// #TODO :: expand to both vertex sets !
  interpolatedSurface.V = Eigen::MatrixXd::Zero(numBoundaryVerticesScan1,3);
  interpolatedSurface.F = Eigen::MatrixXi::Zero(numBoundaryVerticesScan1,3); 

  // construct the sets of boundary vertices
  // we can pre-alloc a ZERO() matrix of a given size, as we know the # of boundary vertices
  // then iteratively replace ith row with corresponding boundary vertex !

  // 	SO THIS IS WHERE THINGS START TO BREAK !

  Eigen::MatrixXd boundaryVertices_scan1 = Eigen::MatrixXd::Zero(numBoundaryVerticesScan1,3);  
  int idx = 0;
  for ( i = 0; i < scan1.V.rows(); i++)
  {
     // #TODO rain check this code
     if ( boundaryVerticesIdxs_scan1_array[i] == i )   // this is a bug ( 0th vertex unavaille ... gaah ! ) 
      {
		boundaryVertices_scan1.row(idx) = (scan1.V).row(i); // this is where THE ISSUE is @ 
        idx++;
        /*
        std::cout << "i = " << i << std::endl;
        std::cout << "Scan data " << (scan1.V).row(i) << std::endl;
        std::cout << "boundary  " << boundaryVertices_scan1.row(idx) << std::endl;
        */
      }
  }
/*
  std::cout << " boundary vertices, scan 1 are " << std::endl;
	//for (auto const& c : boundaryVerticesStatus_scan1)
     //  std::cout << c << std::endl;
  std::cout << boundaryVertices_scan1 << std::endl;
*/

  //std::cout << "the 8thy vertex = " << (scan1.V).row(8) << std::endl; = (0,0,0) 
  Eigen::MatrixXd boundaryVertices_scan2 = Eigen::MatrixXd::Zero(numBoundaryVerticesScan2,3); 
  int idx2 = 0;
  for ( i = 0; i < scan2.V.rows(); i++)
  {
      if ( boundaryVerticesIdxs_scan2_array[i] == i ) 
      {
		boundaryVertices_scan2.row(idx2) = (scan2.V).row(i);
        idx2++;
      }
  }

// CAN CONFIRM :: I seem to be getting the correct set of boundary vertices  !

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
    Eigen::MatrixXd scan1CurrentPoint = boundaryVertices_scan1.row(i);  
  	igl::point_mesh_squared_distance(scan1CurrentPoint,boundaryVertices_scan1,
                                    	Ele_Scan1,
										smallestSquaredDists_Scan1,smallestDistIndxs_Scan1,
										closestBoundaryPointIn_Scan1);

    // construct face data ... output to file
    int scan1_boundaryPoint = boundaryVerticesIdxs_scan1_array[i];
    int scan1_closestPoint = smallestDistIndxs_Scan1(0,0);
    int scan2_closestPoint = smallestDistIndxs(i,0);  // assert if correct
    Eigen::VectorXi newFace = Eigen::Vector3i( scan1_boundaryPoint, scan1_closestPoint, scan2_closestPoint);
	(interpolatedSurface.F).row(j) = newFace.transpose();
    j++; 
  }
  // @ this j, we are now working with the 2nd partial scan

  // CREATE ONE HUGE MESH containing both partial scan pieces ( inspired by example 407 ) 
  // and the interpolated surface

  igl::cat(1,scan1.V,scan2.V,scene.V);
  //igl::cat(1,scene.V,interpolatedSurface.V,scene.V); // technically, we are not adding new vertices to our system !, so why is this even here??

  igl::cat(1,scan1.F, MatrixXi(scan2.F.array() + scan1.V.rows()), scans.F);
  igl::cat(1,scans.F, MatrixXi(scan2.F.array() + scene.V.rows()), scene.F);

 // issue arises with interpolatedSurface.F.array(), but not scan2.F.array() ... I wonder why though! 
  //igl::cat(1,scans.F, MatrixXi(scan2.F.array() + scene.V.rows()), scene.F);  
  //igl::cat(1,scans.F, MatrixXi(scan2.F.array() + scan1.V.rows()), scene.F);  // GETTING ISSUE HERE! ( works with scans.V.rows(), not scene.V.rows() ??? WTF !! ) 

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
      scan1.V, scan1.F,
      scan2.V, scan2.F,
      scan1.O, scan2.O,
      TV,TT,TF,TO);

   generate a set of "pseudo"-temperature values for each vertex, via [ heat-flow ] 
  int slice_no = 0; // huh?   what does this do, exactly? 
  Eigen::VectorXd Temperatures;
  SliceStack::computeLaplace(slice_no,TV,TF,TT,TO,Temperatures);
*/
