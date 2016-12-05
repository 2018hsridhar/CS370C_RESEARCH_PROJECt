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

int countBoundaryVertices(vector<bool> boundaryVerticesStatus);

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
  // count # of boundary vertices ( put this in its own method ! ) 
  int numBoundaryVerticesScan1 = countBoundaryVertices(boundaryVerticesStatus_scan1);
  int numBoundaryVerticesScan2 = countBoundaryVertices(boundaryVerticesStatus_scan2);

  //std::cout << "num Boundary vertices , scan 1 = " << numBoundaryVerticesScan1 << std::endl;
  //std::cout << "num Boundary vertices , scan 2 = " << numBoundaryVerticesScan2 << std::endl;

  /////////////////////////////////////////////
  /////////////////////////////////////////////

  // CONVERT <boolean> vector to <int> vector , as the boundary vertex indices 
  // WILL GENERATE the boundary interpolating surface mesh
  // array , for partial scan 1 , representing idx data ( into scan1 vertices ) for boundary Vertices
  int boundaryVerticesIdxs_scan1_array[numBoundaryVerticesScan1];  
  for ( int i = 0; i < numBoundaryVerticesScan1; i++)
      boundaryVerticesIdxs_scan1_array[i] = -1;

	// this dummy is used to find idxs quickly, ( access in sets is O(n), arrays is O(1))
  set<int> boundaryVerticesIdxs_scan1;
  int i;
  int cur = 0;
  for ( i = 0; i < scan1.V.rows(); ++i)
  {
      if(boundaryVerticesStatus_scan1[i] ) 
	  {
          boundaryVerticesIdxs_scan1_array[cur] = i;
          boundaryVerticesIdxs_scan1.insert(boundaryVerticesIdxs_scan1.end(),i);
          cur++;
      }
  }
  // array , for partial scan 2 , representing idx data ( into scan2 vertices ) for boundary Vertices
  int boundaryVerticesIdxs_scan2_array[numBoundaryVerticesScan2];  
  for ( int i = 0; i < numBoundaryVerticesScan2; i++)
      boundaryVerticesIdxs_scan2_array[i] = -1; 

  std::cout << std::endl;
  for (int i = 0; i < numBoundaryVerticesScan1; i++)
     std::cout << boundaryVerticesIdxs_scan1_array[i] << std::endl;

  set<int> boundaryVerticesIdxs_scan2;
  int j;
  cur = 0;
  for ( j = 0; j < scan2.V.rows(); ++j)
  {
      if(boundaryVerticesStatus_scan2[j] ) {
          boundaryVerticesIdxs_scan2_array[cur] = j;
          boundaryVerticesIdxs_scan2.insert(boundaryVerticesIdxs_scan2.end(),j);
          cur++;
      }
  } 

 // GLARING BUG :: boundaryverticesStatus_scan1 is a bool vector ... u have to count though !
  int totalNumBoundaryVertices = numBoundaryVerticesScan1 + numBoundaryVerticesScan2;
/*
  std::cout << "# boundary vertices , scan 1 = " << numBoundaryVerticesScan1 << std::endl;
  std::cout << "# boundary vertices , scan 2 = " << numBoundaryVerticesScan2 << std::endl;

  std::cout << "# vertices in boundary VertexIdx, set 1 , scan 1 = " << boundaryVerticesIdxs_scan1.size() << std::endl;
  std::cout << "# vertices in boundary VertexIdx, set 2 , scan 2 = " << boundaryVerticesIdxs_scan2.size() << std::endl;
*/


   //std::ostream_iterator< int > output( cout, " " );
   //cout << "set 1, scan 1, contains: ";
   //std::copy( boundaryVerticesIdxs_scan1.begin(), boundaryVerticesIdxs_scan1.end(), output ); // this is right ! missing (8th), which is expected ! 
   //std::copy( boundaryVerticesIdxs_scan2.begin(), boundaryVerticesIdxs_scan2.end(), output ); // this is right ! missing (8th), which is expected ! 
//return 0;


// since we know total # boundary vertices ... we know how many vertice and faces the interpolating surface will have! 
// #TODO :: expand to both vertex sets !
  igl::cat(1,scan1.V,scan2.V,interpolatedSurface.V);
  interpolatedSurface.F = Eigen::MatrixXi::Zero(numBoundaryVerticesScan1,3); 

  // construct the sets of boundary vertices
  // we can pre-alloc a ZERO() matrix of a given size, as we know the # of boundary vertices
  // then iteratively replace ith row with corresponding boundary vertex !

  // 	SO THIS IS WHERE THINGS START TO BREAK !

  Eigen::MatrixXd boundaryVertices_scan1 = Eigen::MatrixXd::Zero(numBoundaryVerticesScan1,3);  
  int idx = 0;


  for ( i = 0; i < numBoundaryVerticesScan1; i++)
  {
     int indexIntoVerticesOfScan1 = boundaryVerticesIdxs_scan1_array[i];
     if ( indexIntoVerticesOfScan1 != -1 ) 
      {
		boundaryVertices_scan1.row(idx) = (scan1.V).row(indexIntoVerticesOfScan1); // this is where THE ISSUE is @ 
        idx++;
      }
  }

  //std::cout << " boundary vertices, scan 1 are " << std::endl;
  //std::cout << boundaryVertices_scan1 << std::endl;

  //std::cout << "the 8thy vertex = " << (scan1.V).row(8) << std::endl; = (0,0,0) 
  Eigen::MatrixXd boundaryVertices_scan2 = Eigen::MatrixXd::Zero(numBoundaryVerticesScan2,3); 
  int idx2 = 0;
  for ( i = 0; i < numBoundaryVerticesScan2; i++)
  {
     int indexIntoVerticesOfScan2 = boundaryVerticesIdxs_scan2_array[i];
     if ( indexIntoVerticesOfScan2 != -1 ) 
      {
		boundaryVertices_scan2.row(idx2) = (scan2.V).row(indexIntoVerticesOfScan2); // this is where THE ISSUE is @ 
        idx2++;
      }
  }

/////////////////////////////////////////////////////////////////////////////
// CAN CONFIRM :: I seem to be getting the correct set of boundary vertices  !
/////////////////////////////////////////////////////////////////////////////

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
    set<int>::iterator iter = boundaryVerticesIdxs_scan1.find(i);
    int setint;
	if (iter != boundaryVerticesIdxs_scan1.end()) {
		setint = *iter;
	}
    int scan1_boundaryPoint = setint;

    Eigen::MatrixXd scan1CurrentPoint = boundaryVertices_scan1.row(i);  
    // note :: you are always closest to yoursefl ... this won't make sense ! take teh 2nd one  
    Eigen::MatrixXd dummyVertices = boundaryVertices_scan1;
    dummyVertices.row(scan1_boundaryPoint) = Eigen::Vector3d(10000,100000,100000); // why does everybody think vertex 15 is the closest? WTF ???
  	//igl::point_mesh_squared_distance(scan1CurrentPoint,dummyVertices,
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

  std::cout << interpolatedSurface.F << std::endl; // this is clearly wrong !

  // @ this j, we are now working with the 2nd partial scan

  // CREATE ONE HUGE MESH containing both partial scan pieces ( inspired by example 407 ) 
  // and the interpolated surface

  igl::cat(1,scan1.V,scan2.V,scans.V);
  igl::cat(1,scans.V,interpolatedSurface.V,scene.V); 
		// technically, we are not adding new vertices to our system !, so why is this even here??
        // also :: thsi really shouldn't do anything, when i think about it !

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




