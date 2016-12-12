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

#include <Eigen/StdVector>

using namespace Eigen;  
using namespace std;
using namespace igl;

///////// FORWARD METHOD DECLARATIONS /////////////////////
int countBoundaryVertices(vector<bool> boundaryVerticesStatus);
void fillBoundaryIndexes(vector<bool> verticesBoundaryStatus, int numBoundaryVertices,int* indexesArray);
template <typename T>
void printcoll (T const& coll);
Eigen::MatrixXd retrieveBoundaryVerticesInScan1(int bndryVerticesIdxs[], int numBndryVertices);
Eigen::MatrixXd retrieveBoundaryVerticesInScan2(int bndryVerticesIdxs[], int numBndryVertices);

std::vector<int> findAdjBndryVertsInScan1(int vertIdx);
std::vector<int> findAdjBndryVertsInScan2(int vertIdx);

int findClosestAdjBndryNodeInScan1(std::vector<int> bndryVertices, int myVertexIndex);
int findClosestAdjBndryNodeInScan2(std::vector<int> bndryVertices, int myVertexIndex);

///////// PREALLOCATION ///////////////////
struct Mesh
{
  Eigen::MatrixXd V; 
  Eigen::MatrixXi F;
  Eigen::MatrixXi E; 
} scan1,scan2,scans,scene,interpolatedSurface;

std::vector<std::vector<int>> Adjacency_Scan1;
std::vector<std::vector<int>> Adjacency_Scan2;

std::vector<bool> boundaryVerticesStatus_scan1;
std::vector<bool> boundaryVerticesStatus_scan2;

Eigen::MatrixXd boundaryVertices_scan1;
Eigen::MatrixXd boundaryVertices_scan2;

int main(int argc, char *argv[])
{
  //if(!readOFF(TUTORIAL_SHARED_PATH "/camelhead.off",scan1.V,scan1.F))
  if(!readOFF(TUTORIAL_SHARED_PATH "/planexy.off",scan1.V,scan1.F)) {
    cout<<"Failed to load partial scan one."<<endl;
  } 
  //if(!readOFF(TUTORIAL_SHARED_PATH "/camelhead2.off",scan2.V,scan2.F))
  if(!readOFF(TUTORIAL_SHARED_PATH "/planexy2.off",scan2.V,scan2.F)) {
    cout<<"Failed to load partial scan two."<<endl;
  }

  // SOLVE for vertex adjacency lists of the two partial scans 
  // USED to determine minimal edges in surface reconstruction algo
  igl::adjacency_list(scan1.F,Adjacency_Scan1);
  igl::adjacency_list(scan2.F,Adjacency_Scan2);

  // discover boundary vertices
  boundaryVerticesStatus_scan1 = igl::is_border_vertex(scan1.V, scan1.F);  
  boundaryVerticesStatus_scan2 = igl::is_border_vertex(scan2.V, scan2.F);  

  // count # of boundary vertices 
  int numBoundaryVerticesScan1 = countBoundaryVertices(boundaryVerticesStatus_scan1);
  int numBoundaryVerticesScan2 = countBoundaryVertices(boundaryVerticesStatus_scan2);
  int totalNumBoundaryVertices = numBoundaryVerticesScan1 + numBoundaryVerticesScan2;

  // sovle for boundary vertex indices
  int boundaryVerticesIdxs_scan1_array[numBoundaryVerticesScan1];  
  fillBoundaryIndexes(boundaryVerticesStatus_scan1,numBoundaryVerticesScan1,boundaryVerticesIdxs_scan1_array);

  int boundaryVerticesIdxs_scan2_array[numBoundaryVerticesScan2];  
  fillBoundaryIndexes(boundaryVerticesStatus_scan2,numBoundaryVerticesScan2,boundaryVerticesIdxs_scan2_array);

  // CONSTRUCT the sets of boundary vertices
  boundaryVertices_scan1 = retrieveBoundaryVerticesInScan1(boundaryVerticesIdxs_scan1_array,numBoundaryVerticesScan1);
  boundaryVertices_scan2 = retrieveBoundaryVerticesInScan2(boundaryVerticesIdxs_scan2_array,numBoundaryVerticesScan2);
 // #TODO :: assert correctness of this!

  //std::cout << " boundary vertices, scan 1 are " << std::endl;
  //std::cout << boundaryVertices_scan1 << std::endl;

/////////////////////////////////////////////////////////////////////////////
// CAN CONFIRM :: I seem to be getting the correct set of boundary vertices  !
// #TODO :: test this 
/////////////////////////////////////////////////////////////////////////////

  std::vector<Eigen::Vector2i> allEdges;
  std::vector<Eigen::Vector3i> newFaces;

/*
  [1] solve for a seed edge :: choose a rand point in scan_1, find closest point in scan_2
  [2] keep alternating edge solving , and use the adjacency lists  
  [3] end once you have the original edge data !
 */

// [1] solve for a seed edge :: choose a rand point in scan_1, find closest point in scan_2

// #TODO :: include a method for getting ( vertex,index ) easily?? seems useful, but l8r 

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

  int scan2ClosestPointToSeedIndex = smallestDistIndxs(0,0);
  Eigen::Vector2i seedEdge = Eigen::Vector2i( scan1SeedPointIndex, scan2ClosestPointToSeedIndex + scan1.V.rows());
  //std::cout << seedEdge << std::endl;
 
  allEdges.push_back(seedEdge);

  // need a method to tell if there is an existing edge ( in the edges vectors ) 

  // [2] Find shortest adjacent edge, to edge e = (v_1,v_2), and use that to construct the new edge

  // actually, just loop, and check vertIsOnBndry ( see ~/include/igl/loop.cpp)
  // #TODO :: make a method for this!

  std::cout << " Going to find the set of adj bndry verts to (v_1,v_2) " << std::endl;
  std::vector<int> bndryVertsNodeOne = findAdjBndryVertsInScan1(scan1SeedPointIndex);
  std::vector<int> bndryVertsNodeTwo = findAdjBndryVertsInScan2(scan2ClosestPointToSeedIndex);

  // find vertices that are closest to {v_1,v_2} in edge e = (v_1,v_2)
  int indexClosestAdjBndryNodeToNodeOne = findClosestAdjBndryNodeInScan1(bndryVertsNodeOne,scan1SeedPointIndex);
  int indexClosestAdjBndryNodeToNodeTwo = findClosestAdjBndryNodeInScan2(bndryVertsNodeTwo,scan2ClosestPointToSeedIndex);
  //  std::cout << "vertiex adj to boundary vertex [" << scan1SeedPointIndex << "] are : "; 
  //  printcoll(bndryVertsNodeOne);
   // std::cout << " The indexClosestAdjBndryNodeToNodeOne is " << indexClosestAdjBndryNodeToNodeOne << std::endl;
   // return 0;

  // determine which edge is minimal :: the one in scan 1, or scan 2
  std::cout << "going to find indices of closest adjancet boundary nodes to (v_1,v_2)" << std::endl;
  Eigen::VectorXd closestAdjBndryNodeToNodeOne = boundaryVertices_scan1.row(indexClosestAdjBndryNodeToNodeOne);
  Eigen::VectorXd closestAdjBndryNodeToNodeTwo = boundaryVertices_scan2.row(indexClosestAdjBndryNodeToNodeTwo);

  Eigen::VectorXd nodeOne = boundaryVertices_scan1.row(scan1SeedPointIndex);
  Eigen::VectorXd nodeTwo = boundaryVertices_scan2.row(scan2ClosestPointToSeedIndex);

  std::cout << " Found closest Adj Bndry Nodes, to (v_1,v_2) " << std::endl;

  int indexOfClosestPoint = -1;
  double scan1Distance =  (nodeOne - closestAdjBndryNodeToNodeOne).norm();
  double scan2Distance =  (nodeTwo - closestAdjBndryNodeToNodeTwo).norm();
  //std::cout << scan1Distance << '\t' << scan2Distance << std::endl;
  bool isItEdgeInScanOne = (scan1Distance < scan2Distance);

  std::cout << "Asserted distnace-norm difference on the scans for (v_1,v_2)" << std::endl;
  
  Eigen::Vector2i newEdge;
  Eigen::Vector3i newFace; // #TODO :: ensure that this is correct! 
  std::vector<int> newTriangleFaces;

  // construct new edge
  if (isItEdgeInScanOne) {
      indexOfClosestPoint = indexClosestAdjBndryNodeToNodeOne;
      newEdge = Eigen::Vector2i(seedEdge(0), indexOfClosestPoint);
 	  newFace = Eigen::Vector3i(seedEdge(1),newEdge(0), newEdge(1));
  }
  else {  
     indexOfClosestPoint = indexClosestAdjBndryNodeToNodeTwo;
     // note need to offset here
     newEdge = Eigen::Vector2i(seedEdge(1), (indexOfClosestPoint + scan1.V.rows()));
 	 newFace = Eigen::Vector3i(seedEdge(0),newEdge(0), newEdge(1));
  }

  // now that new edge is added, continue on with the algorithm ! 
  //std::cout << newEdge << std::endl;
  allEdges.push_back(newEdge); 
  newTriangleFaces.push_back(newFace(0));
  newTriangleFaces.push_back(newFace(1));
  newTriangleFaces.push_back(newFace(2));

  // convert the set of (3*faces) integers, of vertex indices, to a matrix ( for faces data )
  int numOfFaces = newTriangleFaces.size() / 3;
  const Eigen::MatrixXi interpolatedFaces = Eigen::Map<const Eigen::MatrixXi> (&newTriangleFaces[0],numOfFaces,3); //,RowMajor); ( too include ?? not sure ?? ) 
  std::cout << "faces are " << std::endl;
  std::cout << interpolatedFaces << std::endl;
  //return 0;



  // WHERE OLD ALGORITHM USED TO BE ( sectioned off to end of code)

  // CREATE ONE HUGE MESH containing the two partial scans and interpolated surface
  igl::cat(1,scan1.V,scan2.V,interpolatedSurface.V); //# doubel check if this is to be done!
  //igl::cat(1,scan1.V,scan2.V,scans.V);
  igl::cat(1,scans.V,interpolatedSurface.V,scene.V); 
  //igl::cat(1,scan1.V,scan2.V,scene.V);

  //interpolatedSurface.F = Eigen::MatrixXi::Zero(totalNumBoundaryVertices,3);   // #TODO :: update this l8r 
  igl::cat(1,scan1.F, MatrixXi(scan2.F.array() + scan1.V.rows()), scans.F);
  //std::cout << scans.F << std::endl;
  //std::cout << interpolatedFaces << std::endl;
  interpolatedSurface.F = interpolatedFaces;
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
          numBoundaryVertices++; return numBoundaryVertices;
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

Eigen::MatrixXd retrieveBoundaryVerticesInScan1(int bndryVerticesIdxs[], int numBndryVertices)
{
  Eigen::MatrixXd bndryVertices = Eigen::MatrixXd::Zero(numBndryVertices,3);  
  int i;
  for ( i = 0; i < numBndryVertices; i++)
  {
     int indexIntoVertices = bndryVerticesIdxs[i];
     if ( indexIntoVertices != -1 ) {
		bndryVertices.row(i) = (scan1.V).row(indexIntoVertices); 
     }
  }
  return bndryVertices; 
}

Eigen::MatrixXd retrieveBoundaryVerticesInScan2(int bndryVerticesIdxs[], int numBndryVertices)
{
  Eigen::MatrixXd bndryVertices = Eigen::MatrixXd::Zero(numBndryVertices,3);  
  int i;
  for ( i = 0; i < numBndryVertices; i++)
  {
     int indexIntoVertices = bndryVerticesIdxs[i];
     if ( indexIntoVertices != -1 ) {
		bndryVertices.row(i) = (scan2.V).row(indexIntoVertices); 
     }
  }
  return bndryVertices; 
}



std::vector<int> findAdjBndryVertsInScan1(int vertIdx)
{
  const std::vector<int>& localAdjListToVertex = Adjacency_Scan1[vertIdx];
  std::vector<int> bndryVertsToTest; 
  for ( int i = 0; i < localAdjListToVertex.size(); i++)
  {
      int bndryVertIdx = localAdjListToVertex[i];
      if (boundaryVerticesStatus_scan1[bndryVertIdx] == 1) {
			bndryVertsToTest.push_back(bndryVertIdx);
      }
      
  }            
  return bndryVertsToTest;
}

std::vector<int> findAdjBndryVertsInScan2(int vertIdx)
{
  const std::vector<int>& localAdjListToVertex = Adjacency_Scan2[vertIdx];
  std::vector<int> bndryVertsToTest; 
  std::cout << " Found adj list " << std::endl;
  printcoll(localAdjListToVertex); // this is a seg fault ... but why?
  for ( int i = 0; i < localAdjListToVertex.size(); i++)
  {
      int bndryVertIdx = localAdjListToVertex[i];
      if (boundaryVerticesStatus_scan2[bndryVertIdx] == 1) {
			bndryVertsToTest.push_back(bndryVertIdx);
      }
      
  }            
  std::cout << " Found adj bndry vertices " << std::endl;
  return bndryVertsToTest;
}

// note :: does adjacent list, list vertices to be self-adjacent? I hope not ..., else I get a bug here!
// note :: I assume that "myVertexIndex" is guaranteed to be an index to a boundary Vertex!
// note :: I lvoe the cool ( vectors are essentially arrays ) trick! 
int findClosestAdjBndryNodeInScan1(std::vector<int> bndryVertices, int myVertexIndex)
{
    int closestAdjBndryPoint = -1; 
    Eigen::MatrixXd closestAdjBoundaryPoint;
    Eigen::MatrixXd myVertex = scan1.V.row(myVertexIndex);
	Eigen::MatrixXd adjBndryVertices = retrieveBoundaryVerticesInScan1(&bndryVertices[0], bndryVertices.size());

    //std::cout << myVertex << std::endl;
    //std::cout << adjBndryVertices << std::endl; ( you index into this matrix ! ) 

	Eigen::VectorXi Ele = Eigen::VectorXi::LinSpaced(adjBndryVertices.rows(), 0, adjBndryVertices.rows() - 1);
	Eigen::VectorXd smallestSquaredDists;
	Eigen::VectorXi smallestDistIndxs;
	igl::point_mesh_squared_distance(myVertex,adjBndryVertices,
                                 Ele,
                  				smallestSquaredDists,smallestDistIndxs,
                  				closestAdjBoundaryPoint);
    //std::cout << "smallest Dist indx = " << smallestDistIndxs << std::endl;
    //closestAdjBndryPoint = smallestDistIndxs(0,0); 
    closestAdjBndryPoint = bndryVertices[smallestDistIndxs(0,0)]; 
    return closestAdjBndryPoint; 
}

int findClosestAdjBndryNodeInScan2(std::vector<int> bndryVertices, int myVertexIndex)
{
    int closestAdjBndryPoint = -1; 
    Eigen::MatrixXd closestAdjBoundaryPoint;
    Eigen::MatrixXd myVertex = scan2.V.row(myVertexIndex);
	Eigen::MatrixXd adjBndryVertices = retrieveBoundaryVerticesInScan1(&bndryVertices[0], bndryVertices.size());
	Eigen::VectorXi Ele = Eigen::VectorXi::LinSpaced(adjBndryVertices.rows(), 0, adjBndryVertices.rows() - 1);
	Eigen::VectorXd smallestSquaredDists;
	Eigen::VectorXi smallestDistIndxs;
	igl::point_mesh_squared_distance(myVertex,adjBndryVertices,
                                 Ele,
                  				smallestSquaredDists,smallestDistIndxs,
                  				closestAdjBoundaryPoint);
    closestAdjBndryPoint = smallestDistIndxs(0,0); 
    return closestAdjBndryPoint; 
}


template <typename T>
void printcoll (T const& coll)
{
    typename T::const_iterator pos;  // iterator to iterate over coll
    typename T::const_iterator end(coll.end());  // end position

    for (pos=coll.begin(); pos!=end; ++pos) {
        std::cout << *pos << ' ';
    }
    std::cout << std::endl;
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

// useful method, to print elements of any STL container 
// adapted from website http://www.java2s.com/Tutorial/Cpp/0260__template/templatefunctiontoprintelementsofanSTLcontainer.htm

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
