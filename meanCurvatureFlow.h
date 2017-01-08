#ifndef MCF_H
#define MCF_H

// check if all these include files are actually needed?
// also ... do u put include files in header, or the actual code base?
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
#include <igl/is_boundary_edge.h> 
#include <igl/on_boundary.h> 
#include <igl/is_border_vertex.h>
#include <set>

using namespace Eigen; 
using namespace std;

//Eigen::MatrixXd V_one;
//Eigen::MatrixXi F_one;

extern Eigen::MatrixXd V_mcf;
extern Eigen::MatrixXd F_mcf;
extern Eigen::SparseMatrix<double> massMatrix_iterK ;
extern Eigen::SparseMatrix<double> stiffnessMatrix_iterK ; 
const double delta = 0.001;
//extern int k = 0; 					// used to print # iters of MCF, over timesteps ( &key_down  )

bool meshHasBoundary(const Eigen::MatrixXd& V, const Eigen::MatrixXi &F);
//void applyOneTimeStepOfMcfBoundaryCase(const Eigen::MatrixXd& V, const Eigen::MatrixXi &F);
void applyOneTimeStepOfMcfBoundaryCase();
void applyOneTimeStepOfMcfWatertightCase(const Eigen::MatrixXd& V, const Eigen::MatrixXi &F);
void applyMeanCurvatureFlow(const Eigen::MatrixXd& V, const Eigen::MatrixXi &F);
//bool key_down( igl::viewer::Viewer& , unsigned char , int );
#endif
