// purpose of this code :: given a (4x4) transformation matrix/alignemnt, iteratively improve the alignemnt, via stoichastic gradient descent !  try 100 of these things ... take min. it will most likely converge ! 
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


// other include files, at the moment
#include <igl/random_dir.h> // used to generate uniform rand unit dir in 3D ... ret as vector ... needs scaling though !  ... not really going to be used ( jsut focus on using \epsilon e_x,e_y,e_z values ! ) 
#include <igl/axis_angle_to_quat.h> // used to generate a sample rotation matrix 
#include <set>
#include <time.h> // used for rand val generator  ... maybe
//#include <random.h> // c++11 rand val generator  ( note :: do not actually use rand() ! ) 

using namespace Eigen; 
using namespace std;

// Q1 :: What is the order of convergence of this method? not sure atm

// Develop an initial transformation matrix, filled with random data 
Eigen::MatrixXd T;

// REFERNCE THIS FOR RANDOMIZER :: http://stackoverflow.com/questions/288739/generate-random-numbers-uniformly-over-an-entire-range/20136256#20136256
// generate 100 random (3x1) vector of small {e_x,e_y,e_z} values

/*
const double range_from  = -0.001;
const double range_to    = 0.001;
std::random_device                  rand_dev;
std::mt19937                        generator(rand_dev());
std::uniform_real_distribution<double>  distr(range_from, range_to);

int i;
for(i = 0; i < 100; i++)
{

	// generator values for an epsilon vector 
	double e_x = distr(generator);
 	double e_y = distr(generator);
 	double e_z = distr(generator);

    Eigen::Vector3d epsilon_vector;
	epsilon_vector[0] = e_x;
	epsilon_vector[1] = e_y;
	epsilon_vector[2] = e_z;

}
*/





/* initialize random seed: 
randVals = rand(3,1)*7 - 10;
randVals = 10.^randVals;
 * srand (time(NULL));
 * generate secret number between -7 and 2: 
 * int randVal = rand() % 10 - 8;
 * double epsilonVal = 10 ^ randVal;
 */








// generate perturbed ( Rotation, Translation ) matrices. 

// evaluate them, based on surface metric ! 




