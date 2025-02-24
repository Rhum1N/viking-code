#include <iostream>
#include <cstdlib>
#include <math.h>
#include <Eigen/Core>
#include <vector>
#include <iostream>
#include <limits>


using namespace std;
using namespace Eigen;


//Computes the next position of our submarine
Vector4f nextPos(Vector4f X, float angle, float timestep);

struct test_t {
	int r;
	int t;

	test_t(int a, int b);
	void pri();
};
//BERNSTEIN RELATED

//Gives the value of P(x) where P is the Berstein polynomial in matrix form
long double polyValue(float** mat, int d, Vector4f x);

long double density_value(vector<vector<float>> density, int d, float x1, float x2, int index, int a, int b, int e, int f, float sigma_theta);

void printMat(float** mat, int size);
void printVecindex(vector<vector<float>> vec, int index, int size);


// Compute the approximation of h and then compute (theta - h)^2 on [a,b]*[e,f]
void bern(float** coeffCarre, float** coeff,Vector4f z, float theta, int d, int a, int b, int e, int f, float sigma);

void printVec(vector<vector<float>> vec);

//Add the initial distribution of X to the density
void addfv(float** coeffcarre, float mu1, float mu2, float sigma1, float sigma2, float sigma_theta);

//add the density coefficient in dens at index to coeffcarre
void adddens(float** coeffcarre, vector<vector<float>> dens, int index);


long double Fun(float theta, float i1, float i2, Vector4f z, float sigma);

//Normalize the density, works with float**
long double normalisation(float** density, int riemann, int d, float sigma_theta);

//Normalize the density, works with vector<vector<flaot>> 
long double normalization(vector<vector<float>> den, int riemann, int d, int index, float sigma_theta);

//computes the reward 
float reward(Vector4f z, float x1, float x2, float action);

//Comute the reward with density
long double reward_density(Vector4f z, vector<vector<float>> density, vector<long double> norm, int d, float action, int index,float sigma_theta);


//computes the transition probabaility of going to density j from density i
long double transition_theta(float theta, Vector4f z, int d, float a, vector<vector<float>> density, vector<long double> norms, int index,float sigma_theta);

//computes h function 
float h(float zx1, float zx2);