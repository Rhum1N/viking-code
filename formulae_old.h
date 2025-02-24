#include <iostream>
#include <cstdlib>
#include <math.h>
#include <Eigen/Core>
#include <vector>

using namespace std;
using namespace Eigen;


//Computes the next position of our submarine
Vector4f nextPos(Vector4f X, float angle, float timestep);



//BERNSTEIN RELATED

//Gives the value of P(x) where P is the Berstein polynomial in matrix form
double polyValue(float** P, int d, Vector4f x, int a, int b, int e, int f);

double density_value(vector<vector<float>> density, int d, float x1, float x2, int index, int a, int b, int e, int f, float sigma_theta);

void printMat(float** mat, int size);
void printVecindex(vector<vector<float>> vec, int index, int size);


// Compute the approximation of h and then compute (theta - h)^2 on [a,b]*[e,f]
void bern(float** coeffCarre, float** coeff,Vector4f z, float theta, int d, int a, int b, int e, int f, float sigma);

void printVec(vector<vector<float>> vec);

//Add the initial distribution of X to the density
void addfv(float** coeffcarre, float mu1, float mu2, float sigma1, float sigma2, float sigma_theta);

void adddens(float** coeffcarre, vector<vector<float>> dens, int index);

double F(int n, VectorXf theta, float i1, float i2, Vector4f* z);


double normalisation(float** density, int riemann, int d, float sigma_theta);

double normalization(vector<vector<float>> den, int riemann, int d, int index, float sigma_theta);
//computes the reward 
float reward(Vector4f z, float x1, float x2, float action);

//Comute the reward with density
double reward_density(Vector4f z, vector<vector<float>> density, vector<double> norm, int d, float action, int index,float sigma_theta);


//computes the transition probabaility of going to density j from density i
double transition_theta(float theta, Vector4f z, int d, float a, vector<vector<float>> density, vector<double> norms, int index,float sigma_theta);

float h(float zx1, float zx2);