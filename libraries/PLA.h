#include <iostream>
#include <cstdlib>
#include <math.h>
#include <Eigen/Core>
#include <vector>
#include <random>
#include <chrono>


using namespace std;
using namespace Eigen;

void hello();

void PLA(int i, vector<vector<float>> density, int index, int N, float mu, int H);

void MRAS();

// generate a sample from a probability vector
int sample_action(vector<float> probaMat, int index);

