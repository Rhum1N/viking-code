

#include "formulae.h"
#define PI 3.141592
const int d = 1;

const int A = -2;
const int B = 2;

const float vx = 0;
const float vy = 0;

void printVec(vector<vector<float>> vec) {
	for (int i = 0; i < vec.size(); i++) {
		for (int j = 0; j < vec[0].size(); j++) {
			cout << vec[i][j] << ", ";
		}
		cout << endl;
	}
}

//Atan2 function prolonged in 0 ??????
float h(float zx1, float zx2) {
	if ((zx1 == 0) && (zx2 == 0))
		return atan2(0, 0.001);
	return atan2(zx1, zx2);
}

//Next position given the action and a timestep
Vector4f nextPos(Vector4f X, float angle, float timestep) {
	Matrix4f T;
	Vector4f next;
	if (angle != 0.0) {
		T << 1, sin(angle * timestep) / angle, 0, -(1 - cos(angle * timestep)) / angle,
			0, cos(angle * timestep), 0, -sin(angle * timestep),
			0, (1 - cos(angle * timestep)) / angle, 1, sin(angle * timestep) / angle,
			0, sin(angle * timestep), 0, cos(angle * timestep);
	}
	else {
		T << 1, timestep, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, timestep,
			0, 0, 0, 1;
	}
	next = T * X;

	next[0] -= vx * timestep;
	next[2] -= vy * timestep;

	return next;
}

void addfv(float** coeffcarre, float mu1, float mu2, float sigma1, float sigma2, float sigma_theta) {
	float s1 = 1 / (sigma1 * sigma1);
	float s2 = 1 / (sigma2 * sigma2);
	coeffcarre[2][0] = s1 ; // x_1² coeff
	coeffcarre[0][2]  = s1 ; // x_2² coeff
	coeffcarre[1][0] = -2 * mu1 /(sigma1*sigma1) ;
	coeffcarre[0][1] = -2 * mu2 /(sigma2*sigma2) ;
	coeffcarre[0][0] = (mu1 * mu1 * s1 + mu2 * mu2 * s2) ;
	return;
}

void adddens(float** coeffcarre, vector<vector<float>> dens, int index) {
	for (int i = 0; i < 2 * d + 1; i++) {
		for (int j = 0; j < 2 * d + 1; j++) {
			coeffcarre[i][j] = dens[i][j + index*(2*d+1)];
		}
	}
	return;
}


//Normal density with mean mu and standard deviation sigma
double normal(double x, double mu, double sigma) {
	double res = -pow((x - mu) / sigma, 2) / 2;
	return exp(res) / (sigma * sqrt(2 * PI));
}

double transition_theta(float theta, Vector4f z, int d, float a, vector<vector<float>> density, vector<double> norms, int index, float sigma_theta) {
	double res = 0;
	int riemann = 100;
	float theta2 = 0;
	float x0 = 0;
	float x2 = 0;
	float tran;
	for (int i = 0; i < riemann; i++) {
		for (int j = 0; j < riemann; j++) {
			x0 = (B-A) * ((float)i / riemann) +A;
			x2 = (B-A) * ((float)j / riemann) +A;
			theta2 = h(x0 - z(0), x2 - z(2));
			tran = normal(theta - theta2, 0, sigma_theta);
			res += tran * density_value(density, d, x0, x2, index, A,B,A,B, sigma_theta) / norms[index];
		}
	}
	return (B-A)*(B-A)*res / (riemann * riemann);
}

//Binomial coefficient
int comb(int n, int k) {
	if ((n == k) || (k == 0))
		return 1;
	float prod = 1;

	for (int i = 0; i < k; i++) {
		prod *= (float)(n - i + 0.0) / (k - i + 0.0);
	}
	return prod + 0.5;
}

//performs transoformation of intervall x in [c,d] into x in [0,1] 
double phi(double x, double a, double b) {
	return (x - a) / (b - a);
}

//performs the inverse transformation
double invPhi(double x, double a, double b) {
	return x * (b - a) + a;
}


float reward(Vector4f z, float x1, float x2, float action) {
	float r = 0;
	Vector4f next = nextPos(z, action, 1);
	float value = abs(h(x2 - next(2), x1 - next(0)) - h(x2 - z(2), x1 - z(0)));
	return min(value, 2 * (float)PI - value);
}



//integral of the reward against the density
double reward_density(Vector4f z, vector<vector<float>> density, vector<double> norm, int d, float action, int index, float sigma_theta)
{
	int riemann = 100;
	double r = 0;
	for (int i = 0; i < riemann; i++) {
		for (int j = 0; j < riemann; j++) {
			r += reward(z, (B-A) * ((float)i / riemann) +A, (B-A) * ((float)j / riemann) +A, action) * (density_value(density, d, i, j, index, 0, riemann, 0, riemann, sigma_theta) / norm[index]);
		}
	}
	return (B-A)*(B-A)*r / (riemann * riemann);
}


double density_value_pointer(float** density, int d, float x1, float x2, int a, int b, int e, int f, float sigma_theta) {
	double res = 0;
	//float r = sqrt(pow(x(0), 2) + pow(x(2), 2));
	for (int i = 0; i < 2 * d + 1; i++) {
		for (int j = 0; j < 2 * d + 1; j++) {

			res += density[i][j] * pow((B-A) * phi(x1, a, b) + A, i) * pow((B-A) * phi(x2, e, f) +A, j);

		}
	}
	return exp(-res/2);
}

double density_value(vector<vector<float>> density, int d, float x1, float x2, int index, int a, int b, int e, int f,float sigma_theta) {
	double res = 0;
	//float r = sqrt(pow(x(0), 2) + pow(x(2), 2));
	double r = 1;
	for (int i = 0; i < 2 * d + 1; i++) {
		for (int j = 0; j < 2 * d + 1; j++) {
			res += density[i][index * (2 * d + 1) + j] * pow((B-A) * phi(x1 / r, a, b) +A, i) * pow((B-A) * phi(x2 / r, e, f) +A, j);

		}
	}
	return exp(-res / 2);
}

//normalise the numerator give true 2*d
double normalisation(float** density, int riemann, int d, float sigma_theta) {
	double n = 0;
	for (int i = 0; i < riemann; i++) {
		for (int j = 0; j < riemann; j++) {
			n += density_value_pointer(density, d, i, j, 0, riemann, 0, riemann, sigma_theta);
		}
	}

	return (B-A)*(B-A)*n / (riemann * riemann);
}


//normalise the numerator give true 2*d
double normalization(vector<vector<float>> den, int riemann, int d, int index, float sigma_theta) {
	double n = 0;
	for (int i = 0; i < riemann; i++) {
		for (int j = 0; j < riemann; j++) {
			n += density_value(den, d, i, j, index, 0, riemann, 0, riemann, sigma_theta);
		}
	}

	return (B - A) * (B - A)*n / (riemann * riemann);
}




void coefBernHextended(float** coef, Vector4f z, int d, int a, int b, int e, int f) {
	float coefVal;
	int i1;
	int i2;
	for (int i = 0; i < d + 1; ++i) {
		for (int j = 0; j < d + 1; j++) {
			coefVal = 0;
			for (i1 = 0; i1 <= i; i1++) {
				for (i2 = 0; i2 <= j; i2++) {

					coefVal += h(z(0) - ((float)i1 / d) * (b - a) + (b - a) / 2, z(2) - ((float)i2 / d) * (f - e) + (f - e) / 2) * comb(d, i1) * comb(d, i2) * comb(d - i1, i - i1) * comb(d - i2, j - i2) * pow(-1, i + j - i1 - i2);
				}
			}
			coef[i][j] = coefVal;
		}

	}
	return;
}


double polyValue(float** mat, int d, Vector4f x, int a, int b, int e, int f) {
	double res = 0;
	//float r = sqrt(pow(x(0), 2) + pow(x(2), 2));
	double r = 1;
	for (int i = 0; i < d + 1; i++) {
		for (int j = 0; j < d + 1; j++) {
			res += mat[i][j] * pow(phi(x(0) / r, a, b), i) * pow(phi(x(2) / r, e, f), j);
		}
	}
	return res;
}

//hcoef is the approximation of h 
//Returns (theta - h)**2 
void carreCoef(float** coef, float** hcoef, float theta, int d,float sigma) {
	float coefVal;
	for (int i = 0; i < 2 * d + 1; i++) {
		for (int j = 0; j < 2 * d + 1; j++) {
			coefVal = 0;
			/*
			for (int k1 = 0; k1 <= i; k1++) {
				for (int k2 = 0; k2 <= j; k2++) {
					if ((k1 <= d) and (k2 <= d) and (i - k1 <= d) and (j - k2 <= d)) {
						coefVal += hcoef[k1][k2] * hcoef[i - k1][j - k2];
					}
				}
			}*/
			for (int k = max(0, i - d); k <= min(i, d); k++) {
				for (int l = max(0, j - d); l <= min(j, d); l++) {
					coefVal += hcoef[k][l] * hcoef[i - k][j - l];
				}
			}
			if ((i <= d) and (j <= d)) {
				coefVal -= 2 * theta * hcoef[i][j];
			}

			//cout << "i: " << i << " - j: " << j << endl;

			coef[i][j] = coefVal/(sigma*sigma);
		}
	}
	coef[0][0] += pow(theta, 2);
	return;
}

//Function that 
void bern(float** coeffCarre, float** coeff, Vector4f z, float theta, int d, int a, int b, int e, int f,float sigma) {
	coefBernHextended(coeff, z, d, a, b, e, f);
	carreCoef(coeffCarre, coeff, theta, d,sigma);
}




//F function
double F(int n, VectorXf theta, float i1, float i2, Vector4f* z) { //i1 = i1/d and i2 = i2/d in h formula = discretization of x.
	double f = 0;
	for (int i = 0; i < n; i++) {
		f += pow(theta(i) - h(i2 - z[i](2), i1 - z[i](0)), 2);
	}
	return f / (0.01);
}

