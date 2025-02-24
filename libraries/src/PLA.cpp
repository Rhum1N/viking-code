#include "PLA.h"
#include "formulae.h"
#include "struc.h"
#include <math.h>


unsigned seed2 = chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator3(seed2);
uniform_real_distribution<double> distribution3(0.0, 1.0); 


const int d = 1; //approximation parameter
const float sigma = 1;
const float sigma_theta = 0.2;

VectorXf discretizedTheta = VectorXf::LinSpaced(n_theta, -PI, PI);

template <typename T, typename A>
int arg_max(std::vector<T, A> const& vec) {
	return static_cast<int>(std::distance(vec.begin(), max_element(vec.begin(), vec.end())));
}

template <typename T, typename A>
int arg_min(std::vector<T, A> const& vec) {
	return static_cast<int>(std::distance(vec.begin(), min_element(vec.begin(), vec.end())));
}

int index_theta(float theta) {
	int index = 0;
	float best_value = abs(theta - discretizedTheta[0]);
	for (int i = 1; i < n_theta; i++) {
		if (abs(theta - discretizedTheta[i]) < best_value) {
			best_value = abs(theta - discretizedTheta[i]);
			index = i;
		}
	}
	return index;
}

//generate a sample from a probability vector
int sample_action(vector<float> probaMat,int index) {
	int a = 1;
	float randomnum = distribution3(generator3);
	float cumul = probaMat[index*n_A];
	while (randomnum >= cumul) {
		if (a >= n_A) {
			return n_A-1;
		}
		cumul += probaMat[a-1 +index*n_A];
		a++;
	}
	return a - 1;
}

//sample a policy from the probability matrix stored in the tree (policy = a1,a2,...)
vector<int> sample_policy(tree t) {
	vector<int> policy(H, 0);
	return policy;
}
/*

float phiproba(tree t, policy p,int h) {
	if (h == H-1) {
		return 1;
	} 
	float value = 1;
	int cardX = pow(n_theta, h);
	for (int i = 0; i < cardX; i++) {
		value = value * t.probaMatrix[i * (n_A)+p.actions[i]];
	}
	for (int a = 0; a < n_A; a++) {
		value = value * phiproba(*t.childs[a], *p.childs[a], h + 1);
	}
	return value;

}
*/

void PLA(int i, tree* t, int index, int N, float mu) {
	/*
	INPUT : stage i
			state density (index)
	*/
	int ak;
	int a_opti;
	float wk;
	float V;
	vector<float> Ma(n_A, 0);
	vector<int> Na(n_A, 0);
	vector<float> Qa(n_A, 0);

	if (i == H) {
		V = 0;
		return;
	}
	for (int k = 0; k < N; k++) {
		/*
		Sample a(k), wk
		upadate Q estimate for a(k)
		M = M(x,a(k)) + R(x,a(k),wk) + pla(i+1)
		Na(k)(x) = Na(k) +1
		Q = M(x,a(k))/Na(k)

		Update optimal action
		a_opti = argmax Q(x,a) (on regarde si le nouveau est meilleur que l'ancien

		Update proba distribution over action space

		*/

		ak = sample_action(t->probaMatrix,index);
		float wk = distribution3(generator3);
		Ma[ak] += reward(t->position, 0, 0, ak); //+ PLA(zefzf) On simule le theta aussi 
		Na[ak] += 1;
		Qa[ak] = Ma[ak] / Na[ak];
		a_opti = arg_max(Qa);

	}

	return;
}


/*
float path_value(policy* p, tree* t,float theta) {
	int action = p->actions[0];
	float value = reward_density(t->position, t->densities, t->norms, d, action, 0, sigma_theta);
	float theta_rand;
	int index = 0;
	for (int h = 1; h < H; h++) {
		theta_rand = distribution3(generator3);
		t = t->childs[action];
		p = p->childs[action];
		index = index_theta(theta_rand);
		action = p->actions[index];
		value += reward_density(t->position, t->densities, t->norms, d, action, 0, sigma_theta);
	}
	return 0.0f;
}

void ASA(tree* t, vector<float> alpha, vector<float> beta, int Nk, int Mk) {
	
	float randfloat;
	vector<policy*> sample_policies;
	vector<vector<float>> path;
	for (int k=0; k < 10; k++) { //Stopping criterion ????
		for (int n = 0; n < Nk; n++) {
			randfloat = distribution3(generator3);
			if (1 - beta[k] > randfloat) {
				cout << "1-B" << endl;
				sample_policies.push_back(new policy(0));
			}
			else {
				cout << "B" << endl;
				sample_policies.push_back(new policy(0, t));
			}
			for (int j = 0; j < Mk; j++) {
				
			}
		}

		

	}
	
	return;
}
*/
void hello() {
	cout << "Hello Word" << endl;
	vector<float> test = { 0.05,0.02,0.9,0.03 };
	int a = 0;
	for (int i = 0; i <= 10; i++) {
		a = sample_action(test,0);
		cout << a << " - ";
	}
	vector<float> dummy(100,0.1);
	
	return;
}
