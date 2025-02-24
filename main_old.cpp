
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <omp.h>
#include <Eigen/Core>
#include <chrono>
#include <random>
#include <vector>

#include "formulae.h"
#include "PLA.h"

using namespace std;
using namespace Eigen;

#define PI 3.14159265
#define H 3 //The horizon
#define PLA 1

//X space [-1,1]²
#define A -2
#define B 2
#define E -2
#define F 2

const int d = 1; //approximation parameter
const int n_A = 5; //Number of actions
const int n_theta = 10; //Number of thetas...
const float sigma = 0.6;
const float sigma_theta = 0.03;

//Time initialisation process
unsigned seed = chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator(seed);
normal_distribution<float> x_distrib(0.0,sigma);
normal_distribution<float> distribution(0.0, sigma_theta);
uniform_real_distribution<double> distribution2(0.0, 1.0);

//Action space :  n_a actions between pi/2 and -pi/2
VectorXf actions = VectorXf::LinSpaced(n_A, -2*PI / 5, 2*PI / 5);
VectorXf discretizedThetas = VectorXf::LinSpaced(n_theta, -PI, PI);

float** coeff;
float** coeffcarre;

template <typename T, typename M>
int arg_max(std::vector<T, M> const& vec) {
    return static_cast<int>(std::distance(vec.begin(), max_element(vec.begin(), vec.end())));
}

template <typename T, typename M>
int arg_min(std::vector<T, M> const& vec) {
    return static_cast<int>(std::distance(vec.begin(), min_element(vec.begin(), vec.end())));
}

int arg_max_ind(vector<float> M, int index, int size) {
    float max = 0;
    int arg = 0;
    for (int m = 0; m < size; m++) {
        cout << M[size * index + m] << " - ";
        if (M[size * index + m] > max) {
            max = M[size * index + m];
            arg = m;
        }
    }
    return arg;
}

//Tree structure 
struct tree {
    int root;
    Vector4f position; //position Z of the node 
    tree* childs[n_A];
    tree* previous; //pointer to the father of the node
    vector<vector<float>> densities; // matrix with all the densities
    vector<double> value_function; // vector containing the value function result and the action chosen 
    vector<double> norms; // vector containing the denominator of the density
    vector<float> probaMatrix; // vector of size n_A with the probabilities of chosing each action per state 

    tree(int a, Vector4f z) {
        root = a;
        position = z;

        for (int i = 0; i < n_A; i++) {
            childs[i] = NULL;
        }
    }
    
    //Tree creation for the root
    tree(int a, Vector4f z, int depth, float** coeffcarre, float** coeff) {
        root = a;
        previous = NULL;
        position = z;
        int power = pow(n_theta, H - depth);
        densities = vector<vector<float>>(2 * d + 1, vector<float>(power * (2 * d + 1), 0));//2*d+1 rows matrix with n_theta*2 columns (n_theta 2d+1*2d+1 matrices)
        value_function = vector<double>(2 * power, -INFINITY);
        probaMatrix = vector<float>(power*n_A, 1 / (float)n_A);
        //cout << 1 / (float)n_A;
        double norm = 0;
        //cout << z.transpose() << endl;
        //creation of the densities
        for (int i = 0; i < 2 * d + 1; i++) {
            for (int j = 0; j < 2 * d + 1; j++) {
                coeffcarre[i][j] = 0;
            }
        }
        addfv(coeffcarre, 0, 0,sigma,sigma,sigma_theta);
        norm = normalisation(coeffcarre, 100, d, sigma_theta);
	if (std::isinf(norm))
            cout << "ERREUR";
        norms.push_back(norm);
        //cout << norm << "is the norm" << endl;
        for (int i = 0; i < 2 * d + 1; i++) {
            for (int j = 0; j < 2 * d + 1; j++) {
                densities[i][j] = coeffcarre[i][j];
            }
        }
        printVec(densities);

        //cout << densities[0].size() << " test from " << root << endl;
        if (depth != 0) {
            for (int i = 0; i < n_A; i++) {  //recursion to create the next part of the tree
                Vector4f nextZ = nextPos(z, actions(i), 1);
                childs[i] = new tree(i, this, nextZ, depth - 1, coeffcarre, coeff);
            }
        }
        else {
            for (int i = 0; i < n_A; i++) {

                childs[i] = NULL;
            }
        }
    }
    
    //Tree creation depth >1, computes bernstein for the actual step and add the matrix to the previous densities (n_theta new densities for every previous densities)
    tree(int a, tree* prev, Vector4f z, int depth, float** coeffcarre, float** coeff) {
        root = a;
        previous = prev;
        int power = pow(n_theta, H - depth);
        position = z;
        densities = vector<vector<float>>(2 * d + 1, vector<float>(power * (2 * d + 1), 0));//2*d+1 rows matrix with n_theta*2 columns (n_theta 2d+1*2d+1 matrices)
        value_function = vector<double>(2 * power, -10);
        probaMatrix = vector<float>(power * n_A, 1 / (float)n_A);

        double norm = 0.0;
        //cout << z.transpose() << endl;

        //creation of the densities ADD NORM
        for (int teta = 0; teta < n_theta; teta++) { //for every new theta we compute the bernstein polynomial
            bern(coeffcarre, coeff, z, discretizedThetas[teta], d, A, B, E, F,sigma_theta);
            for (int prev_index = 0; prev_index < pow(n_theta, H - depth-1); prev_index++) { //for every previous density we add the new bernstein approximation to the previous densities

                for (int i = 0; i < 2 * d + 1; i++) {
                    for (int j = 0; j < 2 * d + 1; j++) {
                        densities[i][j + teta * (2 * d + 1) + prev_index * n_theta * (2 * d + 1)] = coeffcarre[i][j] + prev->densities[i][j + prev_index * (2 * d + 1)]; //+ l'ancienne densité
                    }
                }
            }
        }
        for (int l = 0; l < densities[0].size() / (2 * d + 1); l++) {
            norm = normalization(densities, 100, d, l,sigma_theta);
	    if (std::isinf(norm))
            	cout << "ERREUR";
            norms.push_back(norm);
            //cout << "reward " << reward_density(position, densities, norms, d, -PI / 2, l, sigma_theta) << " " << reward_density(position, densities, norms, d, 0, l, sigma_theta) << " " << reward_density(position, densities, norms, d, PI / 2, l, sigma_theta) << endl;
        }

        //cout << densities[0].size() << " test from " << root << " with father " << previous->root << endl;
        if (depth != 0) { //recursion to create the next part of the tree
            for (int i = 0; i < n_A; i++) {
                Vector4f nextZ = nextPos(z, actions(i), 1);
                childs[i] = new tree(i, this, nextZ, depth - 1, coeffcarre, coeff);
            }
        }
        else {
            for (int i = 0; i < n_A; i++) {

                childs[i] = NULL;
            }
        }

    }

    tree(int a, Vector4f z, vector<vector<float>> dens, int index ,int depth, float** coeffcarre, float** coeff) {
        root = a;
        previous = NULL;
        position = z;
        int power = pow(n_theta, H - depth);
        densities = vector<vector<float>>(2 * d + 1, vector<float>(power * (2 * d + 1), 0));//2*d+1 rows matrix with n_theta*2 columns (n_theta 2d+1*2d+1 matrices)
        value_function = vector<double>(2 * power, -INFINITY);
        probaMatrix = vector<float>(power * n_A, 1 / (float)n_A);

        double norm = 0;
        //cout << z.transpose() << endl;
        //creation of the densities
        adddens(coeffcarre,dens,index);
        norm = normalisation(coeffcarre, 100, d, sigma_theta);
        norms.push_back(norm);
        for (int i = 0; i < 2 * d + 1; i++) {
            for (int j = 0; j < 2 * d + 1; j++) {
                densities[i][j] = coeffcarre[i][j];
            }
            //cout << "reward " << reward_density(position, densities, norms, d, -PI / 2, teta) << " " << reward_density(position, densities, norms, d, 0, teta) << " " << reward_density(position, densities, norms, d, PI / 2, teta) << endl;
        }

        printVec(densities);


        //cout << densities[0].size() << " test from " << root << endl;
        if (depth != 0) {
            for (int i = 0; i < n_A; i++) {  //recursion to create the next part of the tree
                Vector4f nextZ = nextPos(z, actions(i), 1);
                childs[i] = new tree(i, this, nextZ, depth - 1, coeffcarre, coeff);
            }
        }
        else {
            for (int i = 0; i < n_A; i++) {

                childs[i] = NULL;
            }
        }
    }


};

struct policy {
    vector<int> actions;
    policy* childs[n_A];
    double phi;
    double a_n;
    double b_n;

    //Create policy from q_0 (uniform distribution)
    policy(int h) {
        actions = vector<int>(pow(n_theta, h), 0);
        for (int i = 0; i < pow(n_theta, h); i++) {
            actions[i] = rand() % n_A;
        }
        if (h == H-1) {
            for (int i = 0; i < n_A; i++) {
                childs[i] = NULL;
            }
        }
        else {
            for (int i = 0; i < n_A; i++) {
                childs[i] = new policy(h + 1);
            }
        }
    }
    //Create policy from probaMAt in treed
    policy(int h, tree* t) {
            
        actions = vector<int>(pow(n_theta , h), 0);


        for (int i = 0; i < pow(n_theta, h); i++) {
            actions[i] = sample_action(t->probaMatrix, i);
        }
        if (h == H-1 ) {
            for (int i = 0; i < n_A; i++) {
                childs[i] = NULL;
            }
        }
        else {
            for (int i = 0; i < n_A; i++) {
                childs[i] = new policy(h + 1,t->childs[i]);
            }
        }
    }

    void printpolicy() {
        if (this == NULL) {
            return;
        }
        cout << this->actions[0] << endl;
        this->childs[0]->printpolicy();
    }
};

void print(tree* t) {
    if (t == NULL) {
        //cout << "NULL : ";
        return;
    }
    //cout << t->position << " : ";
    cout << t->value_function[0] << " ; ";
    for (int i = 0; i < n_A; i++) {
        print(t->childs[i]);
    }
}

void printActionD(tree* t, int d) {
    if (t == NULL) {
        //cout << "NULL : ";
        return;
    }
    if (d == 0) {
        cout << "position is " << t->position << endl;
        for (int i = 0; i < t->value_function.size() / 2; i++)
            cout << t->value_function[i * 2 + 1] << ", " << t->value_function[i * 2] << " - ";

        cout << endl;
        return;
    }
    for (int i = 0; i < n_A; i++) {
        printActionD(t->childs[i], d - 1);
    }
}


void printDepthD(tree* t, int d) {
    if (t == NULL) {
        //cout << "NULL : ";
        return;
    }
    if (d == 0) {
        cout << t->position << " : ";

        return;
    }
    for (int i = 0; i < n_A; i++) {
        printDepthD(t->childs[i], d - 1);
    }
}

void freeTree(tree* t) {
    if (t == NULL)
        return;
    for (int i = 0; i < n_A; i++)
        freeTree(t->childs[i]);
    //free(&(t->densities));
    free(t);
    return;
}


//Find the nearest discretized theta
int nearest_theta(float theta) {
    int index = 0;
    float best_value = abs(theta - discretizedThetas[0]);
    for (int i = 1; i < n_theta; i++) {
        if (abs(theta - discretizedThetas[i]) < best_value) {
            best_value = abs(theta - discretizedThetas[i]);
            index = i;
        }
    }
    return index;
}


/// Bellman function.
void bellman(tree* t, int n) { //n is the depth 
    float a;
    double value;
    double transi;
    double bestValue;
    double best_action = 0;
    int indexNextDensity;
    vector<double> value_per_action = vector<double>(n_A, 0);
    int n_densities = t->densities[0].size() / (2 * d + 1);
    for (int j = 0; j < n_densities; j++) { //do the loop for every density in our node
        bestValue = -196848468165168516;
        for (int i = 0; i < n_A; i++) { //compute the value for every action in order to find the min
            a = actions(i);
            value_per_action[i] = reward_density(t->position, t->densities, t->norms, d, a, j,sigma_theta); // reward at position t->position and action a;
            
            if (n != H) { // n == H means that we are on a leaf (no children)
                for (int k = 0; k < n_theta; k++) { //Loop for the k densities we can have from this one
                    transi = transition_theta(discretizedThetas[k], t->position, d, a, t->densities, t->norms, j, sigma_theta);
                    if (transi>0.01)
                        value_per_action[i] += transi * t->childs[i]->value_function[2 * j + 2 * k];
                }
            } 
            if (value_per_action[i] > bestValue) { // we only keep the best value 
                bestValue = value_per_action[i];
                best_action = i;
            }
        }
        t->value_function[2 * j] = bestValue;
        t->value_function[2 * j + 1] = best_action;
    }
    return;
}

//computes bellman on the node of depth d in the tree
void bellmanDepthD(tree* t, int d, int n) {
    if (t == NULL) {
        return;
    }
    if (d == 0) {
        bellman(t, n);
        return;
    }
#pragma omp parallel for
    for (int i = 0; i < n_A; i++) {
        bellmanDepthD(t->childs[i], d - 1, n);
    }
}

void printMat(float** mat, int size) {
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            cout << mat[i][j] << ' ';
        }
        cout << endl;
    }
}
void printDensity(tree* t, int index, int d) {
    for (int i = 0; i < 2 * d + 1; i++) {
        for (int j = 0; j < 2 * d + 1; j++) {
            cout << t->densities[i][index * (2 * d + 1) + j] << " ";
        }
        cout << endl;
    }
}

//Routine créant l'arbre à partie de la distribution initiale 
void routine(tree** t, Vector4f z, float** coeffcarre, float** coeff,vector<vector<float>> dens,int index) {

    auto start = chrono::steady_clock::now();
    //creation of the tree
    if (index == -1)
        *t = new tree(0, z, H, coeffcarre, coeff);
    else
        *t = new tree(0, z, dens, index, H, coeffcarre, coeff);

    auto end = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    //cout << "elapsed time for tree creation: " << elapsed_seconds.count() << "s\n";


    //Main programm
    start = chrono::steady_clock::now();

    for (int h = 0; h <= H; h++) {
        //get index des trucs de profondeur H-h et calcul de bellman dessus


        //cout << "iteration " << h << " with H-h = " << H - h << endl;
        bellmanDepthD(*t, H - h, H - h);
    }

    end = chrono::steady_clock::now();
    elapsed_seconds = end - start;
    //cout << "elapsed time for bellman: " << elapsed_seconds.count() << "s\n";

}

float path_value(policy* p, tree* t, float theta) {
    int action = p->actions[0];
    normal_distribution<float> normal(0, sigma_theta);
    float theta_rand;
    float value = reward_density(t->position, t->densities, t->norms, d, actions[action], 0, sigma_theta);
    int index_theta = 0;
    int index_value;
    for (int h = 1; h < H; h++) {
        t = t->childs[action];
        theta_rand = atan2(0 - t->position[2], 0 - t->position[0])  +normal(generator);
        p = p->childs[action];

        index_theta = nearest_theta(theta_rand);
       
        if (h == 1) {
            index_value = index_theta;
        }
        else if (index_value != 0) {
            index_value = (index_value - 1) * n_theta + index_theta;
        }
        else {
            index_value = index_theta;
        }

        action = p->actions[index_value];
        value += reward_density(t->position, t->densities, t->norms, d, actions[action], index_value, sigma_theta);
    }
    return value;
}

//Compute phi(pi) used in ASA to update the tree
double phiproba(tree *t, policy *p, int h) {
    if (h == H) {
        return 1;
    }
    double value = 1;
    int cardX = pow(n_theta, h);
    for (int i = 0; i < cardX; i++) {

            value = value * t->probaMatrix[i * (n_A)+p->actions[i]];
        
    }
    for (int a = 0; a < n_A; a++) {
        value = value * phiproba(t->childs[a], p->childs[a], h + 1);
    }
    return value;
}

double phizero(int h) {
    if (h == H) {
        return 1;
    }
    double value = 1;
    int cardX = pow(n_theta, h);
    for (int i = 0; i < cardX; i++) {

        value = value/n_A;

    }
    for (int a = 0; a < n_A; a++) {
        value = value * phizero(h + 1);
    }
    return value;

}


double log_phizero(int h) {
    if (h == H) {
        return 1;
    }
    double value = 1;
    int cardX = pow(n_theta, h);
    for (int i = 0; i < cardX; i++) {

        value = value - log(n_A);

    }
    for (int a = 0; a < n_A; a++) {
        value = value + log_phizero(h + 1);
    }
    return value;

}


//Compute the log of phi for numerical purposes
double log_phiproba(tree* t, policy* p, int h) {
    if (h == H) {
        return 0;
    }
    double value = 0;
    int cardX = pow(n_theta, h);
    for (int i = 0; i < cardX; i++) {

        value = value + log(t->probaMatrix[i * (n_A)+p->actions[i]]);

    }
    for (int a = 0; a < n_A; a++) {
        value = value + log_phiproba(t->childs[a], p->childs[a], h + 1);
    }
    return value;
}

//Compute g_{k+1} to update the proba matrix
void boltzmann(vector<policy*> p, tree* t, int k, vector<float> Vk, float Tk, double phi_denominator) {
    float new_proba;
    if (k == H)
        return;
    for (int i = 0; i < pow(n_theta, k); i++) {
        for (int a = 0; a < n_A; a++) {
            new_proba = 0;
            for (int l = 0; l < p.size(); l++) {
                if (p[l]->actions[i] == a)
                    new_proba += exp(-Vk[l] / Tk) / (p[l]->phi * phi_denominator);
            }
            
            t->probaMatrix[i * n_A + a] =pow(0.01,0.501)* new_proba + (1- pow(0.01, 0.501))* t->probaMatrix[i * n_A + a];
        }
    }
    vector<vector<policy*>> temp(n_A,vector<policy*>(p.size(),NULL));
 
    for (int a = 0; a < n_A; a++) {
        for (int l = 0; l < p.size(); l++) {
            p[l]->childs[a]->phi = p[l]->phi;
            temp[a][l] = p[l]->childs[a];
        }
        boltzmann(temp[a], t->childs[a], k + 1, Vk, Tk, phi_denominator);
    }
    return;
}

//Compute g_{k+1} to update the proba matrix
void log_boltzmann(vector<policy*> p, tree* t, int k,int h, vector<float> Vk, float Tk,double max_an,double max_bn) {
    float new_proba;
    double log_denominateur=0;
    if (h == H)
        return;
    for (int l = 0; l < p.size(); l++) {
        log_denominateur += p[l]->a_n * exp(p[l]->b_n - max_bn) / max_an;
    }
   
    for (int i = 0; i < pow(n_theta, h); i++) {
        for (int a = 0; a < n_A; a++) {
            new_proba = 0;
            for (int l = 0; l < p.size(); l++) {
                if (p[l]->actions[i] == a)
                    new_proba += p[l]->a_n * exp(p[l]->b_n - max_bn) / max_an;
            }
            t->probaMatrix[i * n_A + a] = 0.1 * new_proba/log_denominateur + (1 -0.1) * t->probaMatrix[i * n_A + a];
        }
    }
    vector<vector<policy*>> temp(n_A, vector<policy*>(p.size(), NULL));

    for (int a = 0; a < n_A; a++) {
        for (int l = 0; l < p.size(); l++) {
            if (h < H - 1) {

                p[l]->childs[a]->a_n = p[l]->a_n;
                p[l]->childs[a]->b_n = p[l]->b_n;
                temp[a][l] = p[l]->childs[a];
            }
        }
        log_boltzmann(temp[a], t->childs[a],k, h + 1, Vk, Tk, max_an,max_bn);
    }
    return;
}



//compute the ASA algorithm
vector<int> ASA(tree* t, int Nk, int Mk, float theta,int* index) {

    float randfloat;
    float temp_mean = 0;
    float Tk=1.5;
    float beta=0.8;
    float phi = phizero(0);
    float a_n = log_phizero(0);
    vector<int> path_actions(H, 0);
    vector<policy*> sample_policies; //vector with the policies
    vector<float> Vk; // vector containing the mean value of the policy
    double phi_denominator;
    float test_sum = 0;
    int k_max = 20;
    double logphi;
    double max_an = 0;
    double max_bn = 0;


    for (int k = 0; k < k_max; k++) { 
        phi_denominator = 0;
        max_an = 0;
        max_bn = 0;
        //Mk = max(Mk, (int)(1.01 * pow(log(k), 3)));
        //Nk = max(Nk, (int)pow(k, 0.501));
        //cout << Mk << " et " << Nk << endl;
        //beta = 1 / sqrt(k+1);

        for (int n = 0; n < Nk; n++) {
            randfloat = distribution2(generator);
            if (1 - beta > randfloat) {
                //cout << "1-B" << endl;
                sample_policies.push_back(new policy(0));
            }
            else {
                //cout << "B" << endl;
                sample_policies.push_back(new policy(0, t));
            }
            //compute phi(pi,t) in only one time at the creation of the policy
            logphi = log_phiproba(t, sample_policies.back(), 0);
            //sample_policies.back()->phi = phiproba(t, sample_policies.back(), 0);
            sample_policies.back()->a_n = 1/(1+ exp(-abs(beta*logphi - (1-beta)* a_n)));
            //cout << "phi is " <<base->phi << " and log is " <<logphi<< endl;

            test_sum += logphi;
            for (int j = 0; j < Mk; j++) {
                temp_mean += path_value(sample_policies.back(), t, theta);
                //cout << path_value(sample_policies.back(), t, 0) << " - ";
            }
            Vk.push_back(temp_mean / Mk);
            //cout << temp_mean / Mk << " for policy with actions " << sample_policies.back()->actions[0] << endl;

            //phi_denominator += exp(-Vk.back() / Tk) / sample_policies.back()->phi;
            sample_policies.back()->b_n = -Vk.back() / Tk - max(beta * logphi, (1 - beta)* logphi);

            //to normalise
            if (sample_policies.back()->b_n > max_bn) {
                max_an = sample_policies.back()->a_n;
                max_bn = sample_policies.back()->b_n;
            }
            temp_mean = 0;
            //cout << Vk.back() << endl;
        }
        //cout << exp(test_sum) << "sum log" << endl;
        test_sum = 0;
        //cout << t->probaMatrix[0] << endl;
        //cout << "update proba Matrix" << endl;
        //cout << max_an << "  " << max_bn;
       //cout << "phi denominator is " << phi_denominator << endl;
        log_boltzmann(sample_policies, t, k,0, Vk, Tk, max_an,max_bn);
        Vk.clear();
        Tk =Tk / log(k + 0.01); //See red book
        Nk = max(10,min(Nk, (int)(Nk / log(k + 1.01))));
        
        if (k<k_max-1)
            sample_policies.clear();
    }
   
    int arg = arg_max(Vk);
    policy *win = sample_policies[arg];
    int best_a;
    int index_value = 0;
    int index_theta = 0;
    normal_distribution<float> normal(0, sigma_theta );
    for (int h = 0; h < H; h++) {
        index_theta = nearest_theta(atan2(0 - t->position[2], 0 - t->position[0]) + normal(generator));

        if (h == 0) {
            index_value = 0;
        }
        else if (h == 1) {
            index_value = index_theta;
        }
        else if ((h != 0) and (index_value != 0)) {
            index_value = (index_value - 1) * n_theta + index_theta;
        }
        else {
            index_value = index_theta;
        }
           
        //for (int a = 0; a < n_A; a++)
        //    cout << t->probaMatrix[index_value*n_A + a] << " - ";
        //cout << endl;
        best_a = arg_max_ind(t->probaMatrix,index_value,n_A);
        //cout << best_a << "is the best action" << endl;
        path_actions[h] = best_a;
        if (t->childs[best_a] != NULL)
            t = t->childs[best_a];
    }
    *index = index_value;
    //cout << "best is " << endl;
    //win->printpolicy();
    return path_actions;
}

int main(int argc,char** argv)
{
    /* initialize random seed: */
    srand(time(NULL));
    float x1 = x_distrib(generator);
    float x2 = x_distrib(generator);
    //Position and velocity of the target 
    Vector4f x(x1, 0.0, x2, 0);

    //Our position and velocity 
    Vector4f z(-6, 1, -6, 0);

    //Theta history
    vector<float> thetas;


    //Variables initialisation
    float** coeff = new float* [d + 1];
    for (int i = 0; i < d + 1; ++i) {
        coeff[i] = new float[d + 1]; 
    }

    float** coeffcarre = new float* [2 * d + 1];
    for (int i = 0; i < 2 * d + 1; ++i) {
        coeffcarre[i] = new float[2 * d + 1];
    }



    //Distance at which we think the target is
    float initial_dist = 10;
    tree* t = NULL;
    tree* t2 = NULL;

    //create a path with z0 = z
    vector<Vector4f> path = { z };


    //First observation
    int index_theta;
    int action_opti;
    int index_value;


    //Number of simulations
    int N = 30;

    string file_path_a(argv[1]);
    file_path_a.append("/action_file.txt");
    fstream file_a;    //action file
    file_a.open(file_path_a, ios_base::out);

    string file_path_p(argv[1]);
    file_path_p.append("/position_file.txt");
    fstream file_p;     //position file
    file_p.open(file_path_p, ios_base::out);
    file_p << path.back()[0] << "," << path.back()[2] << "," << 0<< "," << 0 << endl;

    string file_path_s(argv[1]);
    file_path_s.append("/start.txt");
   
    fstream file_s;     //start position file 
    file_s.open(file_path_s, ios_base::out | ios_base::in | ios::trunc);
    cout << "hehe";
   
    string file_path_e(argv[1]);
    file_path_e.append("/target.txt");

    fstream file_e;     //start position file
    file_e.open(file_path_e, ios_base::out | ios_base::in | ios::trunc);
    
    file_s << z(0) << "," << z(1) << "," << z(2) << "," << z(3) << endl;
    file_e << x(0) << "," << x(1) << "," << x(2) << "," << x(3) << endl;
    cout << "hoho";

    string file_path_d(argv[1]);
    file_path_d.append("/density_file.txt");
    fstream file_d;     //density file 
    file_d.open(file_path_d, ios_base::out);

    fstream file_v;     //speed file
    file_v.open("valuefunction_file.txt", ios_base::out | ios_base::in);


#if PLA 
        //Do the simulations
        vector<vector<float>> dummy;
        for (int n = 0; n < N; n++) {
            //Create the tree and compute the bellman values and actions associated 
            if (n == 0) {
                routine(&t, path.back(), coeffcarre, coeff, dummy, -1);
            }
            else {
                routine(&t, path.back(), coeffcarre, coeff, dummy, index_value);
            }

            t2 = t;
            for (int h = 0; h < H + 1; h++) {
                thetas.push_back(atan2(x[2] - path.back()[2], x[0] - path.back()[0]) + distribution(generator));
                index_theta = nearest_theta(thetas[n * (H + 1) + h]);
                if (h == 0) {
                    index_value = 0;
                }
                else if (h == 1) {
                    index_value = index_theta;
                }
                else if ((h != 0) and (index_value != 0)) {
                    index_value = (index_value - 1) * n_theta + index_theta;
                }
                else {
                    index_value = index_theta;
                }

                action_opti = t2->value_function[2 * index_value + 1];
                path.push_back(nextPos(path.back(), actions[action_opti], 1));
                cout << "action is " << action_opti << "with value " << t2->value_function[2 * index_value] << " and new position is " << path.back().transpose() << endl;

                //Write to file some results
                dummy = t2->densities;
                file_a << action_opti << endl;
                file_p << path.back()[0] << "," << path.back()[2] << "," << discretizedThetas[index_theta] << "," << thetas.back() << endl;
                for (int i = 0; i < 2 * d + 1; i++) {
                    for (int j = 0; j < 2 * d + 1; j++) {
                        file_d << t2->densities[i][index_value * (2 * d + 1) + j] << ",";
                    }
                    file_d << endl;
                }
                file_v << t2->value_function[2 * index_value] << endl;
                t2 = t2->childs[action_opti];
            }


        }
#else 
    //hello();
    float theta;
    //t = new tree(0, z, H, coeffcarre, coeff);
    file_s << z(0) << "," << z(1) << "," << z(2) << "," << z(3);
    vector<int> path_action(H,0) ;
    int index = -1;
    vector<vector<float>> dens;
    for (int n = 0; n < N; n++) {
        theta = atan2(x[2] - z[2], x[0] - z[0]);
        auto start = chrono::steady_clock::now();
        if (index == -1)
            t = new tree(0, z, H, coeffcarre, coeff);
        else
            t = new tree(0, z, dens, index, H, coeffcarre, coeff);
        auto end = chrono::steady_clock::now();
        chrono::duration<double> elapsed_seconds = end - start;
        cout << "elapsed time for tree creation: " << elapsed_seconds.count() << "s\n";

        start = chrono::steady_clock::now();
        path_action = ASA(t, 30, 10, theta,&index);
        for (int a : path_action) {
            cout << a << " - ";
        }
        //cout << endl;
        end = chrono::steady_clock::now();
        elapsed_seconds = end - start;
        cout << "elapsed time for ASA: " << elapsed_seconds.count() << "s\n";

        for (int h = 0; h < H; h++) {
            //cout << "action " << path_action[h] << endl;
            z = nextPos(z, actions(path_action[h]), 1);
            file_a << path_action[h] << endl;
            file_p << z(0) << ", " << z(2) << endl;
            t = t->childs[path_action[h]];
        }
        dens = t->densities;
        
    }

#endif
    

    //close the stream
    file_a.close();
    file_p.close();
    file_s.close();
    file_d.close();
    file_v.close();


    //Free the allocated spaces.

    for (int i = 0; i < d + 1; ++i)
        delete[] coeff[i];
    delete[] coeff;

    for (int i = 0; i < 2 * d + 1; ++i)
        delete[] coeffcarre[i];
    delete[] coeffcarre;

    freeTree(t2);
    freeTree(t);



}

// Exécuter le programme : Ctrl+F5 ou menu Déboguer > Exécuter sans débogage
// Déboguer le programme : F5 ou menu Déboguer > Démarrer le débogage

// Astuces pour bien démarrer : 
//   1. Utilisez la fenêtre Explorateur de solutions pour ajouter des fichiers et les gérer.
//   2. Utilisez la fenêtre Team Explorer pour vous connecter au contrôle de code source.
//   3. Utilisez la fenêtre Sortie pour voir la sortie de la génération et d'autres messages.
//   4. Utilisez la fenêtre Liste d'erreurs pour voir les erreurs.
//   5. Accédez à Projet > Ajouter un nouvel élément pour créer des fichiers de code, ou à Projet > Ajouter un élément existant pour ajouter des fichiers de code existants au projet.
//   6. Pour rouvrir ce projet plus tard, accédez à Fichier > Ouvrir > Projet et sélectionnez le fichier .sln.

