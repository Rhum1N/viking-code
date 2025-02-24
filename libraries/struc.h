#pragma once
#include <iostream>
#include <cstdlib>
#include <Eigen/Core>
#include <vector>

using namespace std;
using namespace Eigen;

#define PI 3.14159265
#define H 2 //The horizon
const int n_A = 5;
const int n_theta = 10; //Number of thetas...


struct tree {
    int root;
    Vector4f position; //position Z of the node 
    tree* childs[n_A];
    tree* previous; //pointer to the father of the node
    vector<vector<float>> densities; // matrix with all the densities
    vector<double> value_function; // vector containing the value function result and the action chosen 
    vector<double> norms; // vector containing the denominator of the density
    vector<float> probaMatrix; // vector of size |X|*n_A with the probabilities of chosing each action sum(probaMatrix) = 1  

    tree(int a, Vector4f z) {
        root = a;
        position = z;

        for (int i = 0; i < n_A; i++) {
            childs[i] = NULL;
        }
    }
    
};


struct policy {
    vector<int> actions;
    policy* childs[n_A];

    policy(int h) {
        actions = vector<int>(pow(n_theta,h), 0);
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

    policy(int h, tree* t) {

    }
};