#pragma once
#include "newmoves.hpp"
#include "basics.hpp"
#include "DDAG.hpp"
#include "Wsets.hpp"
#include "BD.hpp"
#include "Scorer.hpp"



#include <vector>
#include <algorithm> // for std::shuffle
#include <random>

// global random engine, or pass it in as an argument
// std::random_device rd;
// std::mt19937 gen(rd());

std::vector<int> random_sample(const std::vector<int>& v, int k) {
    std::vector<int> sample = v; // make a copy
    std::shuffle(sample.begin(), sample.end(), gen);
    sample.resize(k); // keep only first k elements
    return sample;
}

void usr(DDAG& dag, Scorer_ptr& scop, int q, int max_parents, int n_nodes,  int mode = 0) {
    std::vector<int> v(n_nodes);
    std::iota(v.begin(), v.end(), 0);

    std::vector<int> w = random_sample(v, std::min(q, n_nodes)); // replace with your own random_sample implementation

    std::vector<std::pair<int, std::vector<int>>> p0;
    for (int i : w) {
        p0.emplace_back(i, dag.P[i]);
    }

    for (int i : w) {
        dag.rep({}, i);
    }

    double Z_forward = 0.0;
    std::vector<std::pair<int, std::vector<int>>> p1;

//vint pred = w;
//cout << endl << " w: "; print(cout, w);
    for (int i : w) {
        std::vector<int> di = dag.descendants(i);
        std::vector<int> pc = sorted_difference(v, di);
//print(cout, pc);
//pred = sorted_difference(pred, {i});	
//pc = sorted_difference(pc, pred);

        double lsi = scop->sum(i, {}, pc);
        std::vector<int> xi = scop->cum(i, {}, pc, lsi + std::log(urnd()));  // urnd() is uniform [0,1]
//cout << "xi: "; print(cout, xi);
        p1.emplace_back(i, xi);
        dag.rep(xi, i);
        Z_forward += lsi;
    }
    
    if (mode == 0) {
        std::reverse(p0.begin(), p0.end());
    } else if (mode == 1) {
        std::shuffle(p0.begin(), p0.end(), gen);
    }

    for (int i : w) {
        dag.rep({}, i);
    }

    double Z_backward = 0.0;
//vint pree = w;    
    for (auto& [i, p0i] : p0) {
        std::vector<int> di = dag.descendants(i);
        std::vector<int> pc = sorted_difference(v, di);
//pree = sorted_difference(pree, {i});
//pc = sorted_difference(pc, pree);
	
        double lsi = scop->sum(i, {}, pc);
        dag.rep(p0i, i);
        Z_backward += lsi;
    }

    double alpha = Z_forward - Z_backward;

//cout << "alpha = " << alpha << endl;
    if (std::log(urnd()) < std::min(0.0, alpha)) {
//cout << "Accept!" << endl;    
        for (auto& [i, pi] : p1) {
            dag.rep(pi, i);
        }
    } else {
        for (auto& [i, pi] : p0) {
            dag.rep(pi, i);
        }
    }
}
