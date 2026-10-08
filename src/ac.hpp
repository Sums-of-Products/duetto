#pragma once
#include "newmoves.hpp"
#include "basics.hpp"
#include "DDAG.hpp"
#include "Wsets.hpp"
#include "BD.hpp"
#include "Scorer.hpp"
#include "Gibby.hpp"

//uniform_real_distribution<> uniform_dist(0.0, 1.0);

void ac(DDAG& dag, Scorer_ptr& scop, int max_parents, int n_nodes) {
    vector<int> v(n_nodes);
    iota(v.begin(), v.end(), 0); 


    uniform_int_distribution<> rand_node(0, n_nodes - 1);
    int i = rand_node(gen), j = rand_node(gen);
    while (i == j) {
        j = rand_node(gen);
    }


    vector<int> pj = dag.P[j], pi = dag.P[i];


    vector<int> dj0_vec = dag.descendants(j);
    vector<int> di0_vec = dag.descendants(i);

    set<int> dj0(dj0_vec.begin(), dj0_vec.end());
    set<int> di0(di0_vec.begin(), di0_vec.end());


    string move;
    if (di0.count(j)) {
        move = "deanc";
    } else if (dj0.count(i)) {
        move = "revanc";
    } else {
        move = "anc";
    }

    dag.rep({}, j);
    dag.rep({}, i);

    vector<int> dj_vec = dag.descendants(j);
    vector<int> di_vec = dag.descendants(i);

    set<int> dj(dj_vec.begin(), dj_vec.end());
    set<int> di(di_vec.begin(), di_vec.end());

    std::unordered_set<int> d_set;
    d_set.insert(dj.begin(), dj.end());
    d_set.insert(di.begin(), di.end());

    std::vector<int> d_(d_set.begin(), d_set.end());

    
    vector<int> lpi; vector<int> lpj; double lsli; double lslj; vector<vector<int>> ssi, ssj; vector<int> xi; vector<int> xj;
    if (move == "revanc") {
        lpi = sorted_difference(v, d_);
        lpj = sorted_difference(v, dj_vec);
        lsli = scop -> sum(i, { }, lpi);
        xi = scop->cum( i, { }, lpi, lsli + log(urnd()) );
        lslj = scop -> sum(j, di_vec, lpj);
        xj = scop->cum( j, di_vec, lpj, lslj + log(urnd()) );

    } else if (move == "deanc") {
        vector<int> lp= sorted_difference(v, d_);
        lsli = scop -> sum(i, { }, lp);
        xi = scop->cum( i, { }, lp, lsli + log(urnd()) );
        lslj = scop -> sum(j, { }, lp);
        xj = scop->cum( j, { }, lp, lslj + log(urnd()) );

    } else {
        lpi = sorted_difference(v, d_);
        lpj = sorted_difference(v, dj_vec);
        lsli = scop -> sum(i, { }, lpi);
        xi = scop->cum( i, { }, lpi, lsli + log(urnd()) );
        lslj = scop -> sum(j, di_vec, lpj);
        xj = scop->cum( j, di_vec, lpj, lslj + log(urnd()) );
    }

    vector<int> pj1 = xj, pi1 = xi;

    // dag.rep({}, j);
    // dag.rep({}, i);

    vector<vector<int>> ssi1, ssj1;
    vector<int> lpj1, lpi1;
    double lsli1, lslj1;
    if (move == "revanc") {
        lpj1 = sorted_difference(v, d_);
        lpi1 = sorted_difference(v, di_vec);
        lsli1 = scop -> sum(i, {dj_vec}, lpi1);
        lslj1 = scop -> sum(j, { }, lpj1);

    } else if (move == "deanc") {
        lpj1 = sorted_difference(v, dj_vec);
        lpi1 = sorted_difference(v, d_);
        lsli1 = scop -> sum(i, { }, lpi1);
        lslj1 = scop -> sum(j, {di_vec}, lpj1);

    } else {
        vector<int> lp1= sorted_difference(v, d_);
        lsli1 = scop -> sum(i, { }, lp1);
        lslj1 = scop -> sum(j, { }, lp1);


    }


    dag.rep(pj1, j);
    dag.rep(pi1, i);

    double logr = lslj + lsli - lslj1 - lsli1;

    if (log(uniform_dist(gen)) > min(0.0, logr)) {
        dag.rep(pj, j);
        dag.rep(pi, i);
    }
}