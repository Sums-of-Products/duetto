#pragma once
#include "newmoves.hpp"
#include "basics.hpp"
#include "DDAG.hpp"
#include "Wsets.hpp"
#include "BD.hpp"
#include "Scorer.hpp"


//uniform_real_distribution<> uniform_dist(0.0, 1.0);

void ec(DDAG& dag, Scorer_ptr& scop, int max_parents, int n_nodes) {
    vector<int> v(n_nodes);
    iota(v.begin(), v.end(), 0); 


    uniform_int_distribution<> rand_node(0, n_nodes - 1);
    int i = rand_node(gen), j = rand_node(gen);
    while (i == j) {
        j = rand_node(gen);
    }


    string move;
    if (dag.arc(i, j)) {
        move = "rem";
    } else if (dag.arc(j, i)) {
        move = "rev";
    } else {
        move = "add";
    }

    vector<int> pj = dag.P[j], pi = dag.P[i];


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
    if (move == "rev") {
        lpj = sorted_difference(v, vector<int>(dj.begin(), dj.end()));
        lpi = sorted_difference(v, d_);
        lsli = scop -> sum(i, { }, lpi);
        xi = scop->cum( i, { }, lpi, lsli + log(urnd()) );
        lslj = scop -> sum(j, {i}, lpj);
        xj = scop->cum( j, {i}, lpj, lslj + log(urnd()) );
        // lpj.erase(remove(lpj.begin(), lpj.end(), i), lpj.end());
        // ssi = subsets(lpi, max_parents);
        // ssj = subsets_with_i(lpj, max_parents, i);
    } else if (move == "rem") {
        lpi = sorted_difference(v, vector<int>(di.begin(), di.end()));
        lpi.erase(remove(lpi.begin(), lpi.end(), j), lpi.end());
        lsli = scop -> sum(i, { }, lpi);
        xi = scop->cum( i, { }, lpi, lsli + log(urnd()) );

        // ssi = subsets(lpi, max_parents);
    } else {
        lpj = sorted_difference(v, vector<int>(dj.begin(), dj.end()));
        lpi = sorted_difference(v, d_);
        lsli = scop -> sum(i, { }, lpi);
        xi = scop->cum( i, { }, lpi, lsli + log(urnd()) );
        lslj = scop -> sum(j, {i}, lpj);
        xj = scop->cum( j, {i}, lpj, lslj + log(urnd()) );
        // lpj.erase(remove(lpj.begin(), lpj.end(), i), lpj.end());
        // ssi = subsets(lpi, max_parents);
        // ssj = subsets_with_i(lpj, max_parents, i);
    }
    // vector<double> lsli;
    // for (const auto& subset : ssi) {
    //     lsli.push_back(sm.get_score(i, subset));
    // }

    // vector<double> gumbels(ssi.size());
    // for (double& g : gumbels) g = sample_gumbel();

    // vector<double> noisy_scores(ssi.size());
    // for (size_t k = 0; k < ssi.size(); ++k) {
    //     noisy_scores[k] = lsli[k] + gumbels[k];
    // }



    //int xi = distance(noisy_scores.begin(), max_element(noisy_scores.begin(), noisy_scores.end()));


    if (move == "rem") {
        //vector<int> pi1 = ssi[xi];
        dag.rep(xi, i);
        vector<int> rdj = dag.descendants(j);
        lpj = sorted_difference(v, rdj);
        lpj.erase(remove(lpj.begin(), lpj.end(), i), lpj.end());
        //ssj = subsets(lpj, max_parents);
        lslj = scop -> sum(j, { }, lpj);
        xj = scop->cum( j, { }, lpj, lslj + log(urnd()) );
    }

    // vector<double> lslj;
    // for (const auto& subset : ssj) {
    //     lslj.push_back(sm.get_score(j, subset));
    // }

    // vector<double> gumbels_j(ssj.size());
    // for (double& g : gumbels_j) g = sample_gumbel();
    

    // vector<double> noisy_scores_j(ssj.size());
    // for (size_t k = 0; k < ssj.size(); ++k) {
    //     noisy_scores_j[k] = lslj[k] + gumbels_j[k];
    // }
    

    // int xj = distance(noisy_scores_j.begin(), max_element(noisy_scores_j.begin(), noisy_scores_j.end()));

    vector<int> pj1 = xj, pi1 = xi;

    dag.rep({}, j);
    dag.rep({}, i);

    vector<vector<int>> ssi1, ssj1;
    vector<int> lpj1, lpi1;
    double lsli1, lslj1;
    if (move == "rev") {
        lpj1 = sorted_difference(v, d_);
        lpi1 = sorted_difference(v, vector<int>(di.begin(), di.end()));
        lsli1 = scop -> sum(i, {j}, lpi1);
        lslj1 = scop -> sum(j, { }, lpj1);

        // lpi1.erase(remove(lpi1.begin(), lpi1.end(), j), lpi1.end());
        // ssj1 = subsets(lpj1, max_parents);
        // ssi1 = subsets_with_i(lpi1, max_parents, j);
    } else if (move == "rem") {
        lpj1 = sorted_difference(v, vector<int>(dj.begin(), dj.end()));
        lpi1 = sorted_difference(v, d_);
        lsli1 = scop -> sum(i, { }, lpi1);
        lslj1 = scop -> sum(j, {i}, lpj1);
        //lpj1.erase(remove(lpj1.begin(), lpj1.end(), i), lpj1.end());
        //ssi1 = subsets(lpi1, max_parents);
        //ssj1 = subsets_with_i(lpj1, max_parents, i);
    } else {
        vector<int> lpi1 = sorted_difference(v, vector<int>(di.begin(), di.end()));
        lpi1.erase(remove(lpi1.begin(), lpi1.end(), j), lpi1.end());
        lsli1 = scop -> sum(i, { }, lpi1);
        //ssi1 = subsets(lpi1, max_parents);
    }

    // vector<double> lsli1;
    // for (const auto& subset : ssi1) {
    //     lsli1.push_back(sm.get_score(i, subset));
    // }

    if (move == "add") {
        dag.rep(pi, i);
        vector<int> rdj = dag.descendants(j);
        vector<int> lpj1 = sorted_difference(v, rdj);
        lpj1.erase(remove(lpj1.begin(), lpj1.end(), i), lpj1.end());
        lslj1 = scop -> sum(j, { }, lpj1);
        //ssj1 = subsets(lpj1, max_parents);
    }

    // vector<double> lslj1;
    // for (const auto& subset : ssj1) {
    //     lslj1.push_back(sm.get_score(j, subset));
    // }

    dag.rep(pj1, j);
    dag.rep(pi1, i);

    //double logr = logsumexp(lslj) + logsumexp(lsli) - logsumexp(lslj1) - logsumexp(lsli1);
    double logr = lslj + lsli - lslj1 - lsli1;

    if (log(uniform_dist(gen)) > min(0.0, logr)) {
        dag.rep(pj, j);
        dag.rep(pi, i);
    }
}