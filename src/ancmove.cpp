#include <iostream>
#include <vector>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <random>
#include <numeric>  // For accumulate()
#include "ScoreManager.hpp"
#include "DDAG.hpp"

using namespace std;




vector<vector<int>> subsets(const vector<int>& vec, int max_parents) {
    vector<vector<int>> t;
    int l = vec.size();

    for (int i = 0; i < min(l + 1, max_parents + 1); i++) {
        if (i == 0) {
            t.push_back({});
        } else {
            vector<int> p(vec.begin(), vec.begin() + i);
            int idx = i - 1;
            int fin = l - 1;
            t.push_back(p);

            vector<int> v(i);
            for (int j = 0; j < i; j++) {
                v[j] = j;
            }

            while (!p.empty() && p[0] != vec[l - i]) {
                bool check = false;

                while (p[idx] == vec[fin]) {
                    check = true;
                    idx--;
                    fin--;
                    if (idx < 0) {
                        return t;
                    }
                }

                if (check) {
                    v[idx]++;
                    p[idx] = vec[v[idx]];

                    for (int j = idx + 1; j < i; j++) {
                        v[j] = v[j - 1] + 1;
                        p[j] = vec[v[j]];
                    }

                    idx = i - 1;
                    fin = l - 1;
                } else {
                    v[idx]++;
                    p[idx] = vec[v[idx]];
                }

                t.push_back(p);
            }
        }
    }

    return t;
}


vector<vector<int>> subsets_with_i(const vector<int>& vec, int max_parents, int i) {
    vector<vector<int>> base_subsets = subsets(vec, max_parents - 1);
    vector<vector<int>> new_subsets;

    for (const auto& subset : base_subsets) {
        vector<int> new_subset;
        bool inserted = false;

        for (int num : subset) {
            if (!inserted && num > i) {
                new_subset.push_back(i);  
                inserted = true;
            }
            new_subset.push_back(num);
        }

        if (!inserted) {
            new_subset.push_back(i);  
        }

        new_subsets.push_back(new_subset);
    }

    return new_subsets;
}

vector<int> sorted_difference(const vector<int>& v1, const vector<int>& v2) {
    unordered_set<int> set_v2(v2.begin(), v2.end());  
    vector<int> result;

    for (int num : v1) {
        if (set_v2.find(num) == set_v2.end()) {  
            result.push_back(num);
        }
    }

    return result;
}

void add_post(DDAG& dag, int n, vector<vector<int>>& adj) {
    for (int j = 0; j < n; ++j) {
        vector<int> p = dag.P[j];  
        for (int i : p) {
            adj[i][j] += 1;  
        }
    }
}

double logsumexp(const vector<double>& vec) {
    if (vec.empty()) return -INFINITY; 

    double max_val = *max_element(vec.begin(), vec.end());  
    double sum_exp = 0.0;

    for (double v : vec) {
        sum_exp += exp(v - max_val);  
    }

    return max_val + log(sum_exp);
}





random_device rd;
mt19937 gen(rd());
uniform_real_distribution<> uniform_dist(0.0, 1.0);


double sample_gumbel() {
    double u = uniform_dist(gen);
    return -log(-log(u));  
}

vector<vector<int>> sub_anc(vector<int> di, int max_parent, int j, const vector<int>& dj, int n) {
    vector<vector<int>> l;
    vector<int> Di = di;  
    vector<int> all(n);
    iota(all.begin(), all.end(), 0);  

    for (int d : dj) {
        all.erase(remove(all.begin(), all.end(), d), all.end());
    }

    for (size_t _ = 0; _ < di.size(); ++_) {
        int a = Di.front();
        Di.erase(Di.begin());

        all.erase(remove(all.begin(), all.end(), a), all.end());

        if (a != j && find(dj.begin(), dj.end(), a) == dj.end()) {
            vector<vector<int>> S = subsets_with_i(all, max_parent, a);
            l.insert(l.end(), S.begin(), S.end());
        }
    }

    return l;
}

vector<vector<int>> sub_deanc(const vector<int>& di, int max_parent, int j, const vector<int>& dj, int n) {
    vector<int> range(n);
    iota(range.begin(), range.end(), 0);

    vector<int> all = sorted_difference(range, di);

    for (int d : dj) {
        all.erase(remove(all.begin(), all.end(), d), all.end());
    }

    vector<vector<int>> S = subsets(all, max_parent);
    return S;
}



void rsampleb(DDAG& dag, ScoreManager& sm, int max_parents, int n_nodes) {
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

    vector<int> d_(dj.begin(), dj.end());
    d_.insert(d_.end(), di.begin(), di.end());

    

    vector<vector<int>> ssi, ssj;
    if (move == "rev") {
        vector<int> lpj = sorted_difference(v, vector<int>(dj.begin(), dj.end()));
        vector<int> lpi = sorted_difference(v, d_);
        lpj.erase(remove(lpj.begin(), lpj.end(), i), lpj.end());
        ssi = subsets(lpi, max_parents);
        ssj = subsets_with_i(lpj, max_parents, i);
    } else if (move == "rem") {
        vector<int> lpi = sorted_difference(v, vector<int>(di.begin(), di.end()));
        lpi.erase(remove(lpi.begin(), lpi.end(), j), lpi.end());
        ssi = subsets(lpi, max_parents);
    } else {
        vector<int> lpj = sorted_difference(v, vector<int>(dj.begin(), dj.end()));
        vector<int> lpi = sorted_difference(v, d_);
        lpj.erase(remove(lpj.begin(), lpj.end(), i), lpj.end());
        ssi = subsets(lpi, max_parents);
        ssj = subsets_with_i(lpj, max_parents, i);
    }

    vector<double> lsli;
    for (const auto& subset : ssi) {
        lsli.push_back(sm.get_score(i, subset));
    }

    vector<double> gumbels(ssi.size());
    for (double& g : gumbels) g = sample_gumbel();

    vector<double> noisy_scores(ssi.size());
    for (size_t k = 0; k < ssi.size(); ++k) {
        noisy_scores[k] = lsli[k] + gumbels[k];
    }


    int xi = distance(noisy_scores.begin(), max_element(noisy_scores.begin(), noisy_scores.end()));


    if (move == "rem") {
        vector<int> pi1 = ssi[xi];
        dag.rep(pi1, i);
        vector<int> rdj = dag.descendants(j);
        vector<int> lpj = sorted_difference(v, rdj);
        lpj.erase(remove(lpj.begin(), lpj.end(), i), lpj.end());
        ssj = subsets(lpj, max_parents);
    }

    vector<double> lslj;
    for (const auto& subset : ssj) {
        lslj.push_back(sm.get_score(j, subset));
    }

    vector<double> gumbels_j(ssj.size());
    for (double& g : gumbels_j) g = sample_gumbel();
    

    vector<double> noisy_scores_j(ssj.size());
    for (size_t k = 0; k < ssj.size(); ++k) {
        noisy_scores_j[k] = lslj[k] + gumbels_j[k];
    }
    

    int xj = distance(noisy_scores_j.begin(), max_element(noisy_scores_j.begin(), noisy_scores_j.end()));

    vector<int> pj1 = ssj[xj], pi1 = ssi[xi];

    dag.rep({}, j);
    dag.rep({}, i);

    vector<vector<int>> ssi1, ssj1;
    if (move == "rev") {
        vector<int> lpj1 = sorted_difference(v, d_);
        vector<int> lpi1 = sorted_difference(v, vector<int>(di.begin(), di.end()));
        lpi1.erase(remove(lpi1.begin(), lpi1.end(), j), lpi1.end());
        ssj1 = subsets(lpj1, max_parents);
        ssi1 = subsets_with_i(lpi1, max_parents, j);
    } else if (move == "rem") {
        vector<int> lpj1 = sorted_difference(v, vector<int>(dj.begin(), dj.end()));
        vector<int> lpi1 = sorted_difference(v, d_);
        lpj1.erase(remove(lpj1.begin(), lpj1.end(), i), lpj1.end());
        ssi1 = subsets(lpi1, max_parents);
        ssj1 = subsets_with_i(lpj1, max_parents, i);
    } else {
        vector<int> lpi1 = sorted_difference(v, vector<int>(di.begin(), di.end()));
        lpi1.erase(remove(lpi1.begin(), lpi1.end(), j), lpi1.end());
        ssi1 = subsets(lpi1, max_parents);
    }

    vector<double> lsli1;
    for (const auto& subset : ssi1) {
        lsli1.push_back(sm.get_score(i, subset));
    }

    if (move == "add") {
        dag.rep(pi, i);
        vector<int> rdj = dag.descendants(j);
        vector<int> lpj1 = sorted_difference(v, rdj);
        lpj1.erase(remove(lpj1.begin(), lpj1.end(), i), lpj1.end());
        ssj1 = subsets(lpj1, max_parents);
    }

    vector<double> lslj1;
    for (const auto& subset : ssj1) {
        lslj1.push_back(sm.get_score(j, subset));
    }

    dag.rep(pj1, j);
    dag.rep(pi1, i);

    double logr = logsumexp(lslj) + logsumexp(lsli) - logsumexp(lslj1) - logsumexp(lsli1);

    if (log(uniform_dist(gen)) > min(0.0, logr)) {
        dag.rep(pj, j);
        dag.rep(pi, i);
    }
}

void ancmove(DDAG& dag, ScoreManager& sm, int max_parents, int n_nodes) {
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

    vector<int> d_(dj.begin(), dj.end());
    d_.insert(d_.end(), di.begin(), di.end());


    vector<vector<int>> ssi, ssj;
    if (move == "revanc") {
        vector<int> lpi = sorted_difference(v, d_);
        ssi = subsets(lpi, max_parents);
        ssj = sub_anc(di_vec, max_parents, j, dj_vec, n_nodes);
    } else if (move == "deanc") {
        ssi = sub_deanc(dj_vec, max_parents, i, di_vec, n_nodes);
        ssj = sub_deanc(di_vec, max_parents, j, dj_vec, n_nodes);
    } else {
        vector<int> lpi = sorted_difference(v, d_);
        ssi = subsets(lpi, max_parents);
        ssj = sub_anc(di_vec, max_parents, j, dj_vec, n_nodes);
    }

    vector<double> lsli;
    for (const auto& subset : ssi) {
        lsli.push_back(sm.get_score(i, subset));
    }

    vector<double> gumbels_i(ssi.size()), noisy_scores_i(ssi.size());
    for (size_t k = 0; k < ssi.size(); ++k) {
        gumbels_i[k] = sample_gumbel();
        noisy_scores_i[k] = lsli[k] + gumbels_i[k];
    }
    int xi = distance(noisy_scores_i.begin(), max_element(noisy_scores_i.begin(), noisy_scores_i.end()));

    vector<double> lslj;
    for (const auto& subset : ssj) {
        lslj.push_back(sm.get_score(j, subset));
    }

    vector<double> gumbels_j(ssj.size()), noisy_scores_j(ssj.size());
    for (size_t k = 0; k < ssj.size(); ++k) {
        gumbels_j[k] = sample_gumbel();
        noisy_scores_j[k] = lslj[k] + gumbels_j[k];
    }
    int xj = distance(noisy_scores_j.begin(), max_element(noisy_scores_j.begin(), noisy_scores_j.end()));

    vector<int> pj1 = ssj[xj];
    vector<int> pi1 = ssi[xi];

    vector<vector<int>> ssi1, ssj1;
    if (move == "revanc") {
        vector<int> lpj1 = sorted_difference(v, d_);
        ssi1 = sub_anc(dj_vec, max_parents, i, di_vec, n_nodes);
        ssj1 = subsets(lpj1, max_parents);
    } else if (move == "deanc") {
        vector<int> lpi1 = sorted_difference(v, d_);
        ssi1 = subsets(lpi1, max_parents);
        ssj1 = sub_anc(di_vec, max_parents, j, dj_vec, n_nodes);
    } else {
        ssi1 = sub_deanc(dj_vec, max_parents, i, di_vec, n_nodes);
        ssj1 = sub_deanc(di_vec, max_parents, j, dj_vec, n_nodes);
    }

    vector<double> lsli1;
    for (const auto& subset : ssi1) {
        lsli1.push_back(sm.get_score(i, subset));
    }

    vector<double> lslj1;
    for (const auto& subset : ssj1) {
        lslj1.push_back(sm.get_score(j, subset));
    }

    dag.rep(pj1, j);
    dag.rep(pi1, i);

    double logr = logsumexp(lslj) + logsumexp(lsli) - logsumexp(lslj1) - logsumexp(lsli1);

    if (log(uniform_dist(gen)) > min(0.0, logr)) {
        dag.rep(pj, j);
        dag.rep(pi, i);
    }
}


double score(DDAG& dag, ScoreManager& sm, int n) {
    double s = 0.0; 

    for (int i = 0; i < n; i++) {
        s += sm.get_score(i, dag.P[i]);  
    }

    return s;
}

void save_vector_to_file(const vector<double>& vec, const string& filename) {
    ofstream file(filename);  // Open file for writing

    if (!file) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    for (double num : vec) {
        file << num << endl;  
    }

    file.close();
}

void save_matrix_to_file(const vector<vector<int>>& adj, int iterations, const string& filename) {
    ofstream file(filename); 

    if (!file) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    int n = adj.size();
    

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            file << (adj[i][j] / static_cast<double>(iterations));  
            if (j < n - 1) file << " ";  
        }
        file << endl; 
    }

    file.close();
}




void print_subsets(const vector<vector<int>>& subsets, const string& label) {
    cout << label << ":\n";
    for (const auto& subset : subsets) {
        cout << "  [";
        for (size_t i = 0; i < subset.size(); ++i) {
            cout << subset[i];
            if (i < subset.size() - 1) cout << ", ";
        }
        cout << "]\n";
    }
}


int main() {
    try {


        int max_parents = 3;
        int n_nodes = 56;
        int iterations = 1000;

        DDAG dag(n_nodes, false);

        ScoreManager sm("/Users/dnikzad/Documents/code/datasets/hail1k3.jkl", max_parents, n_nodes);


        vector<vector<int>> adj(n_nodes, vector<int>(n_nodes, 0));

        vector<double> values;


        for (int iter = 0; iter < iterations; ++iter) {
            for (int _ = 0; _ < 100; ++_) {
                ancmove(dag, sm, max_parents, n_nodes);
            }
            values.push_back(score(dag, sm, n_nodes));
            add_post(dag, n_nodes, adj);
        }
        
        save_vector_to_file(values, "output.txt");
        save_matrix_to_file(adj, iterations, "matrix.txt");

        cout << "Sample complete" << endl;
    } catch (const exception& e) {
        cerr << "Exception caught: " << e.what() << endl;
    }

    return 0;
}