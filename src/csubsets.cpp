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

// Existing subsets function (same as before)
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

// Function to insert 'i' while maintaining order
vector<vector<int>> subsets_with_i(const vector<int>& vec, int max_parents, int i) {
    vector<vector<int>> base_subsets = subsets(vec, max_parents - 1);
    vector<vector<int>> new_subsets;

    for (const auto& subset : base_subsets) {
        vector<int> new_subset;
        bool inserted = false;

        for (int num : subset) {
            if (!inserted && num > i) {
                new_subset.push_back(i);  // Insert 'i' before the first larger number
                inserted = true;
            }
            new_subset.push_back(num);
        }

        if (!inserted) {
            new_subset.push_back(i);  // If 'i' wasn't inserted, add it at the end
        }

        new_subsets.push_back(new_subset);
    }

    return new_subsets;
}

vector<int> sorted_difference(const vector<int>& v1, const vector<int>& v2) {
    unordered_set<int> set_v2(v2.begin(), v2.end());  // Convert v2 into a set for fast lookup
    vector<int> result;

    for (int num : v1) {
        if (set_v2.find(num) == set_v2.end()) {  // If num is not in v2, add to result
            result.push_back(num);
        }
    }

    return result;
}

void add_post(DDAG& dag, int n, vector<vector<int>>& adj) {
    for (int j = 0; j < n; ++j) {
        vector<int> p = dag.P[j];  // Get parents of node j
        for (int i : p) {
            adj[i][j] += 1;  // Increment adjacency matrix
        }
    }
}

double logsumexp(const vector<double>& vec) {
    if (vec.empty()) return -INFINITY;  // Handle empty input safely

    double max_val = *max_element(vec.begin(), vec.end());  // Find max for numerical stability
    double sum_exp = 0.0;

    for (double v : vec) {
        sum_exp += exp(v - max_val);  // Shift to prevent overflow
    }

    return max_val + log(sum_exp);
}




// Random number generators
random_device rd;
mt19937 gen(rd());
uniform_real_distribution<> uniform_dist(0.0, 1.0);

// Function to generate Gumbel-distributed numbers
double sample_gumbel() {
    double u = uniform_dist(gen);
    return -log(-log(u));  // Gumbel distribution sampling
}

void rsampleb(DDAG& dag, ScoreManager& sm, int max_parents, int n_nodes) {
    vector<int> v(n_nodes);
    iota(v.begin(), v.end(), 0); // Fill v with 0 to n_nodes-1

    // Select two different random nodes i and j
    uniform_int_distribution<> rand_node(0, n_nodes - 1);
    int i = rand_node(gen), j = rand_node(gen);
    while (i == j) {
        j = rand_node(gen);
    }

    // Determine the move type: "rem" (remove), "rev" (reverse), or "add" (add)
    string move;
    if (dag.arc(i, j)) {
        move = "rem";
    } else if (dag.arc(j, i)) {
        move = "rev";
    } else {
        move = "add";
    }

    vector<int> pj = dag.P[j], pi = dag.P[i];

    // Reset parent sets
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

    // Find index of the max noisy score
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
    
    // Add Gumbel noise to scores
    vector<double> noisy_scores_j(ssj.size());
    for (size_t k = 0; k < ssj.size(); ++k) {
        noisy_scores_j[k] = lslj[k] + gumbels_j[k];
    }
    
    // Find index of the max noisy score
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

double score(DDAG& dag, ScoreManager& sm, int n) {
    double s = 0.0;  // Sum of scores

    for (int i = 0; i < n; i++) {
        s += sm.get_score(i, dag.P[i]);  // Accumulate score
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
        file << num << endl;  // Write each number on a new line
    }

    file.close();  // Close the file
}

void save_matrix_to_file(const vector<vector<int>>& adj, int iterations, const string& filename) {
    ofstream file(filename);  // Open file for writing

    if (!file) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    int n = adj.size();
    
    // Normalize and save as space-separated values (NumPy-compatible)
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            file << (adj[i][j] / static_cast<double>(iterations));  // Normalize
            if (j < n - 1) file << " ";  // Space separator (NumPy-friendly)
        }
        file << endl;  // New line after each row
    }

    file.close();
}




int main() {
    try {
        // Initialize DDAG with 8 nodes and false (assuming false means unweighted or some property)

        int max_parents = 3;
        int n_nodes = 56;
        int iterations = 1000;

        DDAG dag(n_nodes, false);

        // Initialize ScoreManager with dataset file, 4, and 17 as parameters
        ScoreManager sm("/Users/dnikzad/Documents/code/datasets/hail1k3.jkl", max_parents, n_nodes);

        vector<vector<int>> adj(n_nodes, vector<int>(n_nodes, 0));

        vector<double> values;

        // Run rsampleb 1000 times
        for (int iter = 0; iter < iterations; ++iter) {
            for (int _ = 0; _ < 100; ++_) {
                rsampleb(dag, sm, max_parents, n_nodes);
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