#pragma once
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
#define BD_ID	1	

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

void save_vector_to_file(const vector<double>& vec, const string& filename) {
    ofstream file(filename);  // Open file for writing

    if (!file) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    file << fixed << setprecision(4);  // Ensure 4 decimal places

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
    file << std::fixed << std::setprecision(6); 
    

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


uniform_real_distribution<> uniform_dist(0.0, 1.0);



void append_adj_matrix(std::ofstream& file, const std::vector<std::vector<int>>& adj, int curr_iter) {
    int n = adj.size();
    file << std::fixed << std::setprecision(6);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            file << adj[i][j] / static_cast<double>(curr_iter);
            if (j < n - 1) file << " ";
        }
        file << "\n";
    }
    file << "\n"; // Blank line separating matrix snapshots
}