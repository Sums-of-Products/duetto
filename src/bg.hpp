#pragma once
#include "newmoves.hpp"
#include "basics.hpp"
#include "DDAG.hpp"
#include "Wsets.hpp"
#include "BD.hpp"
#include "Scorer.hpp"



#include <iostream>
#include <vector>
#include <set>
#include <tuple>
#include <algorithm>
#include <random>
#include <numeric>
#include <cmath>
#include <limits>

using namespace std;

int sample_cumsum(const std::vector<double>& a) {
    // Set up random number generator
    static std::random_device rd;
    //static std::mt19937 gen(rd());
    static std::mt19937 gen(13); // [Mikko]

    std::uniform_real_distribution<> dis(0.0, a.back());

    double r = dis(gen);
    int j = 0;
    double s = a[0];

    while (j + 1 < static_cast<int>(a.size()) && r > s) {
        ++j;
        s = a[j];
    }

    return j;
}


int sample_log(const std::vector<double>& log_probs) {
    // Setup random generator and distribution
    // static std::random_device rd;
    // static std::mt19937 gen(rd());
    // static std::uniform_real_distribution<> dis(0.0, 1.0);

    double logr = std::log(uniform_dist(gen)); // log of uniform random in (0,1)
    int j = 0;
    double s = std::exp(log_probs[0]);

    // Cumulative sum of exp(log_probs)
    while (j + 1 < static_cast<int>(log_probs.size()) && logr > std::log(s)) {
        j++;
        s += std::exp(log_probs[j]);
    }

    return j;
}


// Function to check if a graph is a DAG using DFS
bool is_dag(int n, const vector<pair<int, int>>& edges) {
    vector<vector<int>> adj(n);
    for (auto [u, v] : edges)
        adj[u].push_back(v);

    vector<int> visited(n, 0); // 0 = unvisited, 1 = visiting, 2 = visited

    function<bool(int)> dfs = [&](int u) {
        visited[u] = 1;
        for (int v : adj[u]) {
            if (visited[v] == 1)
                return false; // cycle detected
            if (visited[v] == 0 && !dfs(v))
                return false;
        }
        visited[u] = 2;
        return true;
    };

    for (int i = 0; i < n; ++i)
        if (visited[i] == 0 && !dfs(i))
            return false;

    return true;
}

vector<vector<vector<vector<int>>>> gad(int n) {
    vector<pair<int, int>> all_edges;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (i != j)
                all_edges.emplace_back(i, j);

    set<vector<pair<int, int>>> seen;
    vector<vector<vector<vector<int>>>> dags_info;

    int total_edges = all_edges.size();
    for (int mask = 0; mask < (1 << total_edges); ++mask) {
        vector<pair<int, int>> edges;
        for (int i = 0; i < total_edges; ++i)
            if ((mask >> i) & 1)
                edges.push_back(all_edges[i]);

        if (!is_dag(n, edges))
            continue;

        vector<pair<int, int>> edge_key = edges;
        sort(edge_key.begin(), edge_key.end());
        if (seen.count(edge_key))
            continue;
        seen.insert(edge_key);

        // Create parent list
        vector<vector<int>> parents(n);
        for (auto [u, v] : edges)
            parents[v].push_back(u);

        // Create [parents, non-parents] list for each node
        vector<vector<vector<int>>> dag_summary(n);
        for (int i = 0; i < n; ++i) {
            vector<int> p = parents[i];
            vector<int> non_p;
            for (int j = 0; j < n; ++j)
                if (j != i && find(p.begin(), p.end(), j) == p.end())
                    non_p.push_back(j);
            dag_summary[i] = {p, non_p};
        }

        dags_info.push_back(dag_summary);
    }

    return dags_info;
}

void generate_combinations(const vector<int>& elements, int length, int start,
                           vector<int>& current, vector<vector<int>>& results) {
    if (current.size() == length) {
        results.push_back(current);
        return;
    }
    for (int i = start; i < elements.size(); ++i) {
        current.push_back(elements[i]);
        generate_combinations(elements, length, i + 1, current, results);
        current.pop_back();
    }
}

vector<vector<int>> cand(int max_p,
                         const vector<int>& a,
                         const vector<vector<int>>& b,
                         const vector<vector<int>>& c) {
    // Flatten c into a set
    set<int> c_flat;
    for (const auto& group : c)
        for (int x : group)
            c_flat.insert(x);

    // Filter a to exclude elements in c_flat
    vector<int> valid_elements;
    for (int x : a)
        if (c_flat.find(x) == c_flat.end())
            valid_elements.push_back(x);

    vector<vector<int>> result;

    for (int length = 0; length <= max_p; ++length) {
        vector<vector<int>> combinations;
        vector<int> current;
        generate_combinations(valid_elements, length, 0, current, combinations);

        for (const auto& combo : combinations) {
            // Check that each group in b has at least one element in combo
            bool valid = true;
            for (const auto& group : b) {
                bool found = false;
                for (int x : group) {
                    if (find(combo.begin(), combo.end(), x) != combo.end()) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    valid = false;
                    break;
                }
            }
            if (valid)
                result.push_back(combo);
        }
    }

    return result;
}


// The code below added by Mikko on August 1–5, 2026.

int sample_gumbel(vdou &v){
    // Add Gumbel noise and pick index with maximum
    vdou gumbels(v.size()); generate(gumbels.begin(), gumbels.end(), []() { return -log(-log(urnd())); });
    for (size_t i = 0; i < v.size(); ++i) gumbels[i] += v[i];
    return distance(gumbels.begin(), std::max_element(gumbels.begin(), gumbels.end()));
}

// Returns the list of indices to the full list of parent set score for node i = W[v].
vint sBG_at_v(int n, int v, const vint& W, const vector<u8>& H, const vector<dbitset>& psets, const vector<dbitset>& bD){
	int i = W[v], q = W.size();	// Have to "recompute" these, since they are not given as direct arguments.
	
	dbitset pv = psets[i*n+i];  // Default. Possibly overrided next.
	// Compute pvA, and subtract at the end from pv;
	// First compute Qv under H. Qv = N \ (union of Du for u in W \ Hv).
	// Then loop over nodes *not* in Qv. That is, (union of Du for u in W \ Hv).
	// 1. Compute the union of Du's.
	dbitset cQv(n); // Subset of nodes represented as a bitset.
	u8 bits = H[v] ^ ((1 << q) - 1);; // Complement of H[v].
	int u = 0;
	while (bits){
		if (bits & 1) cQv |= bD[u]; // u is *not* a parent of v in H.
		bits >>= 1; u++;
	}
	// Remove i. This is due to the special definition of psets[i*n+i].
	cQv.reset(i);
			
	dbitset pvA(pv.size());	// Of same size as pv but initialized with 0s.
	// 2. Loop over nodes in cQv.
	vint vcQv = make_vint(cQv);
	for (auto j : vcQv) pvA |= psets[i*n+j];
			
	// Now pvA is in place. Next have to possibly change pv to pvB. Actually, just overwrite pv.				
	if (H[v]){  // If in H node v has parents. Have to compute pvB
		pv.set(); // Of same size as pvA but initialized with 1s.
		u8 bits = H[v];
		int u = 0;
		while (bits){ // Loop over the parent u of v in H.
			// Compute Rvu under H. Rvu = Du \ (union of Du' for u' in W \ Hv) = Dx \ cQv.
		 	if (bits & 1){
		 		int j = W[u];
		 		dbitset Rvu = bD[u]- cQv;
				vint vRvu = make_vint(Rvu);
		 		dbitset pvBu(pvA.size()); // To contain the union, initialized to 0s.
		 		for (auto j : vRvu) pvBu |= psets[i*n+j];
		 		pv &= pvBu; // Take intersections.	
		 	}
		 	bits >>= 1; u++;				
		}	
	}
	// Finally, subtract.
	pv -= pvA;
	// Now pv indicates the target parent sets, of which scores we want to add up.
	return make_vint(pv);
}

// Running time:
// 100 000 steps takes about 200 minutes : Hailfinder10k, no score pruning, q = 3
// 100 000 steps takes about 1/3 minutes : Hailfinder10k,    score pruning, q = 3
// 100 000 steps takes about   7 minutes : Hailfinder10k,    score pruning, q = 4
//
void bg (DDAG& dag, Scorer_ptr& scop, int q){
	// 0. Introduce static variables.
	static int count = 0;				// Counter for function calls. 
	static int n = scop->n;				// Number of nodes. Fixed by the first call.
	static vvdo lsc(n);					// Vector of local scores per node, for quick accessing. 
	static vector<dbitset> psets(n*n);	// psets[i*p+j] is a bitmap of parents sets of i that do not contain j.
	
	if (count == 0){ // NOTE: Done only at the first call of the function.
		// Local scores and psets (as bitmaps).
		for (int i = 0; i < n; i++){
			int m = scop->nscores(i);
			lsc[i].resize(m);
			for (int j = 0; j < n; j++) psets[i*n+j].resize(m); // Note: Initially all bits set to 0.
			for (int t = 0; t < m; t++){
				lsc[i][t] = scop->score(i, t);
				dbitset pit = scop->pset(i, t);	// Do not store all these, just use temporarily.
				for (int j = 0; j < n; j++) if (pit.test(j)) psets[i*n+j].set(t); // The t-th set is included if it contains j.
				psets[i*n+i].set(); // Special: encodes the full collection of parent sets for i.	
			}
		}
	}
	count++;
	
	// 1. Draw the node block W = {W[0],...,W[k-1]} of size q.
	q = min(q, n); // Just in case.
	vint W = rperm(n, q);

	// 2. Orphan the current graph at W.
    for (auto i : W) dag.rep({}, i);

	// 3. Store the set of descendants of node W[v] as D[v] for v = 0,...,q-1.
    vvin D; for (auto i : W) D.push_back(dag.descendants(i));	
	// Actually, we will use the bitset representations of these sets.
	vector<dbitset> bD(q);
	for (int v = 0; v < q; v++) bD[v] = make_dbitset(D[v], n);

	// 4. Score each DAG H on W. This step significantly dominates the time complexity; all the rest is peanuts.
	vector<u8> H(q);
	// d2, d3, d4 given in basics.hpp. Point to one of them by d.
	int l; u8* d;
	switch (q){
	case 2: l =   3; d = &d2[0][0]; break;
	case 3: l =  25; d = &d3[0][0]; break;
	case 4: l = 543; d = &d4[0][0]; break;	
	default: cerr << "Block size " << q << " not implemented. Exit now. " << endl; exit(1); 
	}
	vdou val(l); // Here the scores of the DAGs H are stored. Initialized to 0s.

	for (int h = 0; h < l; h++){
		// Handle the h-th DAG H on W. 
		for (int v = 0; v < q; v++) H[v] = d[h*q+v];
				
		// One node i in W in turn.
		for (int v = 0; v < q; v++){
			int i = W[v];		
			vint vpv = sBG_at_v(n, v, W, H, psets, bD);
			vdou si(vpv.size()); // Dense vector of local scores of node i under H.			
			for (int t = 0; t < vpv.size(); t++) si[t] = lsc[i][vpv[t]];
			
			// Add the log of the score sum to the toal for DAG H.  
			val[h] += logsumexp(si);
		}
	}
	
	// 5. Sample a DAG H.
	int h = sample_gumbel(val);
	for (int v = 0; v < q; v++) H[v] = d[h*q+v];
	
	// 6. Sample new parent sets for i in W using H, and replace the current parents sets.
	for (int v = 0; v < q; v++){
		int i = W[v];
		// Not only copy sparse-to-dense local scores but also the full indices to vpv.
		// Then don't call scop->pset(i, t) but scop->pset(i, vpv[t]).
		vint vpv = sBG_at_v(n, v, W, H, psets, bD);
		vdou si(vpv.size()); // Dense vector of local scores of node i under H.			
		for (int t = 0; t < si.size(); t++) si[t] = lsc[i][vpv[t]];

		int t = sample_gumbel(si);
		vint pi = make_vint( scop->pset(i, vpv[t]) );		
		dag.rep(pi, i);
	}
}

// *** An O(4^q) time monotone algorithm for sampling a DAG on q nodes from a modular distribution (Talvitie et al.). ***
// NOTE: Z uses a sparse encoding: all q nodes assume 2^q parent sets, including set that contain the node itself.
//
#define Treal	B2real
#define vTrea	vector<Treal>
#define vvTre	vector<vTrea>
#define vvvTr	vector<vvTre>
// Samples an index t between a and b with pr proportional to w[t].
u32 sample(vTrea& w, u32 a, u32 b){
	Treal sum; sum = 0;
	for (u32 t = a; t <= b; t++) sum += w[t];
	double u = urnd();
	Treal target; target = u * sum;
	u32 t = a; sum = w[t];
	while (sum < target) sum += w[++t];
	return t;
}

// Samples a subset J of U that intersects T with pr proportional to w[S].
u32 isample(vTrea& w, u32 U, u32 T){
	if (T == 0) return 0;
	Treal sum; sum = 0;
	for (u32 J = U; J > 0; J = (J - 1) & U) if (J & T) sum += w[J];
	Treal target; target = urnd() * sum;
	u32 J = U; sum = w[U];
	while (sum < target){ J = (J - 1) & U; if (J & T){ sum += w[J]; } }
	return J;
}

// Fast zeta transform. Assumes that the size of f is a power of 2.
void fzt(const vTrea& f, vTrea& g){
	u32 N = f.size(); int n = popcount32(N-1);
	g = f; for (int i = 0; i < n; i++) for (u32 S = 0; S < N; S++) if (S & (1 << i)) g[S] += g[S ^ (1 << i)];
}

// 119 lines
vector<u32> sample_dag_old(int q, vvdo& Z){ // The older version.
	static int count = 0;
	vector<u32> H(q);
	// Denote V = {0,...,q-1}. However, for convenience, define the following:
	u32 V = (1 << q);

	static vvTre w(q);		// Size q 2^q 
	static vvvTr what(q); 	// Size q 3^q
	static vvTre f(V);		// Size 3^q
	if (count == 0){ // NOTE: This assumes that q will be the same in each call!!!
		for (auto& wi : w) wi.resize(V);
		for (auto& whati : what){ whati.resize(V); for (u32 T = 0; T < V; T++) whati[T].resize(1 << popcount16(T)); } 
		for (u32 U = 0; U < V; U++) f[U].resize(1 << popcount16(U));
	}	
	count++;

	// 0. Compute w[i][J] from Z[][]. 
	for (int i = 0; i < q; i++) for (u32 J = 0; J < V; J++){ w[i][J].set_log(Z[i][J]); if (J & (1 << i)) w[i][J] = 0; }
	
	// 1. Compute w_hat[i][T][R] for all i in V, all T subset of V \ {i} and all R subset of T.
	// That should take time O(q 3^q) using monotone arithmetic.
	vTrea wir(V); vTrea whatir(V); // Auxilliary vectors.
	for (int i = 0; i < q; i++){
		// Compute the base cases where R is empty or a singleton set.
		// Takes time O(q 3^q).
//		fzt(w[i], whatir); for (u32 T = 0; T < V; T++) what[i][T][0] = whatir[T];

		for (u32 T = 0; T < V; T++){
			what[i][T][0] = w[i][0]; for (u32 R = 1; R < V; R <<= 1) if (T & R) what[i][T][proj16(T, R)] = 0; // Slow but enough.
			for (u32 J = T; J > 0; J = (J - 1) & T) what[i][T][0] += w[i][J]; // Could be made faster using FZT, but not the bottleneck.	
		}

		for (int r = 0; r < q; r++){ // Compute what[i][T][{r}] for all T and j. Use FZT. Time O(q^2 2^q) in total.
			// Filter w[i] based on j to wir.
			u32 R = 1 << r;
			for (u32 J = 0; J < V; J++){ wir[J] = 0; if (J & R) wir[J] = w[i][J]; }
			fzt(wir, whatir);
			for (u32 T = 0; T < V; T++) if (T & R) what[i][T][proj16(T, R)] = whatir[T];
		}
		
		// Then apply the recurrence what[i][T][R] = what[i][T][{r}] + what[i][T\{s}][R\{s}], where s is an arbitrary element of R.
		for (u32 T = 0; T < V; T++){
			for (u32 cR = T; cR > 0; cR = (cR - 1) & T){ // Only nonempty subsets cR of T. See Gosper's hack and related literature.
				u32 R = T ^ cR; // Consider the complement to process subsets before their supersets.
				if (popcount16(R) < 2) continue;
				u32 S = R & -R; // The least significant bit of R, raised to the power of 2.
				what[i][T][proj16(T, R)] = what[i][T][proj16(T, S)] + what[i][T ^ S][proj16(T ^ S, R ^ S)];
			}
			// Did not yet process cR = 0, i.e., R = T. Do it now:
			if (popcount16(T) > 1){ u32 S = T & -T; what[i][T][proj16(T, T)] = what[i][T][proj16(T, S)] + what[i][T ^ S][proj16(T ^ S, T ^ S)]; }
		}

	}
		
	// 2. Compute f(R, U) for all U subset of V and all R subset of U.
	// That should take time O(4^q) using monotone arithmetic.

	// Initialize.
	Treal val1; val1 = 1.0;
	for (u32 U = 0; U < V; U++){
		u32 dU = proj16(U, U);
		for (u32 pR = 1; pR <= dU; pR++) f[U][pR] = 0;
		f[U][0] = 0; 
		f[U][dU] = val1; // Perhaps only this line is really needed. Reconsider.
	}

	// Compute f by dynamic programming. See Section 3.1.1 of Talvitie et al.
	vTrea wprod(V);
	vTrea whats(q);
	
	for (u32 U = 0; U < V; U++){
		for (u32 T = U; T > 0; T = (T - 1) & U){
			u32 R = U ^ T; // Consider the complement to process subsets before their supersets.
			if (R == 0) continue; // Require nonempty R.
			// Compute f[U][R] := ...
			// Loop over all nonempty subsets S of T = U \ R, in increasing order by |S|. NOTE: This is the time-critical inner loop!
			// Construct the mapping whats: whats[s] = what[index(s)][V \ T][R].
			u32 S = T; int s = 0; u32 Vl = (V-1) ^ T; u32 RtoVl = proj16(Vl, R);
			while (S > 0){ whats[s++] = what[lsb16(S)][Vl][RtoVl]; S &= (S - 1); } // NOTE: lsb16(S) is in the full range.  
			u32 dT = proj16(T, T); // Dense encoding of T.
			wprod[0] = val1;
			Treal sum; sum = 0;
			for (u32 S = 1; S <= dT; S++){ // We use directly the dense index; no need to use proj16(.,.).
				wprod[S] = whats[lsb16(S)] * wprod[S & (S - 1)]; // NOTE: lsb16(S) is in the projected range.
				sum += f[T][S] * wprod[S];
			}
			f[U][proj16(U, R)] = sum;
		}
	}

	// 3. Loop over U to peel off the parts of a root-partition.
	// That should take time O(q 2^q) using monotone arithmetic.
	u32 U = V-1; u32 Rk = 0; int k = 0;
	vTrea weight(V);
	while (U){
		// Compute the weight of each subset R of U. Will draw R with probability proportional to the weight. 
		u32 dU = proj16(U, U);		
		for (u32 pR = 1; pR <= dU; pR++) weight[pR] = f[U][pR]; // Could simply use weight as a pointer instead of copying.
		u32 S = U; int s = 0; u32 R1tok = (V-1) ^ U; u32 pRk = proj16(R1tok, Rk);
		while (S > 0){ whats[s++] = what[lsb16(S)][R1tok][pRk]; S &= (S - 1); } // For each s-th element i in U store what[i][][].
		for (u32 pR = 1; pR <= dU; pR++){
			u32 S = pR; int s = 0;
			while (S){ if (S & 1){ weight[pR] = weight[pR] * whats[s]; } S >>= 1; s++; } 
		}		
		weight[0] = 0;

		// Draw R, i.e., R(k+1).
		u32 pR = sample(weight, 1, dU);
		// As pR is relative to U, we have to decompress it to be relative to V-1.
		u32 R = U; u32 r = dU;
		while (r > pR){ R = (R - 1) & U; r--; }
		
		// Don't update Rk yet, as we need the current value now.
		// Draw the parent set for each node in R. The parent set must intersect Rk.
		for (int i = 0; i < q; i++){
			// If node i is in R. Consider w[i][Hi] with Hi intersecting Rk to pick H[i] = Hi;
			if (R & (1 << i)) H[i] = isample(w[i], R1tok, Rk);
		}
		// Now update Rk to R and U to U \ R.
		Rk = R; k++; U ^= R; 
	}
	return H;
}

// 145 lines
vector<u32> sample_dag(int q, vvdo& Z){
	static int count = 0;
	vector<u32> H(q);
	// Denote V = {0,...,q-1}. However, for convenience, define the following:
	u32 V = (1 << q);

	static vvTre w(q);		// Size q 2^q 
	//static vvvTr what(q); 	// Size q 3^q
	static vvTre f(V);		// Size 3^q
	static vvvTr wh(q);		// Size q^2 2^q
	static vvTre why(V);	// Size q 2^q
	static vTrea wix(V);	// Size 2^q
	static vTrea wprod(V);
	static vTrea whats(q);
	static vTrea weight(V);

	static vTrea wir(V); 
	static vTrea whatir(V); // Auxilliary vectors.

	if (count == 0){ // NOTE: This assumes that q will be the same in each call!!!
		for (auto& wi : w) wi.resize(V);
		for (u32 U = 0; U < V; U++) f[U].resize(1 << popcount16(U));
		for (auto& whi : wh){ whi.resize(q + 1); for (int x = 0; x < q + 1; x++) whi[x].resize(V); }
		for (auto& whyU : why) whyU.resize(q);
	}	
	count++;

	// 0. Compute w[i][J] from Z[][]. 
	for (int i = 0; i < q; i++) for (u32 J = 0; J < V; J++){ w[i][J].set_log(Z[i][J]); if (J & (1 << i)) w[i][J] = 0; }

for (int rep = 0; rep < 1; rep++){
			
	// 1. Compute w_hat[i][T][R] for all i in V, all T subset of V \ {i} and all R subset of T.
	// That should take time O(q 3^q) using monotone arithmetic.
	for (int i = 0; i < q; i++){
		// Compute the base cases where R is empty or a singleton set.
		// Takes time O(q^2 2^q).
		// Here wix[] is a temporary vector.
		for (int x = 0; x < q; x++){
			for (u32 J = 0; J < V; J++){
				if ( (J & (1U << i)) || ((J & (1U << x)) == 0) ) wix[J] = 0;
				else wix[J] = w[i][J];
			}
			fzt(wix, wh[i][x]);
		}	
		fzt(w[i], wh[i][q]);	
	}
}	
	
	// 2. Compute f(R, U) for all U subset of V and all R subset of U.
	// That should take time O(4^q) using monotone arithmetic.
	// Compute f by dynamic programming. See Section 3.1.1 of Talvitie et al.

	// Initialize.
	Treal val1; val1 = 1.0;

	for (u32 U = 0; U < V; U++){
		u32 dU = proj16(U, U);
		for (u32 pR = 1; pR <= dU; pR++) f[U][pR] = 0;
		f[U][0] = 0; 
		f[U][dU] = val1; // Perhaps only this line is really needed. Reconsider.
	}

	// An alternative way to compute f. 
	// The key idea is the avoid precomputing what[i][T][pR].
	// Only what[i][x][cU] is needed, renamed as wh[i][x][cU].
	// Loop over subsets T of V.
	for (u32 T = 1; T < V - 1; T++){
		u32 cT = (V - 1) ^ T; // The complement of T.
		u32 dcT = proj16(cT, cT); // The complement of T in the form 00111111.
		u32 dT = proj16(T, T);
		
		// Initialize why[pR][s], where 
		// - pR ranges from 0 to the number of subsets of cT = (dcT + 1) ?
		// - s ranges from 0 to the number of 1-bits in T = |T|.
		// In fact, we only have to initialize for pR = 0, i.e., R = emptyset.
		for (int s = 0; s < popcount16(T); s++) why[0][s] = 0;

		// Loop over subsets R of cT in increasing order.
		u32 pR = 0;
		for (u32 cR = cT; cR >= 0; cR = (cR - 1) & cT, pR++){
			if (pR == 0) continue; // R is empty => no action.
			u32 R = cT ^ cR; 
			// Obtain pRx as the dense encoding of Rx := R \ {x} wrt cT for some arbitrary fixed element x of R.
			// Select the x in R that corresponds to the LSB of R. 
			int x = lsb16(R); u32 Rx = R ^ (1U << x);
			u32 pRx = proj16(cT, Rx); 

			// Loop over i in T and for each compute and store what[i](R, V \ T) =: why[pR][s], when i is the s-th element of T.
			// Use the recurrence what[i](R, V \ T) = what[i](R \ {x}, V \ T) + what[i]({x}, V \ T \ (R \ {x})).
			// Note: for what[i](R, U) we may assume that i in V \ U.
			// Note: V \ T \ (R \ {x}) = cT \ Rx
			int s = 0;
			for (u32 Y = T; Y > 0; Y &= (Y - 1), s++) why[pR][s] = why[pRx][s] + wh[lsb16(Y)][x][cT ^ Rx];
			
			// Loop over subsets S of T using the dense encoding dT.
			wprod[0] = val1;
			Treal sum; sum = 0;
			for (u32 S = 1; S <= dT; S++){ // We use directly the dense index; no need to use proj16(.,.).
				wprod[S] = why[pR][lsb16(S)] * wprod[S & (S - 1)]; // NOTE: lsb16(S) is in the projected range.
				sum += f[T][S] * wprod[S];
			}			
			// Store f(R, U). Note: U = T | R.
			u32 U = T | R;
			f[U][proj16(U, R)] = sum;
			if (cR == 0) break;
		}
	}

	// 3. Loop over U to peel off the parts of a root-partition.
	// That should take time O(q 2^q) using monotone arithmetic.
	u32 U = V-1; u32 Rk = 0; int k = 0;
	while (U){
		// Compute the weight of each subset R of U. Will draw R with probability proportional to the weight. 
		u32 dU = proj16(U, U);		
		for (u32 pR = 1; pR <= dU; pR++) weight[pR] = f[U][pR]; // Could simply use weight as a pointer instead of copying.
		u32 S = U; int s = 0; u32 R1tok = (V-1) ^ U; u32 pRk = proj16(R1tok, Rk);
		// For each s-th element i in U store what[i][][].
		if (Rk == 0){
			while (S){ whats[s++] = wh[lsb16(S)][q][R1tok]; S &= (S - 1); }
		} else {
			while (S){
				int i = lsb16(S); Treal tmp; tmp = 0; u32 L = R1tok;
				for (u32 Y = Rk; Y > 0; Y &= (Y - 1)){ int x = lsb16(Y); tmp += wh[i][x][L]; L ^= (1U << x); }
				whats[s++] = tmp; S &= (S - 1); 
			}
		}
		
		for (u32 pR = 1; pR <= dU; pR++){
			u32 S = pR; int s = 0;
			while (S){ if (S & 1){ weight[pR] = weight[pR] * whats[s]; } S >>= 1; s++; } 
		}		
		weight[0] = 0;

		// Draw R, i.e., R(k+1).
		u32 pR = sample(weight, 1, dU);
		// As pR is relative to U, we have to decompress it to be relative to V-1.
		u32 R = U; u32 r = dU;
		while (r > pR){ R = (R - 1) & U; r--; }
		
		// Don't update Rk yet, as we need the current value now.
		// Draw the parent set for each node in R. The parent set must intersect Rk.
		for (int i = 0; i < q; i++){
			// If node i is in R. Consider w[i][Hi] with Hi intersecting Rk to pick H[i] = Hi;
			if (R & (1 << i)) H[i] = isample(w[i], R1tok, Rk);
		}
		// Now update Rk to R and U to U \ R.
		Rk = R; k++; U ^= R; 
	}	
	return H;
}

#undef Treal
#undef vTrea
#undef vvTre
#undef vvvTr

