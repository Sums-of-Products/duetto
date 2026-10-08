
#ifndef AMOC
#define AMOC

#include "AMOs.h"
#include "Graph.h"
#include "UCCG.h"
#include <iostream>
#include <map>
#include <set>
#include <vector>

#define AMO_MAXSIZE 7

using GraphKey = uint32_t;
namespace amos {

const int maxsize = AMO_MAXSIZE;

class AMOcache {
public:
  AMOcache(int maxk) : L(maxk + 1), max_k(maxk) {
    max_k = maxk;
  }
  void gen_amos(){
  	if (all_done) return;
  	using namespace std;
    cerr << "Generating AMOs for k = ";
    for (int k = 2; k <= max_k; ++k) {
    	cerr << k << "...";
    	generate_for_k(k);
    }
    cerr << "AMO generation completed for all k up to " << max_k << "." << endl;
    all_done = true;
  }	

  AmoList& get_amos(int k, GraphKey key) { return L[k][key]; }
  Amo get_amo(int k, GraphKey key, int t) { return L[k][key][t]; }
  int amos_count(int k, GraphKey key) { return L[k][key].size(); }

  static void print_upper_triangle(GraphKey compressed, int n) {
    using namespace std;
    int total_bits = (n * (n - 1)) / 2;
    int bit_position = total_bits - 1;

    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        if (j > i) {
          cout << ((compressed & (1 << bit_position)) ? "1 " : "0 ");
          --bit_position;
        } else {
          cout << "0 ";
        }
      }
      cout << endl;
    }
  }

  void print() const {
    using namespace std;
    for (int k = 0; k < max_k; ++k) {
      cout << "k = " << k << ":\n";

      for (GraphKey G = 0; G < L[k].size(); ++G) {
        auto &amos = L[k][G]; // Access the AmoList at index G

        if (amos.empty())
          continue;

        print_upper_triangle(G, k);

        cout << "AMOs:\n";
        for (const auto &amo : amos) {
          cout << amo << "\n";
        }
      }
    }
  }

  std::vector<int> total_count() {
    std::vector<int> counts(max_k + 1, 0);

    for (int k = 2; k <= max_k; ++k) { 
      for (const auto &amos : L[k]) {
        counts[k] += amos.size(); 
      }
    }

    return counts;
  }

private:
  // AMO encoded to lehmer
  std::vector<std::vector<AmoList>> L;

  int max_k;
  bool all_done = false;

  void generate_for_k(int k) {
    auto graphs = UCCG::generate(k); // Generate UCCGs for k

    std::vector<AmoList> amos_list;
    amos_list.resize(1 << (k * (k - 1)) / 2);

    for (size_t idx = 0; idx < graphs.size(); ++idx) {
      Graph &G = graphs[idx];

      std::vector<u_int16_t> A;
      A.push_back(0); // A[0] initially contains all vertices
      for (int i = 0; i < G.n; ++i) {
        A[0] |= (1 << i);
      }

      AmoList amos;
      std::vector<int> to;

      GraphKey upper_triangle = G.get_upper_triangle();
      AMOs::generate_alt(upper_triangle, amos);

      // Store the AMOs for the current graph
      amos_list[upper_triangle] = std::move(amos);
    }

    // Store the results for this value of k
    L[k] = std::move(amos_list);
  }
};

// Create an AMOcache object. You should not need many of these.
AMOcache amocache(maxsize);


} // namespace amos

#undef AMO_MAXSIZE
#endif

