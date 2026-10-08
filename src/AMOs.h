#ifndef AMOS_H
#define AMOS_H

#include "Graph.h"
// #include "LehmerCode.h"
#include <algorithm>
#include <iostream>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace amos {

using u8 = uint8_t;
using u16 = uint16_t;
using u32 = uint32_t;
using Amo = int;
using AmoList = std::vector<Amo>;

#define SUBSETEQ(A, B)  (A == (A & B))

struct DD {
    u8 P[8] = { 0 };
    u8 A[8] = { 0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40 };
    u32 d = 0;
    u32 code = 0UL;

    bool anc(u8 i, u8 j) {
        return static_cast<bool>(A[j] & (1 << i));
    }

    DD add(u8 i, u8 j) {
        u8 I = (1 << i);
        u8 J = (1 << j);
        P[j] |= I;
        A[j] |= I;

        for (u8 v = 0; v < 7; v++) {
            if (A[v] & J) {
                A[v] |= A[i];
            }
        }

        return *this;
    }
};


class AMOs {
public:
  // Generating AMOs from https://arxiv.org/abs/2301.12212
  // static void generate(const Graph &G, std::vector<u_int16_t> &A,
  //                      std::vector<AmoKey> &amos, std::vector<int> to = {},
  //                      std::vector<int> membership = {}) {
  //   if (membership.empty()) {
  //     membership.resize(G.n, -1);
  //     for (int i = 0; i < A.size(); ++i) {
  //       for (int w = 0; w < G.n; ++w) {
  //         if (A[i] & (1 << w)) { // Check if w is in the subset
  //           membership[w] = i;
  //         }
  //       }
  //     }
  //   }

  //   int n = G.n;
  //   if (to.size() == n) {
  //     amos.push_back(LehmerCode::encode(to));
  //     return;
  //   }

  //   // Find the last non-empty subset in A
  //   int i = A.size() - 1;
  //   while (i >= 0 && A[i] == 0) {
  //     --i;
  //   }

  //   if (i < 0)
  //     return; // No valid subsets left

  //   int v = __builtin_ctz(
  //       A[i]); // Find the first set bit (smallest element in A[i])
  //   int x = v;
  //   u_int16_t R = 0;

  //   while (true) {
  //     A[i] &= ~(1 << x); // Remove x from A[i]
  //     membership[x] = -1;
  //     to.push_back(x);

  //     std::unordered_set<int> to_set(to.begin(), to.end());
  //     // Update neighbors of x
  //     std::vector<int> x_neighbors = G.get_neighbors(x);
  //     for (int w : x_neighbors) {
  //       if (to_set.count(w))
  //         continue;

  //       int j = membership[w];
  //       A[j] &= ~(1 << w); // Remove w from A[j]
  //       membership[w] = j + 1;

  //       if (j + 1 >= static_cast<int>(A.size()) || A[j + 1] == 0) {
  //         A.push_back(0);
  //       }

  //       A[j + 1] |= (1 << w); // Add w to A[j + 1]
  //     }

  //     generate(G, A, amos, to, membership);

  //     to_set.clear(); // Clear the set
  //     to_set.insert(to.begin(), to.end());
  //     // Revert changes for neighbors of x
  //     for (int w : x_neighbors) {
  //       if (to_set.count(w))
  //         continue;

  //       int j = membership[w];
  //       A[j] &= ~(1 << w); // Remove w from A[j]
  //       membership[w] = j - 1;
  //       A[j - 1] |= (1 << w); // Add w back to A[j - 1]
  //     }

  //     A[i] |= (1 << x); // Add x back to A[i]
  //     membership[x] = i;
  //     to.pop_back();

  //     if (x == v) {
  //       R = reachable_in_subset(G, v, A[i]);
  //       R &= ~(1 << v); // Remove v from the subset
  //     }

  //     if (R == 0)
  //       break;
  //     x = __builtin_ctz(R); // Get the next smallest element in R
  //     R &= ~(1 << x);       // Remove x from R
  //   }
  // }

  static void generate_alt(u32 G, AmoList& amos) {
    struct { u8 i, j; } arc[32]; int m = 0;
    u8 N[8] = {0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80}; // Neighbors.
    for (u8 j = 1; j < 7; j++) {
      for (u8 i = 0; i < j; i++, G >>= 1) {
        if (G & 1) {
          arc[m++] = {i, j};
          N[j] |= (1 << i);
          N[i] |= (1 << j);
        }
      }
    }

    DD stck[32];
    int t = 1;

    while (t) {
      DD F = stck[--t];
      u8 d = F.d;

      u8 i = arc[d].i, j = arc[d].j;

      if (!F.anc(j, i) && SUBSETEQ(F.P[j], N[i])) {
        DD L = F;
        L.d++;
        L.code |= 1UL << (i + ((j * (j - 1)) >> 1));
        if (d < m - 1) {
          stck[t++] = L.add(i, j);
        } else {
          amos.push_back(L.code);
        }
      }

      if (!F.anc(i, j) && SUBSETEQ(F.P[i], N[j])) {
        DD R = F;
        R.d++;
        if (d < m - 1) {
          stck[t++] = R.add(j, i);
        } else {
          amos.push_back(R.code);
        }
      }
    }
  }

private:
  // Compute reachable vertices in a subset
  static u_int16_t reachable_in_subset(const Graph &G, int v,
                                       u_int16_t subset) {
    const auto &adj = G.adj;
    u_int16_t visited = 0;
    std::stack<int> stack;
    stack.push(v);

    while (!stack.empty()) {
      int current = stack.top();
      stack.pop();

      if (visited & (1 << current))
        continue;
      visited |= (1 << current); // Mark current as visited

      for (int neighbor = 0; neighbor < G.n; ++neighbor) {
        if ((adj[current][neighbor]) && (subset & (1 << neighbor)) &&
            !(visited & (1 << neighbor))) {
          stack.push(neighbor);
        }
      }
    }

    return visited;
  }
};
} // namespace amos

#endif
