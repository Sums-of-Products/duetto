#ifndef GRAPH_H
#define GRAPH_H

#include <algorithm>
#include <iostream>
#include <queue>
#include <set>
#include <stack>
#include <vector>
#include <bitset>

namespace amos {

/**
 * A simple undirected graph structure with:
 * - adjacency list
 * - is_connected()
 * - is_chordal() using Maximum Cardinality Search (MCS)
 */

using AdjMatrix = std::vector<std::bitset<8>>;

struct Graph {
  int n;
  AdjMatrix adj;

  Graph(int n) : n(n), adj(n) {
    if (n > 8) {
      throw std::invalid_argument(
          "Graph size exceeds 8 vertices; adjust bitset size if needed.");
    }
  }

  std::vector<int> get_neighbors(int u) const {
    std::vector<int> neighbors;
    for (int v = 0; v < n; v++) {
      if (adj[u][v]) {
        neighbors.push_back(v);
      }
    }
    return neighbors;
  }

  void add_edge(int u, int v) {
    adj[u][v] = true;
    adj[v][u] = true;
  }

  // Remove an undirected edge (u, v)
  void remove_edge(int u, int v) {
    adj[u][v] = false;
    adj[v][u] = false;
  }

  bool edge_exist(int u, int v) const { return adj[u][v] == true; }

  // Check if the graph is connected via BFS/DFS
  bool is_connected() const {
    // If there are no vertices or only one vertex, it's trivially connected
    if (n <= 1)
      return true;

    // BFS/DFS to check connectivity
    std::vector<bool> visited(n, false);
    std::queue<int> q;

    // Start from the first vertex
    visited[0] = true;
    q.push(0);

    int countVisit = 1; // Count of visited vertices

    while (!q.empty()) {
      int u = q.front();
      q.pop();

      // Iterate over all vertices to check adjacency
      for (int v = 0; v < n; v++) {
        if (adj[u][v] &&
            !visited[v]) { // Check if there's an edge and if it's unvisited
          visited[v] = true;
          countVisit++;
          q.push(v);
        }
      }
    }

    // The graph is connected if all vertices are visited
    return (countVisit == n);
  }

  // Check if the graph is chordal using Maximum Cardinality Search (MCS).
  // If MCS yields a valid perfect elimination ordering, the graph is chordal.
  bool is_chordal() const {
    std::vector<int> weight(n, 0);
    std::vector<bool> chosen(n, false);
    std::vector<int> peo(n, -1);

    for (int step = n - 1; step >= 0; step--) {
      int maxW = -1, v = -1;
      for (int i = 0; i < n; i++) {
        if (!chosen[i] && weight[i] > maxW) {
          maxW = weight[i];
          v = i;
        }
      }
      if (v == -1)
        return false;

      peo[step] = v;
      chosen[v] = true;

      for (int w = 0; w < n; w++) {
        if (adj[v][w] && !chosen[w]) {
          weight[w]++;
        }
      }
    }

    std::vector<int> posInOrder(n);
    for (int i = 0; i < n; i++) {
      posInOrder[peo[i]] = i;
    }

    for (int i = 0; i < n; i++) {
      int v = peo[i];
      std::vector<int> higherNeighbors;
      for (int w = 0; w < n; w++) {
        if (adj[v][w] && posInOrder[w] > i) {
          higherNeighbors.push_back(w);
        }
      }

      for (size_t a = 0; a < higherNeighbors.size(); a++) {
        for (size_t b = a + 1; b < higherNeighbors.size(); b++) {
          if (!adj[higherNeighbors[a]][higherNeighbors[b]]) {
            return false;
          }
        }
      }
    }

    return true;
  }
  // Compress the upper triangular matrix into a u_int16
  u_int32_t get_upper_triangle() const {
    u_int32_t compressed = 0;
    int total_bits =
        (n * (n - 1)) / 2; 
    // Start from the most significant bit
    int bit_position = total_bits - 1;

    for (int i = 0; i < n; ++i) {
      for (int j = i + 1; j < n; ++j) {
        if (adj[i][j]) {
          compressed |= (1 << bit_position);
        }
        --bit_position; 
      }
    }

    return compressed;
  }

};

struct AdjMatrixHash {
  size_t operator()(const AdjMatrix &matrix) const {
    size_t seed = 0;

    // Hash each row (bitset) in the matrix
    for (const auto &row : matrix) {
      // Convert the bitset to its integer representation and hash it
      seed ^= std::hash<std::bitset<8>>{}(row) + 0x9e3779b9 + (seed << 6) +
              (seed >> 2);
    }

    return seed;
  }
};

} // namespace amos
#endif
