#ifndef SCOREMANAGER_HPP
#define SCOREMANAGER_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <tuple>
#include <algorithm>
#include <regex>
#include <stdexcept>
#include <iomanip>

class ScoreManager {
public:
    // Custom hash function for std::vector<int>
    struct VectorHash {
        size_t operator()(const std::vector<int>& v) const {
            size_t hash_value = 0;
            for (const int& elem : v) {
                hash_value ^= std::hash<int>{}(elem) + 0x9e3779b9 + (hash_value << 6) + (hash_value >> 2);
            }
            return hash_value;
        }
    };

    using CacheType = std::unordered_map<int, std::unordered_map<std::vector<int>, double, VectorHash>>;

    // Constructor to initialize the cache
    ScoreManager(const std::string& path, int max_parents, int n_nodes) {
        cache = cache_f(path, max_parents, n_nodes);
    }

    // Getter for cache
    const CacheType& get_cache() const {
        return cache;
    }

    // Method to get the score for a specific node and parent set
    double get_score(int node, const std::vector<int>& parents) const {
        auto node_it = cache.find(node);
        if (node_it != cache.end()) {
            const auto& parents_map = node_it->second;
            auto parents_it = parents_map.find(parents);
            if (parents_it != parents_map.end()) {
                return parents_it->second; // Return the score
            } else {
                throw std::runtime_error("Parent set not found for node: " + std::to_string(node));
            }
        } else {
            throw std::runtime_error("Node not found: " + std::to_string(node));
        }
    }

private:
    CacheType cache;

    // Function to parse the file and initialize the cache
    CacheType cache_f(const std::string& path, int max_parents, int n_nodes) {
        CacheType cache;
        std::ifstream file(path);
        if (!file.is_open()) {
            throw std::runtime_error("Unable to open file: " + path);
        }

        std::vector<int> vertices;
        std::vector<int> n_parents_set;
        std::string line;
        int t = 1;
        std::vector<std::string> filedata;

        while (std::getline(file, line)) {
            filedata.push_back(line);
        }

        // Parse vertices
        for (int i = 0; i < n_nodes; ++i) {
            line = filedata[t];
            auto numbers = extract_numbers(line);
            vertices.push_back(static_cast<int>(numbers[0])); // First number is the vertex ID
            t += static_cast<int>(numbers[1]) + 1;           // Skip the parent set lines
        }

        // Parse cache
        t = 1;
        for (int i : vertices) {
            line = filedata[t];
            auto numbers = extract_numbers(line);
            int now_n_parents_set = static_cast<int>(numbers[1]);
            n_parents_set.push_back(now_n_parents_set);

            std::unordered_map<std::vector<int>, double, VectorHash> now_dict;

            for (int j = 0; j < now_n_parents_set; ++j) {
                line = filedata[t + j + 1];
                auto numbers = extract_numbers(line);

                double score = numbers[0]; // First number is the score
                int num_parents = static_cast<int>(numbers[1]); // Second number is the number of parents

                // Extract the parent set (starting from the third number)
                std::vector<int> now_parents(numbers.begin() + 2, numbers.end());

                if (now_parents.size() != static_cast<size_t>(num_parents)) {
                    throw std::runtime_error("Mismatch between declared and actual number of parents.");
                }

                std::sort(now_parents.begin(), now_parents.end());
                now_dict[now_parents] = score;
            }

            cache[i] = now_dict;
            t += now_n_parents_set + 1;
        }

        return cache;
    }

    // Helper function to extract numbers from a line (handles both integers and floating-point numbers)
    std::vector<double> extract_numbers(const std::string& line) {
        std::vector<double> numbers;
        std::regex number_regex("-?\\d*\\.?\\d+");
        std::sregex_iterator begin(line.begin(), line.end(), number_regex);
        std::sregex_iterator end;

        for (auto it = begin; it != end; ++it) {
            std::string num_str = it->str();
            
            try {
                if (!num_str.empty()) {
                    double num = std::stod(num_str);
                    numbers.push_back(num);
                }
            } catch (const std::invalid_argument& e) {
                std::cerr << "Invalid argument: " << num_str << " cannot be converted to a number.\n";
            } catch (const std::out_of_range& e) {
                std::cerr << "Out of range: " << num_str << " is too large for a number.\n";
            }
        }

        return numbers;
    }
};

#endif // SCOREMANAGER_HPP
