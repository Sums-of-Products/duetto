#include <iostream>
#include <getopt.h>
#include <vector>
#include <string>
#include <filesystem>
#include <random>
#include <chrono>
#include <cmath>
#include <algorithm>

#include "src/ec.hpp"
#include "src/ac.hpp"
#include "src/usr.hpp"
#include "src/bg.hpp"
#define BD_ID	1	

vint getp(int i, DDAG dagp) { return dagp.P[i]; }

double gets(int i, const vint& P, Scorer_ptr scop) {
    return scop->sfam(i, P);
}

double score(DDAG dagp, Scorer_ptr scop, int n) {
    double S = 0;
    for(int i = 0; i < n; i++) S += gets(i, getp(i, dagp), scop);
    return S;
}


Scorer_ptr initialize_scorer(
    std::string filepath,
    int mid1,
    int mid2,
    int er,
    int prior,
    int prune,
    int maxcan
) {
    Scorer_ptr scop = nullptr;

    switch (1){
 		case BD_ID: scop = (Scorer_ptr)(new BD()); break;
 		default: cerr << "\t ERROR: DAGsampler::init: id = " << 1 << " not implemented; exit.\n"; exit(1); 
    }
    if (filepath.size() >= 4 && filepath.substr(filepath.size() - 4) == ".jkl") {
        scop->rjkl(filepath);
    } else if (filepath.size() >= 4 && filepath.substr(filepath.size() - 4) == ".csv") { 
		scop->dear(filepath);
   	} else {
        scop->read(filepath);
    }

    int pru = prune % 10;
    // double build_param = prune / 10;
    double eps = 1.0 / (1ULL << er);

    scop->set_mid(mid1,mid2);
    scop->set_eps(eps);
    scop->set_pri(prior);
    scop->set_pru(pru);
    scop->set_can(maxcan);
    scop->build();
//	scop->pjkl("mytest.jkl"); // [Mikko]

    return scop;
}

void special_test(){
	double z = 5.432;
	B2real x; x = z;
	cerr << "Testing: double z = " << z << ", Breal x = z after transformations = " << x.get_double() << endl; 

	double w = 1.234;
	B2real y; y = w;
	cerr << "Testing: double w = " << w << ", Breal y = w after transformations = " << y.get_double() << endl; 

	B2real p; p = x * y;
	cerr << "Testing: z * w = " << z * w << ", x * y = " << p.get_double() << endl; 

	exit(1);
}

int main(int argc, char* argv[]) {
	//special_test();

    if (argc < 2) {
        std::cerr << "\nUsage: " << argv[0] << " datafile [options]\n"
                  << "Options:\n"
                  << "  -m <int>   max parents (default 4)\n"
                  << "  -p <int>   prior       (default 1)\n"
                  << "  -P <int>   prune       (default 2)\n"
                  << "  -K <int>   maxcan      (default 64)\n"
                  << "  -i <int>   iterations  (default 10000)\n"
                  << "  -q <int>   q block size (default min(floor(sqrt(n)),10))\n"  //
                  << "  -R <int>   random seed (default 13)\n"
                  << "  -USRrev <n>   MH-block repeat count (default 100)\n"
                  << "  -EC <n>    move1 repeat count    (default 100)\n"
                  << "  -AC <n>    move2 repeat count    (default 100)\n"
                  << "  -iPPE <n>  save posterior matrix every n iterations\n\n";
        return 1;
    }

    std::string filepath = argv[1];
    std::string base_name = filepath.substr(filepath.find_last_of("/\\") + 1);
    size_t lastdot = base_name.find_last_of(".");
    if (lastdot != std::string::npos)
        base_name = base_name.substr(0, lastdot);

    int mid1 = 4, mid2 = 4;   
    int er = 10, prior = 1, prune = 2, maxcan = 64;
    int iterations = 10000, q = -1;
    int burnin = -1;
    int seed = rd();
    //int seed = 13; // [Mikko]

    int USRrev = 100, EC = 25, AC = 25, BG = 0;
    int USRrnd = 50, USRid = 0; 
    int iPPE = 0; 

    for (int a = 2; a < argc; a++) {
        std::string arg = argv[a];
        if (arg == "-m") {
            mid1 = std::atoi(argv[++a]);
        if (a + 1 < argc && argv[a + 1][0] != '-') {
            mid2 = std::atoi(argv[++a]);
        } else {
            mid2 = mid1;
        }
    }
        else if (arg == "-p") prior = atoi(argv[++a]);
        else if (arg == "-P") prune = atoi(argv[++a]);
        else if (arg == "-K") maxcan = atoi(argv[++a]);
        else if (arg == "-i") iterations = atoi(argv[++a]);
        else if (arg == "-q") q = atoi(argv[++a]);
        else if (arg == "-R") seed = atoi(argv[++a]);
        else if (arg == "-USRrev") USRrev = atoi(argv[++a]);
        else if (arg == "-EC") EC = atoi(argv[++a]);
        else if (arg == "-AC") AC = atoi(argv[++a]);
        else if (arg == "-BG") BG = atoi(argv[++a]);       
        else if (arg == "-USRrnd") USRrnd = atoi(argv[++a]);
        else if (arg == "-USRid") USRid = atoi(argv[++a]);
        else if (arg == "-burnin") burnin = atoi(argv[++a]);
        else if (arg == "-a") er = atoi(argv[++a]);
        else if (arg == "-iPPE") iPPE = atoi(argv[++a]); 
        else {
            std::cerr << "Unknown option: " << arg << std::endl;
            return 1;
        }
    }

    if (burnin < 0) burnin = iterations / 10;  

    gen.seed(seed);

    std::filesystem::create_directories("results");

    std::string moves_tag = "";
    if (USRrev > 0)     moves_tag += (moves_tag.empty() ? "USRrev"     : "-USRrev");
    if (USRrnd > 0) moves_tag += (moves_tag.empty() ? "USRrnd" : "-USRrnd");
    if (USRid > 0) moves_tag += (moves_tag.empty() ? "USRid" : "-USRid");
    if (BG  > 0)     moves_tag += (moves_tag.empty() ? "BG"      : "-BG");
    if (EC  > 0)     moves_tag += (moves_tag.empty() ? "EC"      : "-EC");
    if (AC  > 0)     moves_tag += (moves_tag.empty() ? "AC"      : "-AC");
    if (moves_tag.empty()) moves_tag = "NONE";

    std::string score_file = "results/" + base_name + "_s_seed" + std::to_string(seed) + "_" + moves_tag + ".txt";
    std::string post_file  = "results/" + base_name + "_p_seed" + std::to_string(seed) + "_" + moves_tag + ".txt";
    std::string ppe_file_name = "results/" + base_name + "_PPE_seed" + std::to_string(seed) + "_" + moves_tag + ".txt";
    std::string dags_file = "results/" + base_name + "_dags_seed" + std::to_string(seed) + "_" + moves_tag + ".txt";


    std::ofstream score_out(score_file, std::ios::out | std::ios::trunc);
    if (!score_out) {
        std::cerr << "Error opening score file: " << score_file << std::endl;
        return 1;
    }
    score_out << std::fixed << std::setprecision(6);

    std::ofstream ppe_file;
    if (iPPE > 0) {
        ppe_file.open(ppe_file_name, std::ios::out | std::ios::trunc);
        if (!ppe_file) {
            std::cerr << "Error opening PPE file: " << ppe_file_name << std::endl;
            return 1;
        }
    }

    std::ofstream dags_out(dags_file, std::ios::out | std::ios::trunc);
    if (!dags_out) {
        std::cerr << "Error opening dags file: " << dags_file << std::endl;
        return 1;
    }

    auto t_start = std::chrono::high_resolution_clock::now();
    auto t_precomp_end = t_start; 
    auto t_sampling_end = t_start;

    Scorer_ptr scop = initialize_scorer(filepath, mid1, mid2, er, prior, prune, maxcan);
    int n_nodes = scop->n;
    if (q < 0)
        q = std::min(static_cast<int>(std::sqrt(n_nodes)), 10);

    t_precomp_end = std::chrono::high_resolution_clock::now();

    DDAG dag(n_nodes, false);
    std::vector<std::vector<int>> adj(n_nodes, std::vector<int>(n_nodes, 0));

        std::cout << "\n══════════════════════════ Sampling parameters ══════════════════════════\n"
              << "Data file:        " << filepath << "\n"
              << "Nodes:            " << n_nodes << "\n"
              << "Seed:             " << seed << "\n"
              << "Iterations:       " << iterations << "\n"
              << "Burn-in:          " << burnin << "\n"
              << "Mid 1:            " << mid1 << "\n"
              << "Mid 2:            " << mid2 << "\n"
              << "Max candidates:   " << maxcan << "\n"
              << "Prior:            " << prior << "\n"
              << "Prune error:      " << er << "\n"
              << "Prune type:       " << prune << "\n"
              << "Block size:       " << q << "\n"
              << "USRrev repeats:   " << USRrev << "\n"
              << "USRrnd repeats:   " << USRrnd << "\n"
              << "USRid repeats:    " << USRid << "\n"
              << "BG repeats:       " << BG << "\n"
              << "EC repeats:       " << EC << "\n"
              << "AC repeats:       " << AC << "\n"
              << "PPE interval:     " << iPPE << "\n"
              << "═════════════════════════════════════════════════════════════════════════\n";

    if (burnin > 0) {
        for(int it = 0; it < burnin; it++) {
            for (int k = 0; k < USRrev; k++) usr(dag, scop, q, mid1, n_nodes, 0);
            for (int k = 0; k < USRrnd; k++) usr(dag, scop, q, mid1, n_nodes, 1);
            for (int k = 0; k < USRid; k++) usr(dag, scop, q, mid1, n_nodes, 2);
            for(int k = 0; k < BG;  k++) /*bg*/ bg(dag, scop, q);         
            for(int k = 0; k < EC;  k++) ec(dag, scop, mid1, n_nodes);
            for(int k = 0; k < AC;  k++) ac(dag, scop, mid1, n_nodes);

            if (it % 100 == 0) {
                std::cout << "\rburn-in " << it << "/" << burnin << std::flush;
            }
        }
        std::cout << "\rburn-in " << burnin << "/" << burnin << std::endl;
    }

    for(int iter = 0; iter < iterations; iter++) {

        if (iter %  100 == 0) {
            std::cout << "\riter " << iter << "/" << iterations << std::flush;
        }
        for (int k = 0; k < USRrev; k++) usr(dag, scop, q, mid1, n_nodes, 0);
        for (int k = 0; k < USRrnd; k++) usr(dag, scop, q, mid1, n_nodes, 1);
        for (int k = 0; k < USRid; k++) usr(dag, scop, q, mid1, n_nodes, 2);
        for(int k = 0; k < BG;  k++) /*bg*/ bg(dag, scop, q);
        for(int k = 0; k < EC; k++) ec(dag, scop, mid1, n_nodes);
        for(int k = 0; k < AC; k++) ac(dag, scop, mid1, n_nodes);

        double sc = score(dag, scop, n_nodes);
        score_out << sc << "\n";

        add_post(dag, n_nodes, adj);

        if (iPPE > 0 && (iter + 1) % iPPE == 0) {
            append_adj_matrix(ppe_file, adj, iter + 1);
        }
        dag.printDAG(dags_out);
    }

    std::cout << "\riter " << iterations << "/" << iterations << std::endl;

    t_sampling_end = std::chrono::high_resolution_clock::now();

    save_matrix_to_file(adj, iterations, post_file);

    score_out.close();
    dags_out.close();

    if (iPPE > 0) ppe_file.close();

    std::cout << "Saved:\n  " << score_file << "\n  " << post_file << "\n  " << dags_file;  
    if (iPPE > 0) std::cout << "\n  " << ppe_file_name;
    std::cout << "\nDone.\n";

    double precomp_sec = std::chrono::duration<double>(t_precomp_end - t_start).count();
    double sampling_sec = std::chrono::duration<double>(t_sampling_end - t_precomp_end).count();

    std::ofstream post_append(post_file, std::ios::app);
    post_append << "\n# Pre-computation time (seconds): " << std::fixed << std::setprecision(4) << precomp_sec;
    post_append << "\n# Sampling time (seconds): " << std::fixed << std::setprecision(4) << sampling_sec;
    post_append << "\n";
    post_append.close();

    return 0;
}
