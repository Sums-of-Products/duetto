//
// 27 Nov 2024, Mikko Koivisto
//

#ifndef BASICS_HPP
#define BASICS_HPP

#include <iostream>
#include <fstream> 
#include <sstream>
#include <string>
#include <map>
#include <functional> 
#include <iomanip> 
#include <vector>
#include <bit>
#include <bitset>
#include <memory>		// shared_ptr
#include <limits>
#include <random>
#include <boost/dynamic_bitset.hpp>
//#include <immintrin.h>

// MACROS

#define FIXED_F(x, w, p)	std::setw(w) << std::right << std::fixed << std::setprecision(p) << (x)
#define FIXED_FLOAT(x)		std::setw(5) << std::right << std::fixed << std::setprecision(0) << (x)
#define FIXED_INT(x)		std::setw(5) << std::right << (x)

// TYPES ETC. KEYWORD "using" PRFERRED OVER the OLDER "typedef"

using string	= std::string;
using ostream	= std::ostream;
using ifstream	= std::ifstream;
using ofstream	= std::ofstream;
using sstream	= std::stringstream;

using vboo		= std::vector<bool>;
using vint		= std::vector<int>;
using vdou		= std::vector<double>;
using vstr		= std::vector<string>;
using vvdo		= std::vector<vdou>;

using u8		= uint8_t;
using u16		= uint16_t;
using u32		= uint32_t;
using u64		= uint64_t;

using b64		= std::bitset<  64>;
using b256		= std::bitset< 256>;
using b1024		= std::bitset<1024>;

using dbitset	= boost::dynamic_bitset<>;

struct plon { long a; long b; };

using std::cout;	using std::cerr;	using std::endl;

// CONSTANTS

const double infdouble = std::numeric_limits<double>::infinity();

// GLOBAL VARIABLES

bool debug = false;

// RANDOM NUMBERS

std::random_device rd;  					// Will be used to obtain a seed for the random number engine.
std::mt19937 gen(13); //rd()); 			// Standard mersenne_twister_engine seeded with rd().
std::uniform_real_distribution<> myunif(0.0, 1.0);

inline double urnd()			{ return (double) myunif(gen); }
inline int rint(int b) 		{ return std::min( (int) std::floor( urnd() * (b + 1) ), b ); } // Range: 0, 1, ..., b.
inline int rin2(int a, int b){ return std::min( (int) std::floor( a + urnd() * (b - a + 1) ), b ); } // Range: a, a+1, ..., b.
inline int geom(double p)	{ return (int) std::floor(std::log( urnd() ) / std::log(1 - p) ); }

// BITSET MANIPULATION

template<typename T>
T make_bitset(const vint& set){ T b = 0; for (auto x : set) b[x] = 1; return b; }

template<typename T>
vint make_vint(const T b, int n){ vint set; set.reserve(n); for (int i = 0; i < n; i++) if (b[i]) set.push_back(i); return set; }

vint make_vint(u32 b){
	vint v; v.reserve(32); 
	int i = 0; 
	while (b){ if (b & 1){ v.push_back(i); } b >>= 1; i++; } 
	return v; 
}

vint make_vint(const dbitset& b){
	vint v; v.reserve(b.count());
	size_t i = b.find_first();
	while (i != boost::dynamic_bitset<>::npos){ v.push_back(i); i = b.find_next(i); }
	return v;
}

dbitset make_dbitset(vint& v, int n){ dbitset b(n); for (auto i : v) b.set(i); return b; }

inline int bitcount(u64 x){ int c = 0; while (x){ c += x & (u64) 1; x >>= 1;} return c; }	// Not super fast. 

inline u32 popcount32(u32 x){
	return std::popcount(x); // Assumes C++20 available.
}

inline u32 popcount16(u32 x){
	return std::popcount(x); // Assumes C++20 available.
/*
    // Stage 1: Count pairs of bits
    x = x - ((x >> 1) & 0x5555);
    // Stage 2: Count nibbles (4-bit groups)
    x = (x & 0x3333) + ((x >> 2) & 0x3333);
    // Stage 3: Count bytes (8-bit groups) and aggregate
    x = (x + (x >> 4)) & 0x0F0F;  
    // Stage 4: Sum the two lower byte totals together
    return (x + (x >> 8)) & 0x00FF;
*/
}

inline int lsb16(u32 x){
	return std::countr_zero(x); // Assumes C++20 available.
/*	
    if (x == 0) return -1;
    int index = 0;
    // Use the fact that x is under 16 bits to do a quick binary split
    if ((x & 0x00FF) == 0) { x >>= 8; index += 8; }
    if ((x & 0x000F) == 0) { x >>= 4; index += 4; }
    if ((x & 0x0003) == 0) { x >>= 2; index += 2; }
    if ((x & 0x0001) == 0) {          index += 1; }
    return index;
*/
}

// Project S on T, i.e., remove in-between bits not set to 1 in T.
// Example: T = 01001101, S = 01000100 => return 1010.
inline u32 proj16(u32 T, u32 S){
	//return _pext_u32(S, T); // Assumes intrinsics are available.
	 uint32_t res = 0;
    // Iterate only over the set bits of the mask
    for (uint32_t dest_bit = 1; T != 0; dest_bit <<= 1) {
        // Extract the lowest set bit of the mask
        uint32_t low_bit = T & -T; 
        if (S & low_bit) {
            res |= dest_bit;
        }
        T ^= low_bit; // Clear the lowest set bit
    }
    return res;
}

// PRINTING 

void print(ostream& os, const vint& v){
	using namespace std;
	os << " (";
	for (auto x : v) os << " " << x;
	os << ")" << endl;
}
void print(ostream& os, const vint& v, double w){
	using namespace std;
	os << "\t weight: " << FIXED_F(w, 7, 1) << ", vec: ";
	for (auto x : v) os << " " << x;
	os << endl;
}

// PERMUTATIONS

// Shuffles A.
void rperm(vint& A){
	int n = A.size();
	for (int i = 1; i < n; i++){
		int j = rint(i);
		auto tmp = A[i]; A[i] = A[j]; A[j] = tmp;
	}
}

// Return a random k-permutation of {0, 1, ..., n-1}.
vint rperm(int n, int k){
    vint A(n);
    iota(A.begin(), A.end(), 0); // NOTE: I dislike this cumbersome iterator syntax...
    for (int i = 0; i < k; ++i) {
    	int j = rin2(i, n-1);
    	auto tmp = A[i]; A[i] = A[j]; A[j] = tmp;
    }
    A.resize(k);
    return A;
}

// ORDERED ARRAYS

inline void insert_sorted(int item, const vint& A, vint& B){
	int d = A.size(); B.resize(1 + d);
	int j = 0; while (j < d && A[j] < item){ B[j] = A[j]; j++; } B[j] = item; while (j < d){ B[j + 1] = A[j]; j++; } 
}

inline void add_ord_strict_db(int v, vint& L){
	std::cerr << "\t add_ord_strict_db: v = " << v << ", L.size() = " << L.size(); print(std::cerr, L);
	int s = L.size(); 
	int t = s; while (t > 0 && L[t - 1] > v){ t--; }
	std::cerr << "\t add_ord_strict_db: t = " << t; print(std::cerr, L);
	if (t && L[t-1] == v) return; // No change, as v was already in L.
	std::cerr << "\t add_ord_strict_db: t = " << t; print(std::cerr, L);	
	L.resize(s+1);
	while (s > t){ L[s] = L[s-1]; s--; } L[t] = v;
}

inline void add_ord_strict(int v, vint& L){
	//std::cerr << "\t add_ord_strict: v = " << v << ", L.size() = " << L.size(); print(std::cerr, L);
	int s = L.size(); 
	int t = s; while (t > 0 && L[t-1] > v){ t--; }
	if (t && L[t - 1] == v) return; // No change, as v was already in L.
	L.resize(s+1);
	while (s > t){ L[s] = L[s-1]; s--; } L[t] = v;
}

inline void rem_ord_strict(int v, vint& L){
	int s = L.size(); 
	int t = 0; while (t < s && L[t] < v){ t++; }
	if (t == s) return;			// No change, as v was not in L.
	while (t < s-1){ L[t] = L[t+1]; t++; }
	L.resize(s-1);	
}

// HASHING

class vinthas {
	public:
		std::size_t operator()(const vint& vec) const {
  			std::size_t seed = vec.size();
  			for(auto& i : vec){ seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2); }
  			return seed;
		}
};

// BIG DATA

// Dags on 2, 3, and 4 nodes:

u8 d2[3][2] = {
	{	0,	0},		// (empty)
	{	2,	0},		// 0 <- 1
	{	0,	1}};	// 1 <- 0	

u8 d3[25][3] = {
	{	0, 	0, 	0},	// (empty)
	{	2, 	0, 	0}, // 0<-1  
	{	4, 	0, 	0}, // 0<-2  
	{	0, 	1, 	0}, // 1<-0
	{	0, 	4, 	0}, // 1<-2
	{	0, 	0, 	1}, // 2<-0
	{	0, 	0, 	2}, // 2<-1
	{	6, 	0, 	0}, // 0<-1,2 
	{	2, 	4, 	0}, // 0<-1 1<-2
	{	2, 	0, 	1}, // 0<-1 2<-0
	{	2, 	0, 	2}, // 0<-1 2<-1
	{	4, 	4, 	0}, // 0<-2 1<-2
	{	4, 	0, 	2}, // 0<-2 2<-1
	{	4, 	1, 	0}, // 0<-2 1<-0
	{	0, 	5, 	0}, // 1<-0,2
	{	0, 	1, 	1}, // 1<-0 2<-0
	{	0, 	1, 	2}, // 1<-0 2<-1
	{	0, 	4, 	1}, // 1<-2 2<-0
	{	0, 	0, 	3}, // 2<-0,1
	{	6, 	4, 	0}, // 0 < 1 < 2 
	{	6, 	0, 	2}, // 0 < 2 < 1
	{	4, 	5, 	0}, // 1 < 0 < 2
	{	0, 	5, 	1}, // 1 < 2 < 0
	{	2, 	0, 	3}, // 2 < 0 < 1
	{	0, 	1, 	3}};// 2 < 1 < 0

u8 d4[][4] = {
	{ 0, 0, 0, 0 },
	{ 0, 0, 0, 1 },
	{ 0, 0, 0, 2 },
	{ 0, 0, 0, 3 },
	{ 0, 0, 0, 4 },
	{ 0, 0, 0, 5 },
	{ 0, 0, 0, 6 },
	{ 0, 0, 0, 7 },
	{ 0, 0, 1, 0 },
	{ 0, 0, 1, 1 },
	{ 0, 0, 1, 2 },
	{ 0, 0, 1, 3 },
	{ 0, 0, 1, 4 },
	{ 0, 0, 1, 5 },
	{ 0, 0, 1, 6 },
	{ 0, 0, 1, 7 },
	{ 0, 0, 2, 0 },
	{ 0, 0, 2, 1 },
	{ 0, 0, 2, 2 },
	{ 0, 0, 2, 3 },
	{ 0, 0, 2, 4 },
	{ 0, 0, 2, 5 },
	{ 0, 0, 2, 6 },
	{ 0, 0, 2, 7 },
	{ 0, 0, 3, 0 },
	{ 0, 0, 3, 1 },
	{ 0, 0, 3, 2 },
	{ 0, 0, 3, 3 },
	{ 0, 0, 3, 4 },
	{ 0, 0, 3, 5 },
	{ 0, 0, 3, 6 },
	{ 0, 0, 3, 7 },
	{ 0, 0, 8, 0 },
	{ 0, 0, 8, 1 },
	{ 0, 0, 8, 2 },
	{ 0, 0, 8, 3 },
	{ 0, 0, 9, 0 },
	{ 0, 0, 9, 1 },
	{ 0, 0, 9, 2 },
	{ 0, 0, 9, 3 },
	{ 0, 0, 10, 0 },
	{ 0, 0, 10, 1 },
	{ 0, 0, 10, 2 },
	{ 0, 0, 10, 3 },
	{ 0, 0, 11, 0 },
	{ 0, 0, 11, 1 },
	{ 0, 0, 11, 2 },
	{ 0, 0, 11, 3 },
	{ 0, 1, 0, 0 },
	{ 0, 1, 0, 1 },
	{ 0, 1, 0, 2 },
	{ 0, 1, 0, 3 },
	{ 0, 1, 0, 4 },
	{ 0, 1, 0, 5 },
	{ 0, 1, 0, 6 },
	{ 0, 1, 0, 7 },
	{ 0, 1, 1, 0 },
	{ 0, 1, 1, 1 },
	{ 0, 1, 1, 2 },
	{ 0, 1, 1, 3 },
	{ 0, 1, 1, 4 },
	{ 0, 1, 1, 5 },
	{ 0, 1, 1, 6 },
	{ 0, 1, 1, 7 },
	{ 0, 1, 2, 0 },
	{ 0, 1, 2, 1 },
	{ 0, 1, 2, 2 },
	{ 0, 1, 2, 3 },
	{ 0, 1, 2, 4 },
	{ 0, 1, 2, 5 },
	{ 0, 1, 2, 6 },
	{ 0, 1, 2, 7 },
	{ 0, 1, 3, 0 },
	{ 0, 1, 3, 1 },
	{ 0, 1, 3, 2 },
	{ 0, 1, 3, 3 },
	{ 0, 1, 3, 4 },
	{ 0, 1, 3, 5 },
	{ 0, 1, 3, 6 },
	{ 0, 1, 3, 7 },
	{ 0, 1, 8, 0 },
	{ 0, 1, 8, 1 },
	{ 0, 1, 8, 2 },
	{ 0, 1, 8, 3 },
	{ 0, 1, 9, 0 },
	{ 0, 1, 9, 1 },
	{ 0, 1, 9, 2 },
	{ 0, 1, 9, 3 },
	{ 0, 1, 10, 0 },
	{ 0, 1, 10, 1 },
	{ 0, 1, 10, 2 },
	{ 0, 1, 10, 3 },
	{ 0, 1, 11, 0 },
	{ 0, 1, 11, 1 },
	{ 0, 1, 11, 2 },
	{ 0, 1, 11, 3 },
	{ 0, 4, 0, 0 },
	{ 0, 4, 0, 1 },
	{ 0, 4, 0, 2 },
	{ 0, 4, 0, 3 },
	{ 0, 4, 0, 4 },
	{ 0, 4, 0, 5 },
	{ 0, 4, 0, 6 },
	{ 0, 4, 0, 7 },
	{ 0, 4, 1, 0 },
	{ 0, 4, 1, 1 },
	{ 0, 4, 1, 2 },
	{ 0, 4, 1, 3 },
	{ 0, 4, 1, 4 },
	{ 0, 4, 1, 5 },
	{ 0, 4, 1, 6 },
	{ 0, 4, 1, 7 },
	{ 0, 4, 8, 0 },
	{ 0, 4, 8, 1 },
	{ 0, 4, 9, 0 },
	{ 0, 4, 9, 1 },
	{ 0, 5, 0, 0 },
	{ 0, 5, 0, 1 },
	{ 0, 5, 0, 2 },
	{ 0, 5, 0, 3 },
	{ 0, 5, 0, 4 },
	{ 0, 5, 0, 5 },
	{ 0, 5, 0, 6 },
	{ 0, 5, 0, 7 },
	{ 0, 5, 1, 0 },
	{ 0, 5, 1, 1 },
	{ 0, 5, 1, 2 },
	{ 0, 5, 1, 3 },
	{ 0, 5, 1, 4 },
	{ 0, 5, 1, 5 },
	{ 0, 5, 1, 6 },
	{ 0, 5, 1, 7 },
	{ 0, 5, 8, 0 },
	{ 0, 5, 8, 1 },
	{ 0, 5, 9, 0 },
	{ 0, 5, 9, 1 },
	{ 0, 8, 0, 0 },
	{ 0, 8, 0, 1 },
	{ 0, 8, 0, 4 },
	{ 0, 8, 0, 5 },
	{ 0, 8, 1, 0 },
	{ 0, 8, 1, 1 },
	{ 0, 8, 1, 4 },
	{ 0, 8, 1, 5 },
	{ 0, 8, 2, 0 },
	{ 0, 8, 2, 1 },
	{ 0, 8, 3, 0 },
	{ 0, 8, 3, 1 },
	{ 0, 8, 8, 0 },
	{ 0, 8, 8, 1 },
	{ 0, 8, 9, 0 },
	{ 0, 8, 9, 1 },
	{ 0, 8, 10, 0 },
	{ 0, 8, 10, 1 },
	{ 0, 8, 11, 0 },
	{ 0, 8, 11, 1 },
	{ 0, 9, 0, 0 },
	{ 0, 9, 0, 1 },
	{ 0, 9, 0, 4 },
	{ 0, 9, 0, 5 },
	{ 0, 9, 1, 0 },
	{ 0, 9, 1, 1 },
	{ 0, 9, 1, 4 },
	{ 0, 9, 1, 5 },
	{ 0, 9, 2, 0 },
	{ 0, 9, 2, 1 },
	{ 0, 9, 3, 0 },
	{ 0, 9, 3, 1 },
	{ 0, 9, 8, 0 },
	{ 0, 9, 8, 1 },
	{ 0, 9, 9, 0 },
	{ 0, 9, 9, 1 },
	{ 0, 9, 10, 0 },
	{ 0, 9, 10, 1 },
	{ 0, 9, 11, 0 },
	{ 0, 9, 11, 1 },
	{ 0, 12, 0, 0 },
	{ 0, 12, 0, 1 },
	{ 0, 12, 0, 4 },
	{ 0, 12, 0, 5 },
	{ 0, 12, 1, 0 },
	{ 0, 12, 1, 1 },
	{ 0, 12, 1, 4 },
	{ 0, 12, 1, 5 },
	{ 0, 12, 8, 0 },
	{ 0, 12, 8, 1 },
	{ 0, 12, 9, 0 },
	{ 0, 12, 9, 1 },
	{ 0, 13, 0, 0 },
	{ 0, 13, 0, 1 },
	{ 0, 13, 0, 4 },
	{ 0, 13, 0, 5 },
	{ 0, 13, 1, 0 },
	{ 0, 13, 1, 1 },
	{ 0, 13, 1, 4 },
	{ 0, 13, 1, 5 },
	{ 0, 13, 8, 0 },
	{ 0, 13, 8, 1 },
	{ 0, 13, 9, 0 },
	{ 0, 13, 9, 1 },
	{ 2, 0, 0, 0 },
	{ 2, 0, 0, 1 },
	{ 2, 0, 0, 2 },
	{ 2, 0, 0, 3 },
	{ 2, 0, 0, 4 },
	{ 2, 0, 0, 5 },
	{ 2, 0, 0, 6 },
	{ 2, 0, 0, 7 },
	{ 2, 0, 1, 0 },
	{ 2, 0, 1, 1 },
	{ 2, 0, 1, 2 },
	{ 2, 0, 1, 3 },
	{ 2, 0, 1, 4 },
	{ 2, 0, 1, 5 },
	{ 2, 0, 1, 6 },
	{ 2, 0, 1, 7 },
	{ 2, 0, 2, 0 },
	{ 2, 0, 2, 1 },
	{ 2, 0, 2, 2 },
	{ 2, 0, 2, 3 },
	{ 2, 0, 2, 4 },
	{ 2, 0, 2, 5 },
	{ 2, 0, 2, 6 },
	{ 2, 0, 2, 7 },
	{ 2, 0, 3, 0 },
	{ 2, 0, 3, 1 },
	{ 2, 0, 3, 2 },
	{ 2, 0, 3, 3 },
	{ 2, 0, 3, 4 },
	{ 2, 0, 3, 5 },
	{ 2, 0, 3, 6 },
	{ 2, 0, 3, 7 },
	{ 2, 0, 8, 0 },
	{ 2, 0, 8, 1 },
	{ 2, 0, 8, 2 },
	{ 2, 0, 8, 3 },
	{ 2, 0, 9, 0 },
	{ 2, 0, 9, 1 },
	{ 2, 0, 9, 2 },
	{ 2, 0, 9, 3 },
	{ 2, 0, 10, 0 },
	{ 2, 0, 10, 1 },
	{ 2, 0, 10, 2 },
	{ 2, 0, 10, 3 },
	{ 2, 0, 11, 0 },
	{ 2, 0, 11, 1 },
	{ 2, 0, 11, 2 },
	{ 2, 0, 11, 3 },
	{ 2, 4, 0, 0 },
	{ 2, 4, 0, 1 },
	{ 2, 4, 0, 2 },
	{ 2, 4, 0, 3 },
	{ 2, 4, 0, 4 },
	{ 2, 4, 0, 5 },
	{ 2, 4, 0, 6 },
	{ 2, 4, 0, 7 },
	{ 2, 4, 8, 0 },
	{ 2, 8, 0, 0 },
	{ 2, 8, 0, 4 },
	{ 2, 8, 1, 0 },
	{ 2, 8, 2, 0 },
	{ 2, 8, 3, 0 },
	{ 2, 8, 8, 0 },
	{ 2, 8, 9, 0 },
	{ 2, 8, 10, 0 },
	{ 2, 8, 11, 0 },
	{ 2, 12, 0, 0 },
	{ 2, 12, 0, 4 },
	{ 2, 12, 8, 0 },
	{ 4, 0, 0, 0 },
	{ 4, 0, 0, 1 },
	{ 4, 0, 0, 2 },
	{ 4, 0, 0, 3 },
	{ 4, 0, 0, 4 },
	{ 4, 0, 0, 5 },
	{ 4, 0, 0, 6 },
	{ 4, 0, 0, 7 },
	{ 4, 0, 2, 0 },
	{ 4, 0, 2, 1 },
	{ 4, 0, 2, 2 },
	{ 4, 0, 2, 3 },
	{ 4, 0, 2, 4 },
	{ 4, 0, 2, 5 },
	{ 4, 0, 2, 6 },
	{ 4, 0, 2, 7 },
	{ 4, 0, 8, 0 },
	{ 4, 0, 8, 2 },
	{ 4, 0, 10, 0 },
	{ 4, 0, 10, 2 },
	{ 4, 1, 0, 0 },
	{ 4, 1, 0, 1 },
	{ 4, 1, 0, 2 },
	{ 4, 1, 0, 3 },
	{ 4, 1, 0, 4 },
	{ 4, 1, 0, 5 },
	{ 4, 1, 0, 6 },
	{ 4, 1, 0, 7 },
	{ 4, 1, 8, 0 },
	{ 4, 4, 0, 0 },
	{ 4, 4, 0, 1 },
	{ 4, 4, 0, 2 },
	{ 4, 4, 0, 3 },
	{ 4, 4, 0, 4 },
	{ 4, 4, 0, 5 },
	{ 4, 4, 0, 6 },
	{ 4, 4, 0, 7 },
	{ 4, 4, 8, 0 },
	{ 4, 5, 0, 0 },
	{ 4, 5, 0, 1 },
	{ 4, 5, 0, 2 },
	{ 4, 5, 0, 3 },
	{ 4, 5, 0, 4 },
	{ 4, 5, 0, 5 },
	{ 4, 5, 0, 6 },
	{ 4, 5, 0, 7 },
	{ 4, 5, 8, 0 },
	{ 4, 8, 0, 0 },
	{ 4, 8, 0, 1 },
	{ 4, 8, 0, 4 },
	{ 4, 8, 0, 5 },
	{ 4, 8, 2, 0 },
	{ 4, 8, 8, 0 },
	{ 4, 8, 10, 0 },
	{ 4, 9, 0, 0 },
	{ 4, 9, 0, 1 },
	{ 4, 9, 0, 4 },
	{ 4, 9, 0, 5 },
	{ 4, 9, 8, 0 },
	{ 4, 12, 0, 0 },
	{ 4, 12, 0, 1 },
	{ 4, 12, 0, 4 },
	{ 4, 12, 0, 5 },
	{ 4, 12, 8, 0 },
	{ 4, 13, 0, 0 },
	{ 4, 13, 0, 1 },
	{ 4, 13, 0, 4 },
	{ 4, 13, 0, 5 },
	{ 4, 13, 8, 0 },
	{ 6, 0, 0, 0 },
	{ 6, 0, 0, 1 },
	{ 6, 0, 0, 2 },
	{ 6, 0, 0, 3 },
	{ 6, 0, 0, 4 },
	{ 6, 0, 0, 5 },
	{ 6, 0, 0, 6 },
	{ 6, 0, 0, 7 },
	{ 6, 0, 2, 0 },
	{ 6, 0, 2, 1 },
	{ 6, 0, 2, 2 },
	{ 6, 0, 2, 3 },
	{ 6, 0, 2, 4 },
	{ 6, 0, 2, 5 },
	{ 6, 0, 2, 6 },
	{ 6, 0, 2, 7 },
	{ 6, 0, 8, 0 },
	{ 6, 0, 8, 2 },
	{ 6, 0, 10, 0 },
	{ 6, 0, 10, 2 },
	{ 6, 4, 0, 0 },
	{ 6, 4, 0, 1 },
	{ 6, 4, 0, 2 },
	{ 6, 4, 0, 3 },
	{ 6, 4, 0, 4 },
	{ 6, 4, 0, 5 },
	{ 6, 4, 0, 6 },
	{ 6, 4, 0, 7 },
	{ 6, 4, 8, 0 },
	{ 6, 8, 0, 0 },
	{ 6, 8, 0, 4 },
	{ 6, 8, 2, 0 },
	{ 6, 8, 8, 0 },
	{ 6, 8, 10, 0 },
	{ 6, 12, 0, 0 },
	{ 6, 12, 0, 4 },
	{ 6, 12, 8, 0 },
	{ 8, 0, 0, 0 },
	{ 8, 0, 0, 2 },
	{ 8, 0, 0, 4 },
	{ 8, 0, 0, 6 },
	{ 8, 0, 1, 0 },
	{ 8, 0, 1, 2 },
	{ 8, 0, 2, 0 },
	{ 8, 0, 2, 2 },
	{ 8, 0, 2, 4 },
	{ 8, 0, 2, 6 },
	{ 8, 0, 3, 0 },
	{ 8, 0, 3, 2 },
	{ 8, 0, 8, 0 },
	{ 8, 0, 8, 2 },
	{ 8, 0, 9, 0 },
	{ 8, 0, 9, 2 },
	{ 8, 0, 10, 0 },
	{ 8, 0, 10, 2 },
	{ 8, 0, 11, 0 },
	{ 8, 0, 11, 2 },
	{ 8, 1, 0, 0 },
	{ 8, 1, 0, 4 },
	{ 8, 1, 1, 0 },
	{ 8, 1, 2, 0 },
	{ 8, 1, 3, 0 },
	{ 8, 1, 8, 0 },
	{ 8, 1, 9, 0 },
	{ 8, 1, 10, 0 },
	{ 8, 1, 11, 0 },
	{ 8, 4, 0, 0 },
	{ 8, 4, 0, 2 },
	{ 8, 4, 0, 4 },
	{ 8, 4, 0, 6 },
	{ 8, 4, 1, 0 },
	{ 8, 4, 8, 0 },
	{ 8, 4, 9, 0 },
	{ 8, 5, 0, 0 },
	{ 8, 5, 0, 4 },
	{ 8, 5, 1, 0 },
	{ 8, 5, 8, 0 },
	{ 8, 5, 9, 0 },
	{ 8, 8, 0, 0 },
	{ 8, 8, 0, 4 },
	{ 8, 8, 1, 0 },
	{ 8, 8, 2, 0 },
	{ 8, 8, 3, 0 },
	{ 8, 8, 8, 0 },
	{ 8, 8, 9, 0 },
	{ 8, 8, 10, 0 },
	{ 8, 8, 11, 0 },
	{ 8, 9, 0, 0 },
	{ 8, 9, 0, 4 },
	{ 8, 9, 1, 0 },
	{ 8, 9, 2, 0 },
	{ 8, 9, 3, 0 },
	{ 8, 9, 8, 0 },
	{ 8, 9, 9, 0 },
	{ 8, 9, 10, 0 },
	{ 8, 9, 11, 0 },
	{ 8, 12, 0, 0 },
	{ 8, 12, 0, 4 },
	{ 8, 12, 1, 0 },
	{ 8, 12, 8, 0 },
	{ 8, 12, 9, 0 },
	{ 8, 13, 0, 0 },
	{ 8, 13, 0, 4 },
	{ 8, 13, 1, 0 },
	{ 8, 13, 8, 0 },
	{ 8, 13, 9, 0 },
	{ 10, 0, 0, 0 },
	{ 10, 0, 0, 2 },
	{ 10, 0, 0, 4 },
	{ 10, 0, 0, 6 },
	{ 10, 0, 1, 0 },
	{ 10, 0, 1, 2 },
	{ 10, 0, 2, 0 },
	{ 10, 0, 2, 2 },
	{ 10, 0, 2, 4 },
	{ 10, 0, 2, 6 },
	{ 10, 0, 3, 0 },
	{ 10, 0, 3, 2 },
	{ 10, 0, 8, 0 },
	{ 10, 0, 8, 2 },
	{ 10, 0, 9, 0 },
	{ 10, 0, 9, 2 },
	{ 10, 0, 10, 0 },
	{ 10, 0, 10, 2 },
	{ 10, 0, 11, 0 },
	{ 10, 0, 11, 2 },
	{ 10, 4, 0, 0 },
	{ 10, 4, 0, 2 },
	{ 10, 4, 0, 4 },
	{ 10, 4, 0, 6 },
	{ 10, 4, 8, 0 },
	{ 10, 8, 0, 0 },
	{ 10, 8, 0, 4 },
	{ 10, 8, 1, 0 },
	{ 10, 8, 2, 0 },
	{ 10, 8, 3, 0 },
	{ 10, 8, 8, 0 },
	{ 10, 8, 9, 0 },
	{ 10, 8, 10, 0 },
	{ 10, 8, 11, 0 },
	{ 10, 12, 0, 0 },
	{ 10, 12, 0, 4 },
	{ 10, 12, 8, 0 },
	{ 12, 0, 0, 0 },
	{ 12, 0, 0, 2 },
	{ 12, 0, 0, 4 },
	{ 12, 0, 0, 6 },
	{ 12, 0, 2, 0 },
	{ 12, 0, 2, 2 },
	{ 12, 0, 2, 4 },
	{ 12, 0, 2, 6 },
	{ 12, 0, 8, 0 },
	{ 12, 0, 8, 2 },
	{ 12, 0, 10, 0 },
	{ 12, 0, 10, 2 },
	{ 12, 1, 0, 0 },
	{ 12, 1, 0, 4 },
	{ 12, 1, 8, 0 },
	{ 12, 4, 0, 0 },
	{ 12, 4, 0, 2 },
	{ 12, 4, 0, 4 },
	{ 12, 4, 0, 6 },
	{ 12, 4, 8, 0 },
	{ 12, 5, 0, 0 },
	{ 12, 5, 0, 4 },
	{ 12, 5, 8, 0 },
	{ 12, 8, 0, 0 },
	{ 12, 8, 0, 4 },
	{ 12, 8, 2, 0 },
	{ 12, 8, 8, 0 },
	{ 12, 8, 10, 0 },
	{ 12, 9, 0, 0 },
	{ 12, 9, 0, 4 },
	{ 12, 9, 8, 0 },
	{ 12, 12, 0, 0 },
	{ 12, 12, 0, 4 },
	{ 12, 12, 8, 0 },
	{ 12, 13, 0, 0 },
	{ 12, 13, 0, 4 },
	{ 12, 13, 8, 0 },
	{ 14, 0, 0, 0 },
	{ 14, 0, 0, 2 },
	{ 14, 0, 0, 4 },
	{ 14, 0, 0, 6 },
	{ 14, 0, 2, 0 },
	{ 14, 0, 2, 2 },
	{ 14, 0, 2, 4 },
	{ 14, 0, 2, 6 },
	{ 14, 0, 8, 0 },
	{ 14, 0, 8, 2 },
	{ 14, 0, 10, 0 },
	{ 14, 0, 10, 2 },
	{ 14, 4, 0, 0 },
	{ 14, 4, 0, 2 },
	{ 14, 4, 0, 4 },
	{ 14, 4, 0, 6 },
	{ 14, 4, 8, 0 },
	{ 14, 8, 0, 0 },
	{ 14, 8, 0, 4 },
	{ 14, 8, 2, 0 },
	{ 14, 8, 8, 0 },
	{ 14, 8, 10, 0 },
	{ 14, 12, 0, 0 },
	{ 14, 12, 0, 4 },
	{ 14, 12, 8, 0 }};

#endif
