//
// 5 Aug 2026
//
// Class RangeSums implements data structures to support the following queries: 
// - sum	: returns a sum of weights over a constrained collection of sets
// - rnd	: returns a random set drawn from the given distribution
// 

#ifndef RangeSums_HPP
#define RangeSums_HPP

#include "basics.hpp"
#include "BD.hpp"
#include "Breal.hpp"


//#define INTERSECTS(A, B)	(A & B)
#define INTERSECTS(A, B)	((A & B).any())
#define SUBSETEQ(A, B)		(A == (A & B))
#define EMPTY_OR_INTERSECTS(A, B)	(A.none() || (A & B).any())

#define Treal	B2real						// RangeSums operate with this type of objects, only internally.

template<class T>
struct seTreal { T set; Treal weight; };

using w64	= seTreal<b64>;
using w256	= seTreal<b256>;
using w1024	= seTreal<b1024>;

template<class T>
using vsw = std::vector<seTreal<T>>;

template<typename T>
bool wcmp	(T a, T b)	{ return a.weight > b.weight; }

class RangeSums {
	public:
		int		id;
		int		nsorted = 0;	// The number of first elements sorted in decreasing order by weight.

		void	init	(int n0)	{ n = n0; wsco64.clear(); wsco256.clear(); wsco1024.clear(); }
		int		size	()			{ return wsco64.size() + wsco256.size() + wsco1024.size(); }
		void	print	();														// Prints out the score list; currently just weights.
		void 	pjkl	(ofstream& f, vint& can);
		
		void	insert	(const vint& set, double logval);
		void	srt		();														// Sorts the list in decreasing order by weight.
		
		double	weight	(int t);												// The t-th weight.
		vint	set		(int t);												// The t-th set.
		
		double	sum		(const vint& L, const vint& U);							// Log of sum of weights.
		double	sum		(const vint& L, const vint& U, double bsum);			// Log of sum of weights, added to exp(bsum).		
		vint	cum		(const vint& L, const vint& U, double csum);			// The set at quantile exp(csum).		
		vint	cum		(const vint& L, const vint& U, double csum, double bsum);// The set at quantile exp(csum) - exp(bsum).				
		vint	rnd		(const vint& L, const vint& U);							// A random set.

		void	psum	(int v, int k, const vvin& D, vdou& Zv);
		vint	pcum	(int v, int k, const vvin& D, u32 Kv, double csum);		
		vint	pcum	(int v, int k, const vvin& D, u32 Kv, double csum, double bsum);

		void	set_tolerance(double d){ delta = d; }		

	private:
		template<class T> double scansum	(const vsw<T>& s, const vint& L, const vint& U);
		template<class T> double scansum	(const vsw<T>& s, const vint& L, const vint& U, double bsum);
		template<class T> vint   scancum	(const vsw<T>& s, const vint& L, const vint& U, double csum);
		template<class T> vint   scancum	(const vsw<T>& s, const vint& L, const vint& U, double csum, double bsum);
		
		template<class T> void	 scanpsum	(const vsw<T>& s, int v, int k, const vvin& D, vdou& Zv);
		template<class T> vint   scanpcum	(const vsw<T>& s, int v, int k, const vvin& D, u32 Kv, double csum);
		template<class T> vint   scanpcum	(const vsw<T>& s, int v, int k, const vvin& D, u32 Kv, double csum, double bsum);

		int			n;				// Size of the ground set.
		vsw<b64>	wsco64;			// Set operations by simple bitwise operations. Use if n <= 64.	
		vsw<b256>	wsco256;		// Set operations by simple bitwise operations. Use if 64 < n <= 256.	
		vsw<b1024>	wsco1024;		// Set operations by a seuqence of bitwise operations. Use if n > 256.		
		double		delta = 0.00001;	// Tolerance, maximum relative error for sum queries.
};

/////////////////////////////
// RangeSums-implementations:

void RangeSums::print(){
	using namespace std;
	if (n <= 64){
		for (auto x : wsco64) cerr << " " << x.weight.get_log() << " " << std::bitset<27>(x.set.to_ullong()) << endl;
		cerr << endl;
	}
	else if (n <= 1024){
		for (auto x : wsco1024) cerr << " " << x.weight;
		cerr << endl;	
	} 
}
void RangeSums::pjkl(ofstream &f, vint &candi){
	using namespace std;
	//f << id << " " << size() << endl; // Don't print the size here, but in the calling function. 
	if (n <= 64){
		for (auto &x : wsco64){ 
			f << FIXED_F(x.weight.get_log(), 12, 4) << " "; 
			vint v = make_vint(x.set, n);
			f << v.size(); for (auto y : v) f << " " << candi[y];
			f << endl;
		}
	}
	else if (n <= 256){
		for (auto &x : wsco256){ 
			f << FIXED_F(x.weight.get_log(), 12, 4) << " "; 
			vint v = make_vint(x.set, n);
			f << v.size(); for (auto y : v) f << " " << candi[y];
			f << endl;
		}
	}	
	else if (n <= 1024){
		for (auto &x : wsco1024){ 
			f << FIXED_F(x.weight.get_log(), 12, 4) << " "; 
			vint v = make_vint(x.set, n);
			f << v.size(); for (auto y : v) f << " " << candi[y];
			f << endl;
		}
	} 
}

void RangeSums::insert(const vint& set, double logval){
	if (std::isinf(logval)) return; // Do not insert zeros.
	if (n <= 64){
		b64 b = make_bitset<b64>(set);
		Treal w; w.set_log(logval);
		wsco64.push_back({ b, w });
	}
	else if (n <= 256){
		b256 b = make_bitset<b256>(set);
		Treal w; w.set_log(logval);
		wsco256.push_back({ b, w });
	}
	else if (n <= 1024){
		b1024 b = make_bitset<b1024>(set);
		Treal w; w.set_log(logval);
		wsco1024.push_back({ b, w });	
	}
	else {
		std::cerr << "\t ERROR: RangeSums::insert not implemented for n > 1024; n = " << n << "\n";
	}
}

void RangeSums::srt(){
	if 		(n <=   64){ std::sort(  wsco64.begin(),   wsco64.end(), wcmp<  w64>); nsorted =   wsco64.size(); }
	else if (n <=  256){ std::sort( wsco256.begin(),  wsco256.end(), wcmp< w256>); nsorted =  wsco256.size(); }
	else if (n <= 1024){ std::sort(wsco1024.begin(), wsco1024.end(), wcmp<w1024>); nsorted = wsco1024.size(); }
	else { std::cerr << "\t ERROR: RangeSums::srt not implemented for n > 1024; n = " << n << "\n"; exit(1); }
	//cerr << " id = " << id << "  sorted; nsorted = " << nsorted << endl;
}

double RangeSums::weight(int t){
	Treal w;
	if		(n <=   64)	w = wsco64[t].weight;
	else if (n <=  256)	w = wsco256[t].weight; 
	else if (n <= 1024)	w = wsco1024[t].weight; 
	else { std::cerr << "\t ERROR: RangeSums::sum not implemented for n > 1024; n = " << n << "\n"; exit(1); }
	return w.get_log();	
}

vint RangeSums::set(int t){
	if		(n <=   64)	return make_vint<  b64>(  wsco64[t].set, n);
	else if (n <=  256)	return make_vint< b256>( wsco256[t].set, n);
	else if (n <= 1024)	return make_vint<b1024>(wsco1024[t].set, n);
	else { std::cerr << "\t ERROR: RangeSums::sum not implemented for n > 1024; n = " << n << "\n"; exit(1); }	
}

double RangeSums::sum(const vint& L, const vint& U){
	if		(n <=   64)	return scansum<  b64>(  wsco64, L, U);
	else if (n <=  256)	return scansum< b256>( wsco256, L, U); 
	else if (n <= 1024)	return scansum<b1024>(wsco1024, L, U); 
	else { std::cerr << "\t ERROR: RangeSums::sum not implemented for n > 1024; n = " << n << "\n"; exit(1); }
	return 0;
}
double RangeSums::sum(const vint& L, const vint& U, double bsum){
	if		(n <=   64)	return scansum<  b64>(  wsco64, L, U, bsum);
	else if (n <=  256)	return scansum< b256>( wsco256, L, U, bsum); 
	else if (n <= 1024)	return scansum<b1024>(wsco1024, L, U, bsum); 
	else { std::cerr << "\t ERROR: RangeSums::sum not implemented for n > 1024; n = " << n << "\n"; exit(1); }
	return 0;
}
vint RangeSums::cum(const vint& L, const vint& U, double csum){
	if		(n <=   64)	return scancum<  b64>(  wsco64, L, U, csum);
	else if (n <=  256)	return scancum< b256>( wsco256, L, U, csum); 
	else if (n <= 1024)	return scancum<b1024>(wsco1024, L, U, csum); 
	else { std::cerr << "\t ERROR: RangeSums::sum not implemented for n > 1024; n = " << n << "\n"; exit(1); }
}
vint RangeSums::cum(const vint& L, const vint& U, double csum, double bsum){
	if		(n <=   64)	return scancum<  b64>(  wsco64, L, U, csum, bsum);
	else if (n <=  256)	return scancum< b256>( wsco256, L, U, csum, bsum); 
	else if (n <= 1024)	return scancum<b1024>(wsco1024, L, U, csum, bsum); 
	else { std::cerr << "\t ERROR: RangeSums::sum not implemented for n > 1024; n = " << n << "\n"; exit(1); }
}
vint RangeSums::rnd(const vint& L, const vint& U){
	return cum(L, U, sum(L, U) + log(urnd()) );
}

void RangeSums::psum(int v, int k, const vvin& D, vdou& Zv){
	if		(n <=   64)	scanpsum<  b64>(  wsco64, v, k, D, Zv);
	else if (n <=  256)	scanpsum< b256>( wsco256, v, k, D, Zv); 
	else if (n <= 1024)	scanpsum<b1024>(wsco1024, v, k, D, Zv); 
	else { std::cerr << "\t ERROR: RangeSums::psum not implemented for n > 1024; n = " << n << "\n"; exit(1); }
}
vint RangeSums::pcum(int v, int k, const vvin& D, u32 Kv, double csum){
	if		(n <=   64)	return scanpcum<  b64>(  wsco64, v, k, D, Kv, csum);
	else if (n <=  256)	return scanpcum< b256>( wsco256, v, k, D, Kv, csum); 
	else if (n <= 1024)	return scanpcum<b1024>(wsco1024, v, k, D, Kv, csum); 
	else { std::cerr << "\t ERROR: RangeSums::pcum not implemented for n > 1024; n = " << n << "\n"; exit(1); }
}		
vint RangeSums::pcum(int v, int k, const vvin& D, u32 Kv, double csum, double bsum){
	if		(n <=   64)	return scanpcum<  b64>(  wsco64, v, k, D, Kv, csum, bsum);
	else if (n <=  256)	return scanpcum< b256>( wsco256, v, k, D, Kv, csum, bsum); 
	else if (n <= 1024)	return scanpcum<b1024>(wsco1024, v, k, D, Kv, csum, bsum); 
	else { std::cerr << "\t ERROR: RangeSums::pcum not implemented for n > 1024; n = " << n << "\n"; exit(1); }
}


template<class T> double RangeSums::scansum(const vsw<T>& s, const vint& Lv, const vint& Uv){
	T L = make_bitset<T>(Lv), U = make_bitset<T>(Uv);
	if ( !EMPTY_OR_INTERSECTS(L, U) ) return -INFINITY; // No contributions from the list.
	int m = nsorted;
	Treal sum, factor, slack; slack = delta / (m + 1.0);
	int i = 0, count = 0;
	for (; i < m; ++i){
		auto P = s[i].set;
		if ( EMPTY_OR_INTERSECTS(L, P) && SUBSETEQ(P, U) ){ 
			sum = s[i].weight; factor = sum * slack; 
			++i; ++count; break;
		}	
	}
	for (; i < m; ++i){
		auto P = s[i].set;
		if ( EMPTY_OR_INTERSECTS(L, P) && SUBSETEQ(P, U) ){ 
			Treal w = s[i].weight;
			if (w < factor) break; 
			sum += w; ++count;
		}	
	}
	return sum.get_log();
}
template<class T> double RangeSums::scansum(const vsw<T>& s, const vint& Lv, const vint& Uv, double bsum){
	T L = make_bitset<T>(Lv), U = make_bitset<T>(Uv);
	if ( !EMPTY_OR_INTERSECTS(L, U) ) return bsum; // No contributions from the list.
	int m = nsorted;
	Treal sum; sum.set_log(bsum); Treal factor, slack; slack = delta / (m + 1.0);
	int i = 0, count = 0;
	for (; i < m; ++i){
		auto P = s[i].set;
		if ( EMPTY_OR_INTERSECTS(L, P) && SUBSETEQ(P, U) ){ 
			sum += s[i].weight; factor = sum * slack; 
			++i; ++count; break;
		}	
	}
	for (; i < m; ++i){
		auto P = s[i].set;
		if ( EMPTY_OR_INTERSECTS(L, P) && SUBSETEQ(P, U) ){ 
			Treal w = s[i].weight;
			if (w < factor) break; 
			sum += w; ++count;
		}	
	}
	return sum.get_log();
}
template<class T> vint RangeSums::scancum(const vsw<T>& s, const vint& Lv, const vint& Uv, double csum){
	T L = make_bitset<T>(Lv), U = make_bitset<T>(Uv), P, X;
	int m = nsorted;
	Treal sum, target; target.set_log(csum);	
	for (int i = 0; i < m; ++i){
		P = s[i].set;
		if ( EMPTY_OR_INTERSECTS(L, P) && SUBSETEQ(P, U) ){ 
			sum += s[i].weight;
			if (sum > target){ X = P; break; }
		}	
	}
	return make_vint<T>(X, n);
}
template<class T> vint RangeSums::scancum(const vsw<T>& s, const vint& Lv, const vint& Uv, double csum, double bsum){
	T L = make_bitset<T>(Lv), U = make_bitset<T>(Uv), P, X;
	int m = nsorted;
	Treal sum, target; sum.set_log(bsum); target.set_log(csum);
	for (int i = 0; i < m; ++i){
		P = s[i].set;
		if ( EMPTY_OR_INTERSECTS(L, P) && SUBSETEQ(P, U) ){ 
			sum += s[i].weight;
			if (sum > target){ X = P; break; }
		}	
	}
	return make_vint<T>(X, n);
}

template<class T> void RangeSums::scanpsum(const vsw<T>& s, int v, int k, const vvin& D, vdou& Zv){
	// Transform D[.] to type T.
	vector<T> tD(k); for (int u = 0; u < k; u++) tD[u] = make_bitset<T>(D[u]);
	
	vint R(k-1); for (int u = 0; u < k-1; u++) R[u] = (u < v) ? u : u + 1;

	// From Zv to Z.
	vector<Treal> Z(1 << k); for (u32 K = 0; K < (1 << k); K++) Z[K].set_log(Zv[K]);

	int m = nsorted;
	for (int i = 0; i < m; ++i){
		// Obtain set J.
		T J = s[i].set;
		if ( INTERSECTS(J, tD[v]) ) continue;
	
		// Obtain set K.
		u32 K = 0;	
		for (auto u : R){
			K |= ( INTERSECTS(J, tD[u]) ) ? (1 << u) : 0;					
		}					
		Z[K] += s[i].weight;		
	}
	// From Z to Zv.
	for (u32 K = 0; K < (1 << k); K++) Zv[K] = Z[K].get_log();
}

template<class T> vint RangeSums::scanpcum(const vsw<T>& s, int v, int k, const vvin& D, u32 Kv, double csum){
	// Transform D[.] to type T.
	vector<T> tD(k); for (int u = 0; u < k; u++) tD[u] = make_bitset<T>(D[u]);

	vint R(k-1); for (int u = 0; u < k-1; u++) R[u] = (u < v) ? u : u + 1;

	// From Zv to Z.
	Treal ZK; //ZK.set_log(bsum);
	Treal target; target.set_log(csum); 
	int iJ = -1;

	int m = nsorted;
	for (int i = 0; i < m; ++i){
		// Obtain set J.
		auto J = s[i].set;
		if ( INTERSECTS(J, tD[v]) ) continue;

		// Obtain set K.
		u32 K = 0;	
		for (auto u : R){
			K |= ( INTERSECTS(J, tD[u]) ) ? (1 << u) : 0;					
		}					
		if (K == Kv){
			ZK += s[i].weight;
			if (ZK >= target){ iJ = i; break; }	
		}
	}
	return make_vint(s[iJ].set, n);
}

template<class T> vint RangeSums::scanpcum(const vsw<T>& s, int v, int k, const vvin& D, u32 Kv, double csum, double bsum){
	// Transform D[.] to type T.
	vector<T> tD(k); for (int u = 0; u < k; u++) tD[u] = make_bitset<T>(D[u]);

	vint R(k-1); for (int u = 0; u < k-1; u++) R[u] = (u < v) ? u : u + 1;

	// From Zv to Z.
	Treal ZK; ZK.set_log(bsum);
	Treal target; target.set_log(csum); 
	int iJ = -1;

	int m = nsorted;
	for (int i = 0; i < m; ++i){
		// Obtain set J.
		auto J = s[i].set;
		if ( INTERSECTS(J, tD[v]) ) continue;

		// Obtain set K.
		u32 K = 0;	
		for (auto u : R){
			K |= ( INTERSECTS(J, tD[u]) ) ? (1 << u) : 0;					
		}					
		if (K == Kv){
			ZK += s[i].weight;
			if (ZK >= target){ iJ = i; break; }	
		}
	}
	if (iJ < 0){ cerr << "Error in scancum. Exit now." << endl; exit(1); } 
	return make_vint(s[iJ].set, n);
}



#undef INTERSECTS
#undef SUBSETEQ
#undef Treal

#endif
