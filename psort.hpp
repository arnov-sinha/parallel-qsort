#ifndef PSORT_H_INCLUDED
#define PSORT_H_INCLUDED

#include <vector>
#include <stdint.h>
#include <cstdlib>
#include <climits>
#include <omp.h>
#include <string>
#include <ctime>
#include <time.h>
#include <string>
#include <algorithm>
#include <unistd.h>

using namespace std;
using std::string;

template <class T>
uint64_t spartition( T &vec, uint64_t pivot, uint64_t lo, int64_t hi );

int8_t requestBlocks( uint64_t &lo, uint64_t &hi, \
		      uint64_t &LB, uint64_t &RB, uint64_t mode, uint64_t _E, \
		      std::vector<uint64_t> &_lesser, \
		      std::vector<uint64_t> &_greater, uint64_t _blockSize );

template<class T>
void insertionSort( T &vec, uint64_t low, uint64_t high );

template<class T>
void basicquicksort( T &vec, uint64_t lo, uint64_t hi );

template<class T>
void group( T &vec, uint64_t &lo, uint64_t &hi, uint64_t _E, \
	    std::vector<uint64_t> &_lesser, std::vector<uint64_t> &_greater,\
	    uint64_t &_LMAX, uint64_t &_GMIN, uint64_t _blockSize );

template<class T>
void cleanup( T &vec, uint64_t _LMAX, uint64_t _GMIN, \
	      std::vector<uint64_t> _lesser, std::vector<uint64_t> _greater,\
	      uint64_t low, uint64_t high, uint64_t _blockSize, uint64_t _E,\
	      uint64_t &p );

template<class T>
void psortOptimized( T &vec, uint64_t lo, uint64_t hi, uint64_t _blockSize );

#endif /* PSORT_H_INCLUDED */
