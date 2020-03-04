#include <iostream>
#include "psort.hpp"

using namespace std;
using std::string;

/*
 *Inline support functions
 * pivot: Finds a random pivot point
 */
template<class T>
inline T pivot( T min, T max ){
  return lrand48()%(max-min+1)+min;
}

static inline double compute_elapsed( const struct timespec &starttime )
{
  struct timespec endtime;
  clock_gettime( CLOCK_REALTIME, &endtime );
  double elapsed = (( endtime.tv_sec + endtime.tv_nsec / ( double ) 1000000000)\
                    -( starttime.tv_sec + starttime.tv_nsec / \
		       ( double ) 1000000000 )) ;
  return elapsed;
}

template<class T>
uint64_t spartition(T &vec, uint64_t pivot, uint64_t lo, uint64_t hi)
{
  uint64_t i = lo - 1;
  uint64_t j = hi + 1;
  
  while(true){
    do{
      i++;
    }while(vec[i]<pivot);
    do{
      j--;
    }while(vec[j] > pivot);
    if(i >= j){
      return j;
    }
    std::swap(vec[i],vec[j]);
  }
  return -1;
}

int8_t requestBlocks( uint64_t &lo, uint64_t &hi, \
		      uint64_t &LB, uint64_t &RB, uint64_t mode, uint64_t _E, \
		      vector<uint64_t> &_lesser, vector<uint64_t> &_greater, \
		      uint64_t _blockSize)
{
  int8_t fault = 0;
#ifdef _OPENMP
#pragma omp critical (request_block)
#endif
  {
    switch (mode)
      {
      case 2:
	if(lo + ((4 * _blockSize) - 1) < hi) //This is designed so that if one of the queue is empty, 
	  //the other one is still used for allocation
	  {
	    if(!_greater.empty())
	      {
		RB = _greater.back();
		_greater.pop_back();
	      }
	    else
	      {
		hi -= _blockSize;
		RB = hi - (_blockSize - 1);
	      }
	  
	    if(!_lesser.empty())
	      {
		LB = _lesser.back();
		_lesser.pop_back();
	      }
	    else
	      {
		lo += _blockSize;
		LB = lo;
	      }
	    fault = 1;
	  }
	else if(!_lesser.empty() && !_greater.empty())
	  {
	    RB = _greater.back();
	    _greater.pop_back();
	    LB = _lesser.back();
	    _lesser.pop_back();
	    fault = 1;
	  }
	break;
      
      case 1:
	if( (hi - ((2 * _blockSize) - 1)) > (lo + (_blockSize - 1)) )
	  {
	    if( !_greater.empty() )
	      {
		RB = _greater.back();
		_greater.pop_back();
	      }
	    else
	      {
		hi -= _blockSize;
		RB = hi - (_blockSize - 1);
	      }
	    fault = 1; 
	  }
	else if( !_greater.empty() )
	  {
	    RB = _greater.back();
	    _greater.pop_back();
	    fault = 1;
	  }
            
	if(!fault)
	  {
	    _lesser.push_back(LB);
	  }
	break;
      
      case 0:
	if( (lo + ((2 * _blockSize) - 1)) < (hi - (_blockSize - 1)) )
	  {
	    if( !_lesser.empty() )
	      {
		LB = _lesser.back();
		_lesser.pop_back();
	      }
	    else
	      {
		lo += _blockSize;
		LB = lo;
	      }
	    fault = 1;
	  }
	else if( !_lesser.empty() )
	  {
	    LB = _lesser.back();
	    _lesser.pop_back();
	    fault = 1;
	  }

	if(!fault)
	  {
	    _greater.push_back(RB);
	  }
	break;
      }// end of switch
  }// end of critical
  return fault;
}

inline bool partition_check(std::vector<uint64_t> &vec, uint64_t lo, uint64_t hi, uint64_t pivot, uint64_t _E)
{
  uint64_t i = lo;
  for(; i <= pivot; ++i)
    if(vec[ i ] > _E)
      return false;

  for(; i <= hi; ++i)
    if(vec[ i ] < _E)
      return false;

  return true;
}

inline void savelargeorsmall( uint64_t lo, uint64_t hi, uint64_t smalllimit,\
			      vector<uint64_t> &small, vector<uint64_t> &large )
{
  if( lo < hi )
    {
      uint64_t n = hi - lo + 1 ;
      if( n <= smalllimit )
        {
	  small.push_back( lo ) ;
	  small.push_back( hi ) ;
        }
      else
        {
	  large.push_back( lo ) ;
	  large.push_back( hi ) ;
        }
    }
}

template<class T>
void insertionSort(T &vec, uint64_t low, uint64_t high)
{
  uint64_t i = low + 1;
  while( i <= high )
    {
      uint64_t j = i;
      while( ( j > 0 ) && ( vec[ j ] < vec[ j - 1 ] ) )
	{
	  std::swap( vec[ j ], vec[ j - 1 ] );
	  --j;
	}
      ++i;
    }
}
  
template<class T>
void basicquicksort( T &vec, uint64_t lo, uint64_t hi )
{
  if( lo < hi )
    {
      const uint64_t insertionSize = 8;

      uint64_t currsize = hi - lo + 1;
      if( currsize <= insertionSize )
	insertionSort(vec, lo, hi);
      else
	{
	  uint64_t p = spartition(vec, vec[lo], lo, hi);
	  basicquicksort(vec, lo, p);
	  basicquicksort(vec, p + 1, hi);
	}
    }
}

/*
 *Creating a thread safe version of group function, compatible parallel_optimization function
 *
 *@param
 * vec - source vector
 * lo - lower bound
 * hi - higher bound
 * _blockSize - block size to be used for the partition function
 * _lesser & _greater - queues for saving incomplete work
 * _LMAX - the highest position of the lesser element 
 * _GMIN - the least postion of the Greater element 
 */
template<class T>
void group(T &vec, uint64_t &lo, uint64_t &hi, uint64_t _E, \
	   vector<uint64_t> &_lesser, vector<uint64_t> &_greater, uint64_t &_LMAX, \
	   uint64_t &_GMIN, uint64_t _blockSize)
{
  uint64_t LB;
  uint64_t RF;
  uint64_t LF;
  uint64_t RB;
  uint64_t l;
  uint64_t r;
  int8_t haveWork;

  //Assign blocks to the threads in the beginning
  haveWork = requestBlocks(lo, hi, LB, RB, 2, _E, _lesser, _greater, _blockSize);
 
  LF = LB + (_blockSize - 1);
  RF = RB + (_blockSize - 1);
  l  = LB;
  r  = RB;

  while(haveWork)
    {
    
      for(; l <= LF; l++)
	if(vec[l] > _E)
	  break;

      for(; r <= RF; r++)
    	if(vec[r] < _E)
	  break;

      if(l <= LF && r <= RF)
	{
	  std::swap(vec[l], vec[r]);
	  l++;
	  r++;
	}
      else
	{
	  if( (l > LF) && (r > RF) ) //Requesting both blocks
	    {
	      if(LF > _LMAX)
		_LMAX = LF;

	      if(RB < _GMIN)
		_GMIN = RB;

	      haveWork = requestBlocks(lo, hi, LB, RB, 2, _E, _lesser, _greater,\
				       _blockSize);
	    }
	  else if(r > RF) //Requesting right block
	    {
	      if(RB < _GMIN)
		_GMIN = RB;

	      haveWork = requestBlocks(lo, hi, LB, RB, 1, _E, _lesser, _greater, \
				       _blockSize);
	    }
	  else if(l > LF) //Requesting left block
	    {
	      if(LF > _LMAX)
		_LMAX = LF;

	      haveWork = requestBlocks(lo, hi, LB, RB, 0, _E, _lesser, _greater, \
				       _blockSize);

	    }
	  else{
	    cout<<"*********************************I am dead****************************************"<<endl;
	    haveWork = 0;
	  }

	  LF = LB + (_blockSize - 1);
	  RF = RB + (_blockSize - 1);
	  l  = LB;
	  r  = RB;
  	}//end of request block else
    }//end of forever loop
}

/*
 *Creating a thread safe of cleanup compatible with parallel_optimization function
 *
 *@param
 * vec - source vector
 * low - lower bound
 * high - higher bound
 * _blockSize - block size to be used for the partition function
 * _lesser & _greater - queues for saving incomplete work
 * _LMAX - the highest position of the lesser element 
 * _GMIN - the least postion of the Greater element 
 */
template<class T>
void cleanup(T &vec, uint64_t _LMAX, uint64_t _GMIN, \
	     std::vector<uint64_t> _lesser, std::vector<uint64_t> _greater, uint64_t low,\
	     uint64_t high, uint64_t _blockSize, uint64_t _E, uint64_t &p)
{
  uint64_t partition;
  int64_t lsize;
  int64_t gsize;
  int64_t qsize;
  uint64_t newLow;
  uint64_t newHigh;
  
  lsize = _lesser.size();
  gsize = _greater.size();

  qsize = lsize + gsize;
  switch(qsize)
    {
    case 0: // If there is no elements in the queue left
      if(_LMAX == _GMIN - 1) // if the boundary between the LMAX and GMIN is reached
	{
	  if (_LMAX > _GMIN) std::swap(vec[_LMAX], vec[_GMIN]); 
	  partition = _LMAX;
	}
      else
	{
	  if( _GMIN != high ) _GMIN -= 1;
	  if( _LMAX != low ) std::swap(vec[low], vec[_LMAX]);
	}
      partition = spartition(vec, _E, _LMAX, _GMIN);
      break;

    case 1:
      if(lsize == 1) //if lesser queue only has 1 element left
	{
	  newLow  = _lesser[0];
	  newHigh = _GMIN - 1;

	  if(_blockSize == 1) //Off chance that the blocksize is chosen as 1
	    {
	      (vec[newLow] > _E) ? partition = newLow -1 : partition = newLow;
	    }
	  if(newLow != low)
	    {
	      newLow = newLow - 1;
	      std::swap(vec[low], vec[newLow]);
	    }
	}
      else //if greater queue only has 1 element left
	{
	  newHigh = _greater[0] + (_blockSize - 1);
	  newLow = _LMAX;

	  if(_blockSize == 1) partition = newLow;

	  if(newLow != low)
	    std::swap(vec[low], vec[newLow]);
	}
      
      partition = spartition(vec, _E, newLow, newHigh);
      break;

    default:
      if(lsize > 1) //if lesser queue has more than 1 element left
	{
	  std::sort(_lesser.begin(), _lesser.end());
	  newLow = _lesser.front();
	  newHigh = _GMIN - 1;

	  if(_blockSize == 1) //Off chance that the blocksize is chosen as 1
	    (vec[newLow] > _E) ? partition = newLow -1 : partition = newLow;

	  if(newLow != low)
	    {
	      newLow = newLow - 1;
	      std::swap(vec[low], vec[newLow]);
	    }
	}
      else //if greater queue has more than 1 element left
	{
	  std::sort(_greater.begin(), _greater.end());
	  newLow = _LMAX;
	  newHigh = _greater.back() + (_blockSize - 1);
	  
	  if(_blockSize == 1) partition = _greater.front();

	  if(newLow != low)
	    std::swap(vec[low], vec[newLow]);
	}
      
      partition = spartition(vec, _E, newLow, newHigh);
      break;
    }//switch end
  p = partition;
}

template<class T>
void psortOptimized(T &vec, uint64_t lo, uint64_t hi, \
		    uint64_t _blockSize)
{
  if( lo < hi )
    {
      uint64_t p;
      vector<uint64_t> _lesser;
      vector<uint64_t> _greater;
      vector<uint64_t> smallsegs;
      vector<uint64_t> largesegs;
      const uint64_t dosmallsegs = 512;
      
      smallsegs.reserve( 2 * dosmallsegs );
      largesegs.reserve( 128 *2 );

      uint64_t _EIndex;
      uint64_t _E;     
      uint64_t _LMAX;  
      uint64_t _GMIN;  
      uint64_t m_low;  
      uint64_t m_high;
      
      const uint64_t insertionSize = 8;
      const uint64_t basicQuickSize = 16 * 1024;
      uint64_t ompmaxthreads = 1;

      
#ifdef _OPENMP
      ompmaxthreads = omp_get_max_threads();
#endif

      _lesser.reserve( ompmaxthreads ); // Queue size is the same as the no. of threads
      _greater.reserve( ompmaxthreads );// Queue size is the same as the no. of threads
      
      savelargeorsmall( lo, hi, basicQuickSize, smallsegs, largesegs );

#ifdef _OPENMP
#pragma omp parallel
#endif
      while( ( smallsegs.size() + largesegs.size() ) > 0 )
	{
	  if( ( smallsegs.size() > dosmallsegs ) || ( largesegs.size() == 0 ) )
	    {
	      uint64_t nsmall = smallsegs.size();

#ifdef _OPENMP
	      //#pragma omp for schedule ( guided, 4 )
	      //#pragma omp single
#pragma omp for schedule ( guided , 1 )
#endif
	      for( uint64_t i = 0; i < nsmall; i+=2 )
		{
		  uint64_t  h = smallsegs[ i + 1 ];
		  uint64_t  l = smallsegs[ i + 0 ];

		  if( ( h - l + 1 ) <= insertionSize )
		    insertionSort(vec, l, h);
		  else
		    basicquicksort(vec, l, h);
		}
	      
#pragma omp single
	      {
		smallsegs.clear();
	      }
	    }
	  else
	    {
	      
	      uint64_t  nlarge = largesegs.size();
	      uint64_t  h = largesegs[  nlarge - 1 ];
	      uint64_t  l = largesegs[  nlarge - 2 ];

#pragma omp single
	      {
		_lesser.clear();
		_greater.clear();
		_EIndex        = pivot(l, h);
		_E             = vec[_EIndex];
		_LMAX          = l;
		_GMIN          = h;
		m_low          = l;//this is modified by the group function 
		m_high         = h;//this is modified by the group function
	      
		std::swap(vec[_EIndex], vec[m_low]);
  
		m_low -= _blockSize;
		m_high += _blockSize;
	      }

	      group(vec, m_low, m_high, _E, _lesser, _greater, _LMAX, _GMIN, _blockSize);
#pragma omp barrier
#pragma omp single
	      {
		cleanup(vec, _LMAX, _GMIN, _lesser, _greater, l, h, _blockSize, _E, p);
		// cout<<"Partition_Check = "<<partition_check(vec, l, h, p, _E)<<endl;
		largesegs.resize( nlarge - 2 );	     

		savelargeorsmall( l, p, basicQuickSize, smallsegs, largesegs );
		savelargeorsmall( p + 1, h, basicQuickSize, smallsegs, largesegs );    
	      }// end of single
	    }// end of else
	}// end of while
    }// end of if ( lo < hi )
}

int main(int argc, char** argv)
{ 
  cout<<"------------------Code started-------------------"<<std::endl;
  std::vector<uint64_t> vec;
  std::string temp;
  uint64_t n = 100000000;
  if( argc > 1 )
    n = stoi( argv[1] ) ;
  vec.reserve(n);
 
  uint64_t seed = time(NULL); // Seed the time
  
  srand( seed ) ;
  for( uint64_t i = 0 ; i < n ; ++i )
    //vec.push_back( n - i - 1);
    //vec.emplace_back( i % 10000);
    vec.emplace_back( lrand48() );
  //vec.push_back( i );
    
  
  uint64_t size = vec.size();

  uint64_t _blockSize = 1024;

  srand(seed);

  cout<<"Data Size              = "<<size<<endl;
  cout<<"Block Size             = "<<_blockSize<<endl;
#ifdef _OPENMP
  cout<<"No of threads running  = "<<omp_get_max_threads()<<endl;
#else
  cout<<"No of threads running  = 1"<<endl;
#endif
  cout<<"Seed                   = "<<seed<<endl;

  struct timespec t0 ;
  clock_gettime( CLOCK_REALTIME, &t0 ) ;
  psortOptimized(vec, 0, size , _blockSize);
  cout<<"Time for Parallel Sort = "<<compute_elapsed(t0)<<endl;

  bool check = true;
  for(uint64_t i = 1; i<=size; ++i)
    {
      if(vec[i] < vec[i-1])
	{
	  cout<<"Error Check: FAILED at "<<i<<endl;
	  check = false;
	  break;
	}
    }

  if(check == true)
    cout<<"Error Check: PASSED"<<endl;
  
  vec.clear(); //cleaning the vector

  for( uint64_t i = 0 ; i < n ; ++i )
    //vec.push_back( n - i - 1);
    //vec.emplace_back( i % 10000);
    vec.emplace_back( lrand48() );
  //vec.push_back( i );
  
  cout<<"Running std::sort..."<<endl ;
  clock_gettime( CLOCK_REALTIME, &t0 ) ;
    std::sort(vec.begin(), vec.end());
    cout<<"Time for std::sort: "<<compute_elapsed( t0 )<<endl;
  
  cout<<"-------------------Code ended--------------------"<<std::endl;

  return 0;

}
