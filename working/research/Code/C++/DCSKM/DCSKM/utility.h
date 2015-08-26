#ifndef __UTILITY
#define __UTILITY
// handle
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <vector>
#include <deque>
#include <map>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <time.h>

using namespace std;


//typedef map<string,string> ConfigType;
#define FastArray	vector
#define FastList	deque
#define BitStore	vector<bool>

#define DEFAULT_CACHESIZE 1024
#define DEFAULT_BLOCKLENGTH	4096


const float FLOAT_MAX=(float)INT_MAX;

#define printIntSet(Obj)\
{\
	for (int __counter=0;__counter<Obj.size();__counter++)\
		printf("%d ",Obj[__counter]);\
	printf("\n");\
}

// "InitClock" also initialize the seeds for 2 random generators
void InitClock();
void PrintElapsed();
void CheckFile(FILE* fp,const char* filename);

//Modified by Qin Xu
//#define min(a, b) (((a) < (b))? (a) : (b)  )
//#define max(a, b) (((a) > (b))? (a) : (b)  )
#define ValAbs(x) (((x) >  0 )? (x) : -(x)  )

//void TrimSpace(char* str);
//void AddConfigFromFile(ConfigType &cr,const char* filename);
//void AddConfigFromCmdLine(ConfigType &cr,int argc,char** argv);
//void ListConfig(ConfigType &cr);
//
//float getConfigFloat(const char* key,ConfigType &cr,bool required=true,float _default=0);
//int getConfigInt(const char* key,ConfigType &cr,bool required=true,int _default=0);
//const char* getConfigStr(const char* key,ConfigType &cr,bool required=true,const char* _default=NULL);

int Cardinality(BitStore* bs);	// cardinality of a set
void printPattern(BitStore* bs); // print the position of bitstore where is ture
int getSetIndex(vector<string> &vec,char* str); // get the index of str, if no then put str into vector

void AddAll(FastList<int> &S,FastList<int> &C); // add all the element of C to S
void RemoveAll(FastList<int> &S,FastList<int> &C); // remove all the elements of C from S

void InitZipfMaxVal(int maxnum,double theta); // compute the zipfMaxValue with maxnum and theta.
int zipf(double theta); // given zipfMaxVal to compute the i between [0, maxnum]; 
double gaussian(double mean, double sigma); // generate the data follow guassian distribution
long poisson(long lambda);
float AvgAbsMeanDev(vector<float> &S); // return ((sum==0)?(1):(diff/sum));???
float variance(vector<float> &S); //??

//-----------
// LRU buffer
//-----------

enum AccessMode {isData,isIndex};
int getBlockLength(); // return the length of block
void InitCache(int csize); // construct a cache space with size csize
void RefreshCache(); // reset the cache 
void DestroyCache(); // free the dynamic parameter
bool getCacheBlock(char* buffer,int UserId,int BlockId); // get the block from cache to buffer
void storeCacheBlock(char* buffer,int UserId,int BlockId,AccessMode mode);	// store the buffer to cache, may invoke the LRU strategy
//void printPageAccess();
int64_t printPageAccess();
void RefreshStat(); // reset the state

//----------
// FreqCache
//----------

struct FreqCache
{
    char* buffer;
    int UserId,BlockId;
};

void InitFreqCache(FreqCache& fc); // initial the freqCache
void storeFreqCache(FreqCache& fc,char* buf,int uid,int bid); // store the buffer to fc
bool inFreqCache(FreqCache& fc,int uid,int bid); // verify whether uid and bid in fc or not
void DestroyFreqCache(FreqCache& fc); // destroy the freqCache

#endif //__UTILITY
