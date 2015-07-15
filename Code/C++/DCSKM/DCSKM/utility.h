#ifndef __UTILITY
#define __UTILITY

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <math.h>
#include <assert.h>

using namespace std;
#include <vector>
#include <deque>
#include <map>
#include <iostream>
#include <algorithm>
#include <fstream>

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
void printPattern(BitStore* bs);
int getSetIndex(vector<string> &vec,char* str);

void AddAll(FastList<int> &S,FastList<int> &C);
void RemoveAll(FastList<int> &S,FastList<int> &C);

void InitZipfMaxVal(int maxnum,double theta);
int zipf(double theta);
double gaussian(double mean, double sigma);
long poisson(long lambda);
float AvgAbsMeanDev(vector<float> &S);
float variance(vector<float> &S);

//-----------
// LRU buffer
//-----------

enum AccessMode {isData,isIndex};
int getBlockLength();
void InitCache(int csize);
void RefreshCache();
void DestroyCache();
bool getCacheBlock(char* buffer,int UserId,int BlockId);
void storeCacheBlock(char* buffer,int UserId,int BlockId,AccessMode mode);	// user's responsibility
//void printPageAccess();
int64_t printPageAccess();
void RefreshStat();

//----------
// FreqCache
//----------

struct FreqCache
{
    char* buffer;
    int UserId,BlockId;
};

void InitFreqCache(FreqCache& fc);
void storeFreqCache(FreqCache& fc,char* buf,int uid,int bid);
bool inFreqCache(FreqCache& fc,int uid,int bid);
void DestroyFreqCache(FreqCache& fc);

#endif //__UTILITY
