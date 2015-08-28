#ifndef __NET_SHARED
#define __NET_SHARED

#include <queue>
#include <bitset>
#include <vector>
#include <set>
using namespace std;
// handle
#define MAX_DIST 9999999
#define ATTRIBUTE_DIMENSION 6
int MAX_KEYWORDS = 64;

int NodeNum;
int EdgeNum;
// define the query point
struct QueryPoint
{
	int Ni,Nj;
    float dist_Ni; // location
    float distCnst; // distrance constraint	
    //unsigned long long keywords; // ?? 和关键字有什么关系，利用二进制记录关键字信息，有的用0，没有用1
    bitset<ATTRIBUTE_DIMENSION> subSpace; // query subspace
	int nOfKwd;
	set<int> kwd; // keyword information
};

// output the query point to ostream
ostream& operator<<(ostream& os,const QueryPoint& Q)
{
    os<<"("<<Q.Ni<<","<<Q.Nj<<")	";
	os<<Q.dist_Ni<<"	"<<Q.distCnst<<"	";
	os<<bitset<ATTRIBUTE_DIMENSION>(Q.subSpace).to_string()<<"	";
	os<<Q.nOfKwd;
	set<int>::iterator it = Q.kwd.begin();
	for (; it != Q.kwd.end(); ) {
		os << "	" << *it;
	}
    return os;
    
}
// define the struct of POI point 
struct POI
{
    //int Ni,Nj;
	int poid;
    //int pre;//Refine use only ?? 
    //unsigned long long keywords;
	float dist_toquery; 
	float attr[ATTRIBUTE_DIMENSION];
    //float dist_Ni;
    //float dist_Nj;
   
};
// redefine the () to compare the distance of query point to POI
struct POIComparison
{
    bool operator () (const POI& left, const POI& right) const
    {
        return left.dist_toquery > right.dist_toquery;
    }
};

/*  ---------------M--not used-------------
// output the POI to ostream
ostream& operator<<(ostream& os, const POI& poi)
{
    os<<poi.Ni<<"	"<<poi.dist_toquery;
	for(int i=0; i<ATTRIBUTE_DIMENSION; i++) {
		os<<"	"<<poi.attr[i];
	}
    //os<<bitset<64>(poi.keywords).to_string();
    //os<<" distNi:"<<poi.dist_Ni;
    //os<<" distQ:"<<poi.dist_toquery;
    //os<<" preNode:"<<poi.pre<<endl;
    
    return os;
}
*/
struct DStepEvent
{
    double dist;
    int node;
};

struct DStepComparison
{
    bool operator () (const DStepEvent& left, const DStepEvent& right) const
    {
        return left.dist > right.dist;
    }
};

typedef	priority_queue<DStepEvent,vector<DStepEvent>,DStepComparison> DStepQueue;

// for Dijkstra
struct edgePair
{
    int Ni;
    int Nj;
};

struct edgeState 
{
	int vState;
	float iDisToQuery;
	float jDisToQuery;
};

struct eSComparison
{
    bool operator () (const edgePair& left, const edgePair& right) const
    {
		//return left.disToQuery > right.disToQuery;
		if(left.Ni == right.Ni && left.Nj == right.Nj) return true;
		return false;
    }
};

struct dijkVisit 
{
	int N;
	float disTQ;
};
//?????是不是按照最小顺序排列
struct dVComparison
{
    bool operator () (const dijkVisit& left, const dijkVisit& right) const
    {
		//return left.disToQuery > right.disToQuery;
		return left.disTQ > right.disTQ;
    }
};
typedef	priority_queue<dijkVisit ,vector<dijkVisit>,dVComparison> dVQueue;


struct point
{
    int Ni,Nj,pos;
};

//Modified by Qin Xu
//Denote POI on RoadEdge

// record the information of POIcc
struct InerNode
{
	int poid;
    float dis;
	float attr[ATTRIBUTE_DIMENSION];
    //unsigned long long vct;//Vector of keywords denoted by 64-bit
	int nOfK;
	set<int> kwd;
};
// record the edge information, identical to point file
struct edge
{
    int FirstRow;
    int Ni,Nj;
    float dist;
    FastArray<InerNode> pts;
    //float dLB,dUB;		// for gendata use only
	//-----------for extend egtree
	set<int> kwds;
	float attrBound[ATTRIBUTE_DIMENSION][2];
};

typedef map<int,edge*> EdgeMapType; // map node to edges

// build AdjList on fly
FastArray<int>* AdjList;
FastArray<point> PtList;
EdgeMapType EdgeMap;	// key: i0*NodeNum+j0

// get the key of edge
inline int getKey(int i,int j)
{
    int i0=i<j?i:j,j0=j<i?i:j;	// be careful of the order in other places
    return (i0*NodeNum+j0); // map this edge to unique key
}
// break the index of key of edge
inline void breakKey(int key,int& Ni,int& Nj)
{
    Ni=key/NodeNum;
    Nj=key%NodeNum;
}

void printEdgeKeys()
{
    int Ni,Nj;
    printf("EdgeKeys:\n");
    EdgeMapType::iterator p=EdgeMap.begin();
    while (p!=EdgeMap.end())
    {
        edge* e=p->second;
        breakKey(p->first,Ni,Nj);
        printf("%d %d %d\n",Ni,Nj,(e==NULL));
        p++;
    }
}

#endif // __NET_SHARED


