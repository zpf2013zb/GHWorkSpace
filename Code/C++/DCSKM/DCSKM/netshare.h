#ifndef __NET_SHARED
#define __NET_SHARED

#include <queue>
#include <bitset>

#define MAX_DIST 9999999
int MAX_KEYWORDS = 64;

int NodeNum;
int EdgeNum;
struct QueryPoint
{
    int k;
    unsigned long long keywords;
    int Ni,Nj;
    float dist_Ni;
};


ostream& operator<<(ostream& os,const QueryPoint& Q)
{
    os<<"("<<Q.Ni<<","<<Q.Nj<<")    ";
    os<<bitset<64>(Q.keywords).to_string();
    os<<"   "<<Q.k<<" "<<Q.dist_Ni;
    return os;
    
}

struct POI
{
    int Ni,Nj;
    int pre;//Refine use only
    unsigned long long keywords;
    float dist_Ni;
    float dist_Nj;
    float dist_toquery;
    
};


struct POIComparison
{
    bool operator () (const POI& left, const POI& right) const
    {
        return left.dist_toquery> right.dist_toquery;
    }
};


ostream& operator<<(ostream& os, const POI& poi)
{
    os<<"("<<poi.Ni<<","<<poi.Nj<<") ";
    os<<bitset<64>(poi.keywords).to_string();
    os<<" distNi:"<<poi.dist_Ni;
    os<<" distQ:"<<poi.dist_toquery;
    os<<" preNode:"<<poi.pre<<endl;
    
    return os;
}

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

struct point
{
    int Ni,Nj,pos;
};

//Modified by Qin Xu
//Denote POI on RoadEdge
struct InerNode
{
    float dis;
    unsigned long long vct;//Vector of keywords denoted by 64-bit
};

struct edge
{
    int FirstRow;
    int Ni,Nj;
    float dist;
    FastArray<InerNode> pts;
    //float dLB,dUB;		// for gendata use only
};

typedef map<int,edge*> EdgeMapType;

// build AdjList on fly
FastArray<int>* AdjList;
FastArray<point> PtList;
EdgeMapType EdgeMap;	// key: i0*NodeNum+j0

inline int getKey(int i,int j)
{
    int i0=i<j?i:j,j0=j<i?i:j;	// be careful of the order in other places
    return (i0*NodeNum+j0);
}

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


