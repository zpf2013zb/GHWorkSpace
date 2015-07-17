/*
 * =====================================================================================
 *
 *       Filename:  querytest.cc
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  05/06/2013 03:15:53 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Qin Xu
 *   Organization:
 *
 * =====================================================================================
 */


#include <iostream>
#include "diskbased.h"
#include "ConfigType.h"
#include "KeywordsGenerator.h"


using namespace std;

// handle
//Edge must existed on Road Network
void getRandEdge(int &Ni,int &Nj,float &Edgedist)
{
    int adjgrpaddr,adjsize=0;
    do
    {
		// get the random edge 
        Ni=random()%NodeNum;
        adjgrpaddr=getAdjListGrpAddr(Ni);
		// adjsize record the number of adjacency edges
        getFixedF(SIZE_A,Ref(adjsize),adjgrpaddr);
        if(adjsize>0)
        {
            int i=random()%adjsize;
            getVarE(ADJNODE_A,Ref(Nj),adjgrpaddr,i);
            getVarE(DIST_A,Ref(Edgedist),adjgrpaddr,i);
        }
    }
    while(adjsize==0);

    if(Ni>Nj)
    {
        int temp=Ni;
        Ni=Nj;
        Nj=temp;
    }

}
// generate the random query Q
void genRandQ(struct QueryPoint &Q,unsigned long long keywordsSet,int topk)
{
    float edgedist=0.0;
    getRandEdge(Q.Ni,Q.Nj,edgedist);
    Q.dist_Ni=drand48()*edgedist;
    Q.k=topk;
    Q.keywords= keywordsSet;
}
// put all the query points to file
void WriteQuery2File(vector<QueryPoint> Qset,const char* filename)
{
    char tmpfilename[255];
    sprintf(tmpfilename,"%s",filename);
	// delete the tmpfilename file
    remove(tmpfilename);
    FILE *f=fopen(tmpfilename,"w");
    if (f) {
        for(QueryPoint Q:Qset)
        {
            char stringline[255];
            sprintf(stringline,"%d %d %f %llu %d\n",Q.Ni,Q.Nj,Q.dist_Ni,Q.keywords,Q.k);
            fputs(stringline,f);
        }
        fclose(f);
    }
    else
    {
        cerr<<"File open error.";
    }
    
}
// open and close the disk to get the basic information and write query to file 
int main(int argc,char** argv)
{
//  Query File default to map_queryPoints_querykeywordsnumber_k_cachepages
    string configFileName = "config.prop";
    ConfigType cr(configFileName,argc, argv);
    cr.ListConfig();
    InitClock();
    OpenDiskComm(cr.getDataFileName().c_str(),cr.getParameterCachePages());

    QueryPoint Q;
    vector<QueryPoint> Qset;
    
    string queryFileName  = cr.getQueryFileName();
    int numQueryPoint = cr.getParameterNumberOfQueryPoints();
    int numQueryKeywords = cr.getParameterQueryKeywordNumbers();
    int k = cr.getParameterK();
    vector<unsigned long long> key = KeywordsGenerator::Instance().getConstantKeywords(numQueryPoint, numQueryKeywords);
    for(unsigned i = 0; i < numQueryPoint; i++)
    {
        genRandQ(Q,key[i],k);
        cout<<Q<<endl;
        Qset.push_back(Q);
    }

    WriteQuery2File(Qset,queryFileName.c_str());
    CloseDiskComm();
    PrintElapsed();
    return 0;
}


