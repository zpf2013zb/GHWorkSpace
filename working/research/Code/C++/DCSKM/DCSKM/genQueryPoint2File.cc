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
#include <cstdlib>
#include "diskbased.h"
#include "ConfigType.h"
#include "KeywordsGenerator.h"
using namespace std;

#define DIS_CST 1000 //定义距离约束
#define MAX_KWDN 5
#define MAX_KWD 10000 //定义关键字编号
// handle
//Edge must existed on Road Network
//----------------------M--midify random to rand-----------
void getRandEdge(int &Ni,int &Nj,float &Edgedist)
{
    int adjgrpaddr,adjsize=0;
    do
    {
		// get the random edge 
        //Ni=random(INT_MAX)%NodeNum;
		Ni = rand() % NodeNum + 1;
        adjgrpaddr=getAdjListGrpAddr(Ni);
		// adjsize record the number of adjacency edges
        getFixedF(SIZE_A,Ref(adjsize),adjgrpaddr);
        if(adjsize>0)
        {
            int i=rand()%adjsize;
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
//----------------------M--midify random to rand-----------
void genSubspace(bitset<ATTRIBUTE_DIMENSION> &ss) {

	int seq = 0;
	int seqMax = (int)pow(2, ATTRIBUTE_DIMENSION);
	while (seq == 0) {
		seq = rand() % seqMax;
	}
	int i = 0;
	while (seq)
	{
		ss[i] = seq % 2;
		seq = seq / 2;
		i++;
	}
}
//-------------------------M--修改了Q中的数据格式
void genRandQ(struct QueryPoint &Q)
{
    float edgedist=0.0;
    getRandEdge(Q.Ni,Q.Nj,edgedist);
    Q.dist_Ni=(rand()+1)*edgedist/RAND_MAX;
    Q.distCnst= (rand() + 1)*DIS_CST / RAND_MAX;
	Q.nOfKwd = rand()% MAX_KWDN + 1;
	for (int i = 0; i < Q.nOfKwd; i++) {
		Q.kwd.insert(rand() % MAX_KWD + 1);
	}
	genSubspace(Q.subSpace);
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
            sprintf(stringline,"%d %d %f %f %d",Q.Ni,Q.Nj,Q.dist_Ni,Q.distCnst,Q.nOfKwd);
			set<int> ::iterator it = Q.kwd.begin();
			for (; it != Q.kwd.end(); it++) {
				sprintf(stringline, " %d", *it);
			}
			unsigned long ul = Q.subSpace.to_ulong;
			sprintf(stringline, " %ul\n", ul);
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
        genRandQ(Q,key[i]);
        cout<<Q<<endl;
        Qset.push_back(Q);
    }

    WriteQuery2File(Qset,queryFileName.c_str());
    CloseDiskComm();
    PrintElapsed();
    return 0;
}


