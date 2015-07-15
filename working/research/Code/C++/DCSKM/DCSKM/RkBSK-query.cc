/*
 * =====================================================================================
 *
 *       Filename:  RkBSK-query.cc
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  05/06/2013 03:15:53 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Qin Xu
 *   Organization:  ZJU
 *
 * =====================================================================================
 */


#include <iostream>
#include <unordered_map>
#include <queue>
#include <bitset>
#include "diskbased.h"
#include "ConfigType.h"
#include <fstream>

using namespace std;

#define BasicCountMap unordered_map<unsigned long long,int> //map<keywords,count>
#define EnhancedCountMap unordered_map<unsigned long long,pair<int,int> > //map<keywords,<count,specialcount> >
#define SubsetMap unordered_map<unsigned long long,vector<unsigned long long> > //map<keywords,subset of keywords vector>

#define CandidateSet vector<POI>
#define ResultSet vector<POI>


enum Querytype {baseline1,baseline2,enhanced};


SubsetMap subsmap;
unordered_map<int ,BasicCountMap > NodeBCM_map;//NodeBCL_map[Ni]==BasicCountlist of Ni
unordered_map<int ,EnhancedCountMap > NodeECM_map;
FastArray<float> NodeDisttoQ;
FastArray<bool> isVisted;
FastArray<int> NodePre;


//Experimental Evaluation use variable
int64_t base1datapageaccessed=0,base2datapageaccessed=0,enhanceddatapageaccessed=0;
int64_t base1poiverified=0,base2poiverified=0,enhancedpoiverified=0;
int64_t l1prunned=0,l2prunned=0,l3prunned=0,l4prunned=0;
int64_t l3prunnedperquery=0,l4prunnedperquery=0;
double totalbase1=0,totalbase2=0,totalenhanced=0;
int64_t l1appliednum=0,l2appliednum=0,l3appliednum=0,l4appliednum=0;
int64_t base2numnodeexpand=0,base2numedgeexpand=0,enhancednumnodeexpand=0,enhancednumedgeexpand=0;

void initCount()
{
    base1datapageaccessed=0,base2datapageaccessed=0,enhanceddatapageaccessed=0;
    base1poiverified=0,base2poiverified=0,enhancedpoiverified=0;
    l1prunned=0,l2prunned=0,l3prunned=0,l4prunned=0;
    l3prunnedperquery=0,l4prunnedperquery=0;
    totalbase1=0,totalbase2=0,totalenhanced=0;
    l1appliednum=0,l2appliednum=0,l3appliednum=0,l4appliednum=0;
    base2numnodeexpand=0,base2numedgeexpand=0,enhancednumnodeexpand=0,enhancednumedgeexpand=0;
}

CandidateSet cS;//Store Candidate POI during the query procedure
ResultSet rS;//Store Final result POI

void RefineDist(CandidateSet &cS,const QueryPoint &Q);//Refine distance of POI

int getCount(unsigned long long &keyword)
{
    int count=0;
    int64_t key=1;
    for(int i=0; i<64; i++)
    {
        if((key<<i)&keyword)
            count++;
    }
    return count;
}

void genSubset(unsigned long long keywordset,vector<unsigned long long > &subset)
{
    bool flag=false;
    for(auto iter=subset.begin(); iter!=subset.end(); iter++)
        if(*iter==keywordset) flag=true;
    if(!flag)subset.push_back(keywordset);
    if(getCount(keywordset)==1)
    {

        return ;
    }
    else
    {
        for(int i=0; i<sizeof(unsigned long long)*8; i++)
        {
            if((int64_t(1)<<i)&keywordset)
            {
                unsigned long long tmpkeywordset=keywordset&(~(int64_t(1)<<i));
                genSubset(tmpkeywordset,subset);
            }
        }
    }
}

void printSubset(vector<unsigned long long > &subset)
{
    auto iter=subset.begin();
    for(; iter!=subset.end(); iter++)
    {
        cout<<*iter<<endl;
    }
}

void genSubsetMap(SubsetMap &subsmap,const QueryPoint &Q)
{
    vector<unsigned long long> subsetofQ;
    genSubset(Q.keywords,subsetofQ);
    subsmap.insert(pair<unsigned long long,vector<unsigned long long> >(Q.keywords,subsetofQ));
    vector<unsigned long long>::iterator iter=subsetofQ.begin();
    while(iter!=subsetofQ.end())
    {
        vector<unsigned long long> tmp;
        genSubset(*iter,tmp);
        subsmap.insert(pair<unsigned long long,vector<unsigned long long> >(*iter,tmp));
        iter++;
    }
    //cout<<"genSubsetMap Done!"<<endl;
    return ;
}



void showAllAdj(int NodeID)
{

    int Ni=NodeID,AdjGrpAddr,AdjListSize,NewNodeID,PtGrpKey;
    float Egdist;

    AdjGrpAddr=getAdjListGrpAddr(Ni);
    getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);
    cout<<Ni<<"'s Adj # = "<<AdjListSize<<endl;

    for (int i=0; i<AdjListSize; i++)
    {
        getVarE(ADJNODE_A,Ref(NewNodeID),AdjGrpAddr,i);
        getVarE(DIST_A,Ref(Egdist),AdjGrpAddr,i);
        getVarE(PTKEY_A,Ref(PtGrpKey),AdjGrpAddr,i);
        cout<<"AdjNodeID: "<<NewNodeID<<" AdjDist: "<<Egdist<<" AdjPtKey: "<<PtGrpKey<<endl;
    }

}

bool LcontianR(unsigned long long L,unsigned long long R)
{
    //High bit set to zero
    //R&=tmp;
    return  (L&R)==R;
}

bool LintersectR(unsigned long long L,unsigned long long R,unsigned long long &keywords)
{
    if((L&R)!=0)
    {
        keywords=(L&R);
        return true;
    }
    else
        return false;

}


bool INEBasedVerify(POI p,float dist_toQ,int k)//Need to remove p itself on edge where p located
{
    //cout<<"Verifing :"<<endl;
    //printPOI(p);

    vector<bool> visit;
    visit.assign(NodeNum,false);
    vector<float> distance_p;
    distance_p.assign(NodeNum,MAX_DIST);
    DStepQueue sQ;
    int cnt=0;//Current Top-k #
    bool flag=false;
    //Data needed
    int Ni=p.Ni,Nj=p.Nj,AdjGrpAddr,AdjListSize,NewNodeID,PtGrpKey,PtNumOnEdge;
    unsigned long long keywords;
    float EdgeDist,PtDist;

    //Precess sub edge of p.Ni,p.Nj divided by p
    AdjGrpAddr=getAdjListGrpAddr(Ni);
    getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);
    for (int i=0; i<AdjListSize; i++)
    {
        getVarE(ADJNODE_A,Ref(NewNodeID),AdjGrpAddr,i);
        //getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
        if(NewNodeID == Nj)
        {
            getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
            getVarE(PTKEY_A,Ref(PtGrpKey),AdjGrpAddr,i);
            //cout<<"("<<Ni<<","<<Nj<<") EdgeDist:"<<EdgeDist;
            if(PtGrpKey==-1)
            {
                //cout<<"No POI existed on Edge where p located"<<endl;
            }
            else
            {
                getFixedF(SIZE_P,Ref(PtNumOnEdge),PtGrpKey);
                //cout<<" PtNumOnEdge:"<<PtNumOnEdge<<endl;
                for(int j=0; j<PtNumOnEdge; j++)
                {
                    getVarE(PT_P,Ref(PtDist),PtGrpKey,j);
                    if((!flag)&&fabs(PtDist-p.dist_Ni)<1e-6)
                    {
                        flag=true;    //Remove p itself
                        continue;
                    }
                    getVarE(PT_VCT,Ref(keywords),PtGrpKey,j);
                    if(fabs(PtDist-p.dist_Ni)<dist_toQ&&LcontianR(keywords,p.keywords)&&cnt<k)
                    {
                        cnt++;
                        //cout<<"One result found on Edge p located"<<endl;
                        //cout<<"POI dist to p:"<<PtDist-p.dist_Ni<<"POI keywords:";
                        //printBinary(keywords);
                        //cout<<endl;
                    }
                    //cout<<"InerNode "<<j<<" Dist:"<<PtDist<<endl;
                    //cout<<"KeywordsVct:";
                    //printBinary(keywords);
                    //cout<<"KeywordsVct#:"<<keywords<<endl;
                }
                if(cnt>=k)
                    return false;
            }

            if(p.dist_Ni<dist_toQ)
            {
                DStepEvent event;
                event.node=Ni;
                event.dist=p.dist_Ni;
                sQ.push(event);


            }
            visit[p.Ni]=true;
            distance_p[Ni]=p.dist_Ni;
            if((EdgeDist-p.dist_Ni)<dist_toQ)
            {
                DStepEvent event;
                event.node=Nj;
                event.dist=EdgeDist-p.dist_Ni;
                sQ.push(event);


            }
            visit[p.Nj]=true;
            distance_p[Nj]=EdgeDist-p.dist_Ni;
            break;

        }
    }

//Subedge process ended here!
//Begin INE if sQ.empty()!=Null
//temppoi to handle those poi satisfy keywords constraint but tempdist > dist_toQ
    vector<POI> temppoi;
    while(!sQ.empty())
    {

        DStepEvent event=sQ.top();
        sQ.pop();

        if(visit[event.node]&&event.node!=p.Ni&&event.node!=p.Nj) continue;
        visit[event.node]=true;

        //for each non-visited adj. do
        AdjGrpAddr=getAdjListGrpAddr(event.node);
        getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);
        for (int i=0; i<AdjListSize; i++)
        {

            getVarE(ADJNODE_A,Ref(NewNodeID),AdjGrpAddr,i);
            //getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
            if(!visit[NewNodeID])
            {
                getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
                getVarE(PTKEY_A,Ref(PtGrpKey),AdjGrpAddr,i);
                //cout<<"("<<Ni<<","<<Nj<<") EdgeDist:"<<EdgeDist;
                if(PtGrpKey!=-1)
                {
                    getFixedF(SIZE_P,Ref(PtNumOnEdge),PtGrpKey);
                    //cout<<" PtNumOnEdge:"<<PtNumOnEdge<<endl;
                    for(int j=0; j<PtNumOnEdge; j++)
                    {
                        getVarE(PT_P,Ref(PtDist),PtGrpKey,j);
                        getVarE(PT_VCT,Ref(keywords),PtGrpKey,j);
                        float poitop;
                        if(event.node<NewNodeID) poitop=event.dist+PtDist;
                        else poitop=EdgeDist-PtDist+event.dist;

                        //poitop may not be correct dist
                        //correct answer must less than current poitop
                        if(LcontianR(keywords,p.keywords))
                        {
                            if(poitop<dist_toQ)
                            {
                                cnt++;//another lemma???Exist one path to p less than dist(p,Q)
                                //cout<<"One result found on Edge:"<<"("<<event.node<<","<<NewNodeID<<")"<<endl;
                                //cout<<"POI dist to p:"<<poitop<<endl;
                                //cout<<"POI keywords:";
                                //printBinary(keywords);
                                //cout<<endl;
                            }
                            else
                            {
                                //cout<<"Special Case Happens!!!!"<<endl;
                                POI tmp;
                                tmp.Ni=event.node<NewNodeID?event.node:NewNodeID;
                                tmp.Nj=event.node<NewNodeID?NewNodeID:event.node;
                                tmp.dist_Ni=PtDist;
                                tmp.dist_Nj=EdgeDist-PtDist;
                                tmp.dist_toquery=poitop;
                                temppoi.push_back(tmp);
                            }

                        }

                        //cout<<"InerNode "<<j<<" Dist:"<<PtDist<<endl;
                        //cout<<"KeywordsVct:";
                        //printBinary(keywords);
                        //cout<<"KeywordsVct#:"<<keywords<<endl;
                    }
                    if(cnt>=k)
                        return false;
                }
                DStepEvent newevent=event;
                newevent.node=NewNodeID;
                newevent.dist=event.dist+EdgeDist;
                if (newevent.dist<distance_p[NewNodeID])
                {
                    if(newevent.dist<dist_toQ)
                        sQ.push(newevent);
                    distance_p[NewNodeID]=newevent.dist;
                    //cout<<newevent.dist;
                }
            }
        }

    }

    while(!temppoi.empty())//get exact dist of poi to p in this set
    {
        POI tmp;
        tmp=temppoi.back();
        if(tmp.dist_Ni+distance_p[tmp.Ni]<dist_toQ)
        {
            cnt++;
            //cout<<"special case!!!"<<endl;
            //cout<<tmp.dist_Ni+distance_p[tmp.Ni]<<endl;
        }
        if(tmp.dist_Nj+distance_p[tmp.Nj]<dist_toQ)
        {
            cnt++;
            //cout<<"special case!!!"<<endl;
            //cout<<tmp.dist_Nj+distance_p[tmp.Nj]<<endl;
        }
        temppoi.pop_back();
        if(cnt>=k) return false;
    }

    //cout<<"True Hit!"<<endl;
    return true;
}


//InitCountMap of Ni,Nj where Q located,must be Invoked after genSubsetMap
void InitCountMap(const struct QueryPoint &Q,Querytype type)
{

    if(type==baseline2)
    {
        BasicCountMap bcmNi,bcmNj;
        vector<unsigned long long>::iterator iter=subsmap[Q.keywords].begin();
        while(iter!=subsmap[Q.keywords].end())
        {
            bcmNi.insert(pair<unsigned long long,int>(*iter,Q.k));//keywords,count
            //cout<<"bcmNi "<<*iter<<" "<<Q.k<<"Inserted."<<endl;
            bcmNj.insert(pair<unsigned long long,int>(*iter,Q.k));
            //cout<<"bcmNj "<<*iter<<" "<<Q.k<<"Inserted."<<endl;
            iter++;
        }
        NodeBCM_map.insert(pair<int,BasicCountMap>(Q.Ni,bcmNi));
        NodeBCM_map.insert(pair<int,BasicCountMap>(Q.Nj,bcmNj));
        //cout<<"InitCountMap Done!"<<endl;

    }
    else if(type==enhanced)
    {
        EnhancedCountMap ecmNi,ecmNj;
        vector<unsigned long long>::iterator iter=subsmap[Q.keywords].begin();
        while(iter!=subsmap[Q.keywords].end())
        {
            ecmNi.insert(pair<unsigned long long,pair<int,int> >(*iter,pair<int,int>(Q.k,0)));
            ecmNj.insert(pair<unsigned long long,pair<int,int> >(*iter,pair<int,int>(Q.k,0)));
            iter++;
        }
        NodeECM_map.insert(pair<int,EnhancedCountMap>(Q.Ni,ecmNi));
        NodeECM_map.insert(pair<int,EnhancedCountMap>(Q.Nj,ecmNj));
        //cout<<"InitCountMap Done!"<<endl;
    }
    else
    {
        //cout<<"This query type do not need InitCountMap."<<endl;
    }

    return;
}

void InitQuery(const QueryPoint& Q,Querytype type)//forBaseline method 2 and Enhanced method
{
    RefreshStat();

    if(type==baseline2||type==enhanced)
    {
        //Must clear Attention!!!!
        l3prunnedperquery=0;
        l4prunnedperquery=0;
        subsmap.clear();
        NodeBCM_map.clear();
        NodeECM_map.clear();
        genSubsetMap(subsmap,Q);
        InitCountMap(Q,type);
    }
    NodePre.clear();    //Lemma 4 use only
    NodePre.assign(NodeNum,-1);
    cS.reserve(PtNum);
    rS.reserve(PtNum);
    cS.clear();
    rS.clear();
    NodeDisttoQ.assign(NodeNum,MAX_DIST);
    isVisted.assign(NodeNum,false);
    //cout<<"Init Query time:"<<clock()-t1<<endl;
}


//Spread Road Network using INE
//Add poi satisfy keywords constraint to cS
void INEBasedKeywordFilter(const QueryPoint &Q)
{
    DStepQueue sQ;
    int poicnt=0;
    int roadnodecnt=0;

    //Data
    int Ni=Q.Ni,Nj=Q.Nj;
    int AdjGrpAddr,AdjListSize,NewNodeID,PtGrpKey,PtNumOnEdge;
    unsigned long long keywords;
    float EdgeDist,PtDist;

    //Precess sub edge of Q.Ni,Q.Nj divided by Q
    AdjGrpAddr=getAdjListGrpAddr(Ni);
    getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);
    for (int i=0; i<AdjListSize; i++)
    {
        getVarE(ADJNODE_A,Ref(NewNodeID),AdjGrpAddr,i);
        //getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
        if(NewNodeID == Nj)
        {
            getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
            getVarE(PTKEY_A,Ref(PtGrpKey),AdjGrpAddr,i);
            //cout<<"("<<Ni<<","<<Nj<<") EdgeDist:"<<EdgeDist;
            if(PtGrpKey==-1)
            {
                //cout<<"No POI existed on Edge where Q located."<<endl;
            }
            else
            {
                getFixedF(SIZE_P,Ref(PtNumOnEdge),PtGrpKey);
                //cout<<" PtNum on Edge where Q located:"<<PtNumOnEdge<<endl;
                for(int j=0; j<PtNumOnEdge; j++)
                {
                    poicnt++;
                    getVarE(PT_P,Ref(PtDist),PtGrpKey,j);
                    getVarE(PT_VCT,Ref(keywords),PtGrpKey,j);
                    if(LcontianR(Q.keywords,keywords))
                    {
                        POI tmp;
                        tmp.dist_toquery=fabs(PtDist-Q.dist_Ni);
                        tmp.keywords=keywords;
                        tmp.dist_Ni=PtDist;
                        tmp.dist_Nj=EdgeDist-PtDist;
                        tmp.Ni=Ni;
                        tmp.Nj=Nj;
                        tmp.pre=-2;
                        /*
                        if(PtDist<Q.dist_Ni)
                        {
                            tmp.dist_Ni=fabs(PtDist-Q.dist_Ni);
                            tmp.dist_Nj=PtDist;
                            tmp.Ni=-2;
                            tmp.Nj=Ni;
                            tmp.pre=-2;
                        }
                        else
                        {
                            tmp.dist_Ni=fabs(PtDist-Q.dist_Ni);
                            tmp.dist_Nj=EdgeDist-PtDist;
                            tmp.Ni=-2;
                            tmp.Nj=Nj;
                            tmp.pre=-2;
                        }
                        */
                        //cout<<"POI dist to Q:"<<tmp.dist_toquery<<"  POI keywords:";
                        //printBinary(keywords);
                        //cout<<endl;
                        cS.push_back(tmp);
                    }
                }

            }

            DStepEvent event;
            event.node=Ni;
            event.dist=Q.dist_Ni;
            sQ.push(event);
            NodeDisttoQ[Ni]=Q.dist_Ni;
            event.node=Nj;
            event.dist=EdgeDist-Q.dist_Ni;
            sQ.push(event);
            NodeDisttoQ[Nj]=event.dist;
            isVisted[Q.Ni]=true;
            isVisted[Q.Nj]=true;
            NodePre[Q.Ni]=-2;//-2 denote Q
            NodePre[Q.Nj]=-2;
            break;

        }
    }
    while(!sQ.empty())
    {
        DStepEvent event=sQ.top();
        sQ.pop();
        if(event.node!=Q.Ni&&event.node!=Q.Nj&&isVisted[event.node]) continue;
        isVisted[event.node]= true;
        roadnodecnt++;
        //for each non-visited adj. do
        AdjGrpAddr=getAdjListGrpAddr(event.node);
        getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);
        for (int i=0; i<AdjListSize; i++)
        {

            getVarE(ADJNODE_A,Ref(NewNodeID),AdjGrpAddr,i);
            //getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
            if(!isVisted[NewNodeID])
            {

                getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
                getVarE(PTKEY_A,Ref(PtGrpKey),AdjGrpAddr,i);
                //cout<<"("<<Ni<<","<<Nj<<") EdgeDist:"<<EdgeDist;
                if(PtGrpKey!=-1)
                {
                    getFixedF(SIZE_P,Ref(PtNumOnEdge),PtGrpKey);
                    //cout<<" PtNumOnEdge:"<<PtNumOnEdge<<endl;
                    for(int j=0; j<PtNumOnEdge; j++)
                    {
                        poicnt++;
                        getVarE(PT_P,Ref(PtDist),PtGrpKey,j);
                        getVarE(PT_VCT,Ref(keywords),PtGrpKey,j);
                        //for all poi who satisfy keywords constraint added it to CandidateSet
                        if(LcontianR(Q.keywords,keywords))
                        {
                            POI tmp;
                            tmp.keywords=keywords;
                            tmp.Ni=event.node<NewNodeID?event.node:NewNodeID;
                            tmp.Nj=event.node<NewNodeID?NewNodeID:event.node;
                            tmp.dist_Ni=PtDist;
                            tmp.dist_Nj=EdgeDist-PtDist;
                            tmp.pre=event.node;
                            cS.push_back(tmp);
                        }

                    }

                }
                DStepEvent newevent=event;
                newevent.node=NewNodeID;
                newevent.dist=event.dist+EdgeDist;
                if (newevent.dist<NodeDisttoQ[NewNodeID])
                {
                    NodePre[NewNodeID]=event.node;
                    NodeDisttoQ[NewNodeID]=newevent.dist;
                    sQ.push(newevent);
                }
            }
        }
    }
    //cout<<"POI Read #"<<poicnt<<endl;
    //cout<<"RoadNode Read #"<<roadnodecnt<<endl;
}

//Show all count or special # of Node Ni
//For test purpose
void showCountMap(int NodeID,Querytype type)
{
    cout<<"NodeID:"<<NodeID<<endl;
    if(type==baseline2)
    {
        if(NodeBCM_map.find(NodeID)==NodeBCM_map.end())
        {
            cout<<"No CountMap existed!"<<endl;
            return;
        }
        BasicCountMap::iterator iter=NodeBCM_map[NodeID].begin();
        while(iter!=NodeBCM_map[NodeID].end())
        {
            cout<<"Subset of Keywords:"<<(*iter).first<<"  Count #:"<<(*iter).second<<endl;
            iter++;
        }
    }
    else if(type==enhanced)
    {
        if(NodeECM_map.find(NodeID)==NodeECM_map.end())
        {
            cout<<"No CountMap existed!"<<endl;
            return;
        }
        EnhancedCountMap::iterator iter=NodeECM_map[NodeID].begin();
        while(iter!=NodeECM_map[NodeID].end())
        {
            cout<<"Subset of Keywords:"<<(*iter).first;
            cout<<"  Count #:"<<(*iter).second.first<<"  Special #:"<<(*iter).second.second<<endl;
            iter++;
        }
    }
    else
    {
        cout<<"This Query type do not need CountMap."<<endl;
    }
}

//Only satisfy keywords constraint should this fun be invoked
//for Baseline method 2 must contain
//for Enhanced method keywords intersection with Q must be not null
void updateBasicCountMap(int NodeID,unsigned long long keywords)
{
    NodeBCM_map[NodeID][keywords]--;
    if(NodeBCM_map[NodeID][keywords]==0)
    {
        NodeBCM_map[NodeID].erase(keywords);
        //showCountMap(NodeID,Querytype(1));
    }
    /*
    if(!NodeBCM_map[NodeID].size())
        cout<<"End Expansion in this direction."<<endl;
    */
}

void updateEnhancedCountMap(int NodeID,unsigned long long keywords,bool flag)//flag donate the dist to event.node true or false
{
    vector<unsigned long long>::iterator iter=subsmap[keywords].begin();
    if(flag)
    {
        while(iter!=subsmap[keywords].end())
        {
            if(NodeECM_map[NodeID].find(*iter)!=NodeECM_map[NodeID].end())
            {
                NodeECM_map[NodeID][*iter].first--;
                NodeECM_map[NodeID][*iter].second++;
                if(NodeECM_map[NodeID][*iter].first<=0)
                    NodeECM_map[NodeID].erase(*iter);
            }
            iter++;
        }
    }
    else
    {
        while(iter!=subsmap[keywords].end())
        {
            if(NodeECM_map[NodeID].find(*iter)!=NodeECM_map[NodeID].end())
            {
                NodeECM_map[NodeID][*iter].first--;
                if(NodeECM_map[NodeID][*iter].first<=0)
                    NodeECM_map[NodeID].erase(*iter);
            }
            iter++;
        }
    }
}



//using ecm0 to effect ecm1
void updateEnhancedCountMap(EnhancedCountMap &ecm0,EnhancedCountMap &ecm1)
{
    //l3applied++;
    for(auto &key:ecm0)
    {
        if((key.second.second)&&(ecm1.find(key.first)!=ecm1.end()))
        {
            //Update Keywords Counts using special count

            ecm1[key.first].first-=key.second.second;
            l3prunnedperquery+=key.second.second;

            if(ecm1[key.first].second)
            {
                key.second.first-=ecm1[key.first].second;
                l3prunnedperquery+=ecm1[key.first].second;
            }

        }

    }
    /*
        for(auto &key:ecm1)

            if((key.second.second)&&(ecm0.find(key.first)!=ecm0.end()))
            {
                //Update Keywords Counts using special count
                ecm0[key.first].first-=key.second.second;
            }
    */
}

//Baseline method 1
void Baselinemethod1(const QueryPoint &Q)
{
    InitQuery(Q,Querytype(0));
    INEBasedKeywordFilter(Q);
    RefineDist(cS,Q);
    CandidateSet::iterator iter=cS.begin();
    //int PtVerified=0;
    while(iter!=cS.end())
    {
        //PtVerified++;
        //Must refine dist of POI here!!
        /*
        if((*iter).dist_toquery>(*iter).dist_Ni+NodeDisttoQ[(*iter).Ni])
        {
            specialcnt++;
        }
        if((*iter).dist_toquery>(*iter).dist_Nj+NodeDisttoQ[(*iter).Nj])
        {
            specialcnt++;
        }

        float refineddist=(*iter).dist_Ni+NodeDisttoQ[(*iter).Ni]<(*iter).dist_Nj+NodeDisttoQ[(*iter).Nj]?(*iter).dist_Ni+NodeDisttoQ[(*iter).Ni]:(*iter).dist_Nj+NodeDisttoQ[(*iter).Nj];
        */

        if(INEBasedVerify(*iter,(*iter).dist_toquery,Q.k))
        {
            // (*iter).dist_toquery=refineddist;
            rS.push_back(*iter);
        }
        iter++;
    }
    //cout<<"Total # of POI Verified:"<<PtVerified<<endl;
}


void INEBasedSpatialKeywordFilter(const QueryPoint& Q)
{
    DStepQueue sQ;
    //float Dmax=0.0;
    int poicnt=0;
    int roadnodecnt=0;
    int edgeexpanded=1;
    //Data
    int Ni=Q.Ni,Nj=Q.Nj;//Ni must less than Nj
    int AdjGrpAddr,AdjListSize,NewNodeID,PtGrpKey,PtNumOnEdge;
    unsigned long long keywords;
    float EdgeDist,PtDist;
    //Precess sub edge of Q.Ni,Q.Nj divided by Q
    AdjGrpAddr=getAdjListGrpAddr(Ni);
    getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);
    for (int i=0; i<AdjListSize; i++)
    {
        getVarE(ADJNODE_A,Ref(NewNodeID),AdjGrpAddr,i);
        //getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
        if(NewNodeID == Nj)
        {
            getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
            getVarE(PTKEY_A,Ref(PtGrpKey),AdjGrpAddr,i);
            vector<POI> poinj;
            vector<POI> poini;
            //cout<<"("<<Ni<<","<<Nj<<") EdgeDist:"<<EdgeDist;
            if(PtGrpKey==-1)
            {
                //cout<<"No POI existed on Edge where Q located."<<endl;
            }
            else
            {
                getFixedF(SIZE_P,Ref(PtNumOnEdge),PtGrpKey);
                //cout<<" PtNum on Edge where Q located:"<<PtNumOnEdge<<endl;
                //Notice the order POI visited on edge
                for(int j=0; j<PtNumOnEdge; j++)
                {
                    poicnt++;
                    getVarE(PT_P,Ref(PtDist),PtGrpKey,j);
                    getVarE(PT_VCT,Ref(keywords),PtGrpKey,j);
                    if(LcontianR(Q.keywords,keywords))
                    {
                        POI tmp;
                        tmp.Ni=Ni;
                        tmp.Nj=Nj;
                        tmp.dist_toquery=fabs(PtDist-Q.dist_Ni);
                        tmp.dist_Ni=PtDist;
                        tmp.dist_Nj=EdgeDist-PtDist;
                        tmp.keywords=keywords;
                        tmp.pre=-2;

                        if(PtDist<Q.dist_Ni&&NodeBCM_map[Ni].find(keywords)!=NodeBCM_map[Ni].end())
                        {
                            poini.push_back(tmp);

                            //updateBasicCountMap(Ni,keywords);
                            //showCountMap(Ni,Querytype(1));
                            //if(tmp.dist_toquery>Dmax) Dmax=tmp.dist_toquery;
                            //cS.push_back(tmp);
                        }

                        //showCountMap(Nj,Querytype(1));
                        if(PtDist>Q.dist_Ni&&NodeBCM_map[Nj].find(keywords)!=NodeBCM_map[Nj].end())
                        {
                            poinj.push_back(tmp);

                            //updateBasicCountMap(Nj,keywords);
                            //showCountMap(Nj,Querytype(1));
                            //if(tmp.dist_toquery>Dmax) Dmax=tmp.dist_toquery;
                            //cS.push_back(tmp);
                        }
                    }
                }

            }
            sort(poini.begin(),poini.end(),POIComparison());
            sort(poinj.begin(),poinj.end(),POIComparison());

            while(!poinj.empty())
            {
                POI tmp=poinj.back();
                if( NodeBCM_map[Nj].find(tmp.keywords)!=NodeBCM_map[Nj].end())
                {
                    updateBasicCountMap(Nj,tmp.keywords);
                    cS.push_back(tmp);

                }

                poinj.pop_back();
            }
            while(!poini.empty())
            {
                POI tmp=poini.back();
                if( NodeBCM_map[Ni].find(tmp.keywords)!=NodeBCM_map[Ni].end())
                {
                    updateBasicCountMap(Ni,tmp.keywords);
                    cS.push_back(tmp);

                }

                poini.pop_back();
            }
            DStepEvent event;
            event.node=Ni;
            event.dist=Q.dist_Ni;

            if(NodeBCM_map[Ni].size())
                sQ.push(event);//exist count != 0
            else
                NodeBCM_map.erase(Ni);
            NodeDisttoQ[Ni]=Q.dist_Ni;
            //Dmax=event.dist;
            event.node=Nj;
            event.dist=EdgeDist-Q.dist_Ni;
            //if(event.dist>Dmax)Dmax=event.dist;
            if(NodeBCM_map[Nj].size())
                sQ.push(event);
            else
                NodeBCM_map.erase(Nj);
            NodeDisttoQ[Nj]=event.dist;
            isVisted[Q.Ni]=true;//Remove Dual visit of edge where Q located
            isVisted[Q.Nj]=true;
            NodePre[Q.Ni]=-2;//-2 denote Q
            NodePre[Q.Nj]=-2;
            break;
        }
    }


    while(!sQ.empty())
    {
        if(NodeBCM_map.size()==0) break;
        DStepEvent event;
        event=sQ.top();
        sQ.pop();
        if(isVisted[event.node]&&event.node!=Q.Ni&&event.node!=Q.Nj) continue;
        isVisted[event.node] = true;
        roadnodecnt++;
        //for each non-visited adj. do
        AdjGrpAddr=getAdjListGrpAddr(event.node);
        getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);
        for (int i=0; i<AdjListSize; i++)
        {
            getVarE(ADJNODE_A,Ref(NewNodeID),AdjGrpAddr,i);

            //getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
            if(!isVisted[NewNodeID])
            {
                edgeexpanded++;
                getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
                getVarE(PTKEY_A,Ref(PtGrpKey),AdjGrpAddr,i);
                if(NodeBCM_map.find(event.node)!=NodeBCM_map.end()&&NodeBCM_map[event.node].size())
                {
                    if(event.dist+EdgeDist<NodeDisttoQ[NewNodeID])//Only update countmap satisfy some constraint
                    {
                        BasicCountMap tmpbcm=NodeBCM_map[event.node];
                        NodeBCM_map.insert(pair<int,BasicCountMap>(NewNodeID,tmpbcm));
                        if(PtGrpKey!=-1)
                        {
                            getFixedF(SIZE_P,Ref(PtNumOnEdge),PtGrpKey);
                            //cout<<" PtNumOnEdge:"<<PtNumOnEdge<<endl;
                            if(event.node<NewNodeID)//For purpose to visited in order
                            {
                                for(int j=PtNumOnEdge-1; j>=0; j--)
                                {
                                    poicnt++;
                                    getVarE(PT_P,Ref(PtDist),PtGrpKey,j);
                                    getVarE(PT_VCT,Ref(keywords),PtGrpKey,j);
                                    //for all poi who satisfy keywords constraint added it to CandidateSet
                                    if(LcontianR(Q.keywords,keywords)&&NodeBCM_map[NewNodeID].find(keywords)!=NodeBCM_map[NewNodeID].end())
                                    {
                                        updateBasicCountMap(NewNodeID,keywords);
                                        POI tmp;
                                        tmp.keywords=keywords;
                                        tmp.Ni=event.node<NewNodeID?event.node:NewNodeID;
                                        tmp.Nj=event.node<NewNodeID?NewNodeID:event.node;
                                        tmp.dist_Ni=PtDist;
                                        tmp.dist_Nj=EdgeDist-PtDist;
                                        //if(event.node==tmp.Ni&&NodeDisttoQ[tmp.Ni]+PtDist>Dmax) Dmax=NodeDisttoQ[tmp.Ni]+PtDist;
                                        //if(event.node==tmp.Nj&&NodeDisttoQ[tmp.Nj]+EdgeDist-PtDist>Dmax) Dmax=NodeDisttoQ[tmp.Nj]+EdgeDist-PtDist;
                                        cS.push_back(tmp);
                                    }

                                }
                            }
                            else
                            {
                                for(int j=0; j<PtNumOnEdge; j++)
                                {
                                    poicnt++;
                                    getVarE(PT_P,Ref(PtDist),PtGrpKey,j);
                                    getVarE(PT_VCT,Ref(keywords),PtGrpKey,j);
                                    //for all poi who satisfy keywords constraint added it to CandidateSet
                                    if(LcontianR(Q.keywords,keywords)&&NodeBCM_map[NewNodeID].find(keywords)!=NodeBCM_map[NewNodeID].end())
                                    {
                                        updateBasicCountMap(NewNodeID,keywords);
                                        POI tmp;
                                        tmp.keywords=keywords;
                                        tmp.Ni=event.node<NewNodeID?event.node:NewNodeID;
                                        tmp.Nj=event.node<NewNodeID?NewNodeID:event.node;
                                        tmp.dist_Ni=PtDist;
                                        tmp.dist_Nj=EdgeDist-PtDist;
                                        //if(event.node==tmp.Ni&&NodeDisttoQ[tmp.Ni]+PtDist>Dmax) Dmax=NodeDisttoQ[tmp.Ni]+PtDist;
                                        //if(event.node==tmp.Nj&&NodeDisttoQ[tmp.Nj]+EdgeDist-PtDist>Dmax) Dmax=NodeDisttoQ[tmp.Nj]+EdgeDist-PtDist;
                                        cS.push_back(tmp);
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        if(PtGrpKey!=-1)
                        {
                            getFixedF(SIZE_P,Ref(PtNumOnEdge),PtGrpKey);

                            for(int j=0; j<PtNumOnEdge; j++)
                            {
                                poicnt++;
                                getVarE(PT_P,Ref(PtDist),PtGrpKey,j);
                                getVarE(PT_VCT,Ref(keywords),PtGrpKey,j);
                                //for all poi who satisfy keywords constraint added it to CandidateSet
                                if(LcontianR(Q.keywords,keywords))
                                {
                                    POI tmp;
                                    tmp.keywords=keywords;
                                    tmp.Ni=event.node<NewNodeID?event.node:NewNodeID;
                                    tmp.Nj=event.node<NewNodeID?NewNodeID:event.node;
                                    tmp.dist_Ni=PtDist;
                                    tmp.dist_Nj=EdgeDist-PtDist;
                                    //if(event.node==tmp.Ni&&NodeDisttoQ[tmp.Ni]+PtDist>Dmax) Dmax=NodeDisttoQ[tmp.Ni]+PtDist;
                                    //if(event.node==tmp.Nj&&NodeDisttoQ[tmp.Nj]+EdgeDist-PtDist>Dmax) Dmax=NodeDisttoQ[tmp.Nj]+EdgeDist-PtDist;
                                    cS.push_back(tmp);
                                }
                            }

                        }
                    }
                }
                DStepEvent newevent=event;
                newevent.node=NewNodeID;
                newevent.dist=event.dist+EdgeDist;
                if (newevent.dist<NodeDisttoQ[NewNodeID])
                {
                    NodePre[NewNodeID]=event.node;
                    NodeDisttoQ[NewNodeID]=newevent.dist;
                    //if(NodeBCM_map.find(NewNodeID)!=NodeBCM_map.end()&&NodeBCM_map[NewNodeID].size()||newevent.dist<Dmax)//Has been changed on 05.29
                    sQ.push(newevent);

                }
                if(NodeBCM_map.find(NewNodeID)!=NodeBCM_map.end()&&NodeBCM_map[NewNodeID].size()==0)
                    NodeBCM_map.erase(NewNodeID);
            }
        }
        NodeBCM_map.erase(event.node);

    }
    //cout<<"POI Read           #:"<<poicnt<<endl;

    //cout<<"Road Node Tested   #:"<<roadnodecnt<<endl;
    //cout<<"Road Edge Expanded #:"<<edgeexpanded<<endl;
    base2numedgeexpand+=edgeexpanded;
    base2numnodeexpand+=roadnodecnt;

}

void Baselinemethod2(const QueryPoint &Q)
{
    InitQuery(Q,Querytype(1));
    INEBasedSpatialKeywordFilter(Q);
    RefineDist(cS,Q);
    CandidateSet::iterator iter=cS.begin();
    int PtVerified=0;
    while(iter!=cS.end())
    {
        PtVerified++;
        //Must refine dist of POI here!!
        /*
        float refineddist=(*iter).dist_Ni+NodeDisttoQ[(*iter).Ni]<(*iter).dist_Nj+NodeDisttoQ[(*iter).Nj]?(*iter).dist_Ni+NodeDisttoQ[(*iter).Ni]:(*iter).dist_Nj+NodeDisttoQ[(*iter).Nj];
        (*iter).dist_toquery=refineddist;
        */
        if(INEBasedVerify(*iter,(*iter).dist_toquery,Q.k))
        {
            rS.push_back(*iter);
        }
        iter++;
    }
    //cout<<"Total POI Verified #:"<<PtVerified<<endl;
}
void setSpecialZero(EnhancedCountMap &tempecm)
{
    auto iter=tempecm.begin();
    while(iter!=tempecm.end())
    {
        (*iter).second.second=0;
        iter++;
    }
}
void INEBasedEnhancedFilter(const QueryPoint &Q)
{

    DStepQueue sQ;
    //float Dmax=0.0;
    int poicnt=0;
    int roadnodecnt=0;
    int edgeexpanded=1;

    //Data
    int Ni=Q.Ni,Nj=Q.Nj;//Ni must less than Nj
    int AdjGrpAddr,AdjListSize,NewNodeID,PtGrpKey,PtNumOnEdge;
    unsigned long long keywords;
    float EdgeDist,PtDist;

    //Precess sub edge of Q.Ni,Q.Nj divided by Q
    AdjGrpAddr=getAdjListGrpAddr(Ni);
    getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);
    for (int i=0; i<AdjListSize; i++)
    {
        getVarE(ADJNODE_A,Ref(NewNodeID),AdjGrpAddr,i);
        //getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
        if(NewNodeID == Nj)
        {
            getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
            getVarE(PTKEY_A,Ref(PtGrpKey),AdjGrpAddr,i);
            //cout<<"("<<Ni<<","<<Nj<<") EdgeDist:"<<EdgeDist;
            if(PtGrpKey==-1)
            {
                //cout<<"No POI existed on Edge where Q located."<<endl;
            }
            else
            {
                getFixedF(SIZE_P,Ref(PtNumOnEdge),PtGrpKey);

                vector<POI> poinj;
                vector<POI> poini;

                //cout<<" PtNum on Edge where Q located:"<<PtNumOnEdge<<endl;
                //Notice the order POI visited on edge
                for(int j=0; j<PtNumOnEdge; j++)
                {
                    poicnt++;
                    getVarE(PT_P,Ref(PtDist),PtGrpKey,j);
                    getVarE(PT_VCT,Ref(keywords),PtGrpKey,j);


                    if(LcontianR(Q.keywords,keywords))
                    {
                        POI tmp;
                        tmp.Ni=Ni;
                        tmp.Nj=Nj;
                        tmp.dist_toquery=fabs(PtDist-Q.dist_Ni);
                        tmp.dist_Ni=PtDist;
                        tmp.dist_Nj=EdgeDist-PtDist;
                        tmp.keywords=keywords;
                        tmp.pre=-2;
                        //showCountMap(Ni,Querytype(2));
                        //if(keywords==16) cout<<" "<<endl;
                        if(PtDist<Q.dist_Ni&&NodeECM_map[Ni].find(keywords)!=NodeECM_map[Ni].end())
                        {


                            //updateEnhancedCountMap(Ni,keywords,false);
                            //showCountMap(Ni,Querytype(1));
                            //if(tmp.dist_toquery>Dmax) Dmax=tmp.dist_toquery;
                            //cS.push_back(tmp);
                            poini.push_back(tmp);
                        }

                        //showCountMap(Nj,Querytype(2));
                        else if(PtDist>Q.dist_Ni&&NodeECM_map[Nj].find(keywords)!=NodeECM_map[Nj].end())
                        {
                            poinj.push_back(tmp);

                            //updateEnhancedCountMap(Nj,keywords,false);

                            //showCountMap(Nj,Querytype(1));
                            //if(tmp.dist_toquery>Dmax) Dmax=tmp.dist_toquery;
                            //cS.push_back(tmp);
                        }
                        else {}

                    }
                }
                sort(poini.begin(),poini.end(),POIComparison());
                sort(poinj.begin(),poinj.end(),POIComparison());

                while(!poinj.empty())
                {
                    POI tmp=poinj.back();
                    if( NodeECM_map[Nj].find(tmp.keywords)!=NodeECM_map[Nj].end())
                    {
                        updateEnhancedCountMap(Nj,tmp.keywords,false);
                        cS.push_back(tmp);

                    }

                    poinj.pop_back();
                }
                while(!poini.empty())
                {
                    POI tmp=poini.back();
                    if( NodeECM_map[Ni].find(tmp.keywords)!=NodeECM_map[Ni].end())
                    {
                        updateEnhancedCountMap(Ni,tmp.keywords,false);
                        cS.push_back(tmp);

                    }

                    poini.pop_back();
                }
            }

            DStepEvent event;
            event.node=Ni;
            event.dist=Q.dist_Ni;

            if(NodeECM_map[Ni].size())
                sQ.push(event);//exist count != 0
            else
                NodeECM_map.erase(Ni);//2013.5.30
            NodeDisttoQ[Ni]=Q.dist_Ni;
            event.node=Nj;
            event.dist=EdgeDist-Q.dist_Ni;

            if(NodeECM_map[Nj].size())
                sQ.push(event);
            else
                NodeECM_map.erase(Nj);
            NodeDisttoQ[Nj]=event.dist;
            isVisted[Q.Ni]=true;//Remove Dual visit of edge where Q located
            isVisted[Q.Nj]=true;
            NodePre[Q.Ni]=-2;//-2 denote Q
            NodePre[Q.Nj]=-2;
            break;
        }
    }

    while(!sQ.empty())
    {
        if(NodeECM_map.size()==0) break;
        DStepEvent event=sQ.top();
        sQ.pop();


        if(isVisted[event.node]&&event.node!=Q.Ni&&event.node!=Q.Nj) continue;
        if(NodeECM_map.find(event.node)!=NodeECM_map.end())
            setSpecialZero(NodeECM_map[event.node]);
        //showCountMap(event.node,Querytype(2));
        isVisted[event.node]= true;
        roadnodecnt++;
        //for each non-visited adj. do
        AdjGrpAddr=getAdjListGrpAddr(event.node);
        getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);
        vector<int> adjnode;//For mutual effect
        adjnode.reserve(AdjListSize);
        for (int i=0; i<AdjListSize; i++)
        {
            getVarE(ADJNODE_A,Ref(NewNodeID),AdjGrpAddr,i);
            //getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
            if(!isVisted[NewNodeID])
            {
                edgeexpanded++;
                getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
                getVarE(PTKEY_A,Ref(PtGrpKey),AdjGrpAddr,i);

                if(NodeECM_map.find(event.node)!=NodeECM_map.end()&&NodeECM_map[event.node].size())//Notice here is event.node!!
                {
                    if(event.dist+EdgeDist<NodeDisttoQ[NewNodeID])//Only update countmap satisfy some constraint
                    {
                        EnhancedCountMap tmpecm=NodeECM_map[event.node];
                        //setSpecialZero(tmpecm);
                        NodeECM_map.insert(pair<int,EnhancedCountMap>(NewNodeID,tmpecm));
                        adjnode.push_back(NewNodeID);//Store all adj. has CountMap
                        if(PtGrpKey!=-1)
                        {
                            getFixedF(SIZE_P,Ref(PtNumOnEdge),PtGrpKey);
                            //cout<<" PtNumOnEdge:"<<PtNumOnEdge<<endl;
                            if(event.node<NewNodeID)//For purpose to visited in order
                            {
                                for(int j=PtNumOnEdge-1; j>=0; j--)
                                {
                                    poicnt++;
                                    getVarE(PT_P,Ref(PtDist),PtGrpKey,j);
                                    getVarE(PT_VCT,Ref(keywords),PtGrpKey,j);
                                    unsigned long long tmpkw=0;
                                    //condition of lemma3
                                    bool flag=false;
                                    if(PtDist<NodeDisttoQ[event.node]) flag=true;

                                    //for all poi who satisfy keywords constraint added it to CandidateSet
                                    if(LcontianR(Q.keywords,keywords)&&NodeECM_map[NewNodeID].find(keywords)!=NodeECM_map[NewNodeID].end())
                                    {
                                        updateEnhancedCountMap(NewNodeID,keywords,flag);
                                        POI tmp;
                                        tmp.keywords=keywords;
                                        tmp.Ni=event.node;
                                        tmp.Nj=NewNodeID;
                                        tmp.dist_Ni=PtDist;
                                        tmp.dist_Nj=EdgeDist-PtDist;
                                        //if(event.node==tmp.Ni&&NodeDisttoQ[tmp.Ni]+PtDist>Dmax) Dmax=NodeDisttoQ[tmp.Ni]+PtDist;
                                        //if(event.node==tmp.Nj&&NodeDisttoQ[tmp.Nj]+EdgeDist-PtDist>Dmax) Dmax=NodeDisttoQ[tmp.Nj]+EdgeDist-PtDist;
                                        cS.push_back(tmp);
                                    }
                                    else if(LintersectR(Q.keywords,keywords,tmpkw)&&NodeECM_map[NewNodeID].find(tmpkw)!=NodeECM_map[NewNodeID].end())
                                    {
                                        updateEnhancedCountMap(NewNodeID,tmpkw,flag);
                                    }
                                    else {}
                                }
                            }
                            else
                            {
                                for(int j=0; j<PtNumOnEdge; j++)
                                {
                                    poicnt++;
                                    getVarE(PT_P,Ref(PtDist),PtGrpKey,j);
                                    getVarE(PT_VCT,Ref(keywords),PtGrpKey,j);
                                    unsigned long long tmpkw=0;
                                    bool flag=false;
                                    if((EdgeDist-PtDist)<NodeDisttoQ[event.node]) flag=true;
                                    //for all poi who satisfy keywords constraint added it to CandidateSet
                                    if(LcontianR(Q.keywords,keywords)&&NodeECM_map[NewNodeID].find(keywords)!=NodeECM_map[NewNodeID].end())
                                    {
                                        updateEnhancedCountMap(NewNodeID,keywords,flag);
                                        POI tmp;
                                        tmp.keywords=keywords;
                                        tmp.Ni=NewNodeID;
                                        tmp.Nj=event.node;
                                        tmp.dist_Ni=PtDist;
                                        tmp.dist_Nj=EdgeDist-PtDist;
                                        //if(event.node==tmp.Ni&&NodeDisttoQ[tmp.Ni]+PtDist>Dmax) Dmax=NodeDisttoQ[tmp.Ni]+PtDist;
                                        //if(event.node==tmp.Nj&&NodeDisttoQ[tmp.Nj]+EdgeDist-PtDist>Dmax) Dmax=NodeDisttoQ[tmp.Nj]+EdgeDist-PtDist;
                                        cS.push_back(tmp);
                                    }
                                    else if(LintersectR(Q.keywords,keywords,tmpkw)&&NodeECM_map[NewNodeID].find(tmpkw)!=NodeECM_map[NewNodeID].end())
                                    {
                                        updateEnhancedCountMap(NewNodeID,tmpkw,flag);
                                    }
                                    else {}

                                }
                            }
                        }
                    }
                    else
                    {
                        if(PtGrpKey!=-1)
                        {
                            getFixedF(SIZE_P,Ref(PtNumOnEdge),PtGrpKey);
                            for(int j=0; j<PtNumOnEdge; j++)
                            {
                                poicnt++;
                                getVarE(PT_P,Ref(PtDist),PtGrpKey,j);
                                getVarE(PT_VCT,Ref(keywords),PtGrpKey,j);
                                //for all poi who satisfy keywords constraint added it to CandidateSet
                                if(LcontianR(Q.keywords,keywords))
                                {
                                    POI tmp;
                                    tmp.keywords=keywords;
                                    tmp.Ni=event.node<NewNodeID?event.node:NewNodeID;
                                    tmp.Nj=event.node<NewNodeID?NewNodeID:event.node;
                                    tmp.dist_Ni=PtDist;
                                    tmp.dist_Nj=EdgeDist-PtDist;
                                    //if(event.node==tmp.Ni&&NodeDisttoQ[tmp.Ni]+PtDist>Dmax) Dmax=NodeDisttoQ[tmp.Ni]+PtDist;
                                    //if(event.node==tmp.Nj&&NodeDisttoQ[tmp.Nj]+EdgeDist-PtDist>Dmax) Dmax=NodeDisttoQ[tmp.Nj]+EdgeDist-PtDist;
                                    cS.push_back(tmp);
                                }
                            }

                        }
                    }
                }
                DStepEvent newevent=event;
                newevent.node=NewNodeID;
                newevent.dist=event.dist+EdgeDist;

                if (newevent.dist<NodeDisttoQ[NewNodeID])
                {
                    NodePre[NewNodeID]=event.node;
                    NodeDisttoQ[NewNodeID]=newevent.dist;
                    ///////Problem find a correct terminate condition!!!
                    //if(NodeECM_map.find(NewNodeID)!=NodeECM_map.end()&&NodeECM_map[NewNodeID].size()||newevent.dist<Dmax)//NodeECM_map[NewNodeID].size()||
                    sQ.push(newevent);

                }

                if(NodeECM_map.find(NewNodeID)!=NodeECM_map.end()&&(NodeECM_map[NewNodeID].size()==0))
                {
                    NodeECM_map.erase(NewNodeID);
                }

            }
        }
        NodeECM_map.erase(event.node);
        //For mutual effect
        //May exist some false hit

        for(unsigned i=0; i<adjnode.size(); i++)
        {
            if(NodeECM_map.find(adjnode[i])!=NodeECM_map.end()&&(NodeECM_map[adjnode[i]].size()!=0))
                for(unsigned j=i+1; j<adjnode.size(); j++)
                    if(NodeECM_map.find(adjnode[j])!=NodeECM_map.end()&&(NodeECM_map[adjnode[j]].size()!=0) )
                    {
                        updateEnhancedCountMap(NodeECM_map[adjnode[i]],NodeECM_map[adjnode[j]]);
                        l3appliednum++;
                    }
        }



    }
    //cout<<"POI Read           #:"<<poicnt<<endl;
    //cout<<"Dmax               #:"<<Dmax<<endl;
    //cout<<"Road Node Tested   #:"<<roadnodecnt<<endl;
    //cout<<"Road Edge Expanded #:"<<edgeexpanded<<endl;
    enhancednumedgeexpand+=edgeexpanded;
    enhancednumnodeexpand+=roadnodecnt;
}

void RefineDist(CandidateSet &cS,const QueryPoint &Q)
{
    CandidateSet::iterator iter=cS.begin();
    while(iter!=cS.end())
    {
        if((!((*iter).Ni==Q.Ni&&(*iter).Nj==Q.Nj)))
            if((*iter).dist_Ni+NodeDisttoQ[(*iter).Ni]<(*iter).dist_Nj+NodeDisttoQ[(*iter).Nj])
            {
                //if((*iter).pre!=(*iter).Ni) cout<<"S";
                (*iter).pre=(*iter).Ni;
                (*iter).dist_toquery=(*iter).dist_Ni+NodeDisttoQ[(*iter).Ni];
            }
            else
            {
                //if((*iter).pre!=(*iter).Nj) cout<<"S";
                (*iter).pre=(*iter).Nj;
                (*iter).dist_toquery=(*iter).dist_Nj+NodeDisttoQ[(*iter).Nj];

            }
        else
        {
            (*iter).pre=-2;
        }
        iter++;
    }
}

bool preTest(POI &poi,int pre)
{
    int tmppre=poi.pre;
    while(tmppre!=pre&&tmppre!=-2)
    {
        tmppre=NodePre[tmppre];
    }
    if(tmppre==pre)
        return true;
    else
        return false;
}


void EnhancedMethod(const QueryPoint &Q)
{
    InitQuery(Q,Querytype(2));
    INEBasedEnhancedFilter(Q);
    RefineDist(cS,Q);
    CandidateSet::iterator iter=cS.begin();
    int PtVerified=0;
    int PtPruned=0;
    while(iter!=cS.end())
    {
        PtVerified++;
        POI tmp=(*iter);
        if(INEBasedVerify(tmp,tmp.dist_toquery,Q.k))
        {
            rS.push_back(tmp);
            iter++;
        }
        //else iter++;

        //Lemma 4 Remove POI Invalidated from CandidateSet
        else if(!((*iter).Ni==Q.Ni&&(*iter).Nj==Q.Nj))
        {
            l4appliednum++;
            iter++;

            while(iter!=cS.end()&&!((*iter).Ni==Q.Ni&&(*iter).Nj==Q.Nj))
            {

                if(tmp.dist_toquery<(*iter).dist_toquery&&LcontianR(tmp.keywords,(*iter).keywords)&&preTest(*iter,tmp.pre==tmp.Ni?tmp.Nj:tmp.Ni))
                {

                    PtPruned++;
                    iter++;
                }
                /*
                else if(tmp.dist_toquery<(*iter).dist_toquery&&LcontianR(tmp.keywords,(*iter).keywords)&&((*iter).Ni==tmp.Ni&&(*iter).Nj==tmp.Nj&&(*iter).pre==tmp.pre))
                {
                    PtPruned++;
                    iter++;
                }
                */

                else
                {
                    break;
                }

            }

        }

        else
        {
            iter++;
        }


    }
    l4prunnedperquery=PtPruned;
    //cout<<"Total POI Verified #:"<<PtVerified<<endl;
    //cout<<"Lemma3 successfully applied #:"<<l3prunnedperquery<<endl;
    //cout<<"POI Pruned by Lemma4:"<<PtPruned<<endl;
}

bool isPOIEqual(POI p0,POI p1)
{
    if((p0.Ni==p1.Ni)&&(p0.Nj==p1.Nj)&&fabs(p0.dist_Ni-p1.dist_Ni)<1e-6)
    {
        return true;
    }
    else
        return false;
}


//For each POI in set rS0 find if it exist in set rS1
//Put poi do not exist in rS1 to tmpcS
int compareResult(ResultSet rS0,ResultSet rS1,CandidateSet &tmpcS)
{
    int r=0;
    for(POI p:rS0)
    {
        bool flag=false;
        auto iter=rS1.begin();
        while(iter!=rS1.end())
        {
            if(isPOIEqual(p,(*iter)))
            {
                flag=true;
                break;
            }
            iter++;
        }
        if(!flag)
        {
            r++;
            tmpcS.push_back(p);
        }
    }
    return r;
}

void compareNodedist(vector<float> &v1,vector<float> &v2)
{
    for(int i=0; i<NodeNum; i++)
        if(fabs(v1[i]-v2[i])>1e-6&&fabs(v1[i]-MAX_DIST)>1e-6)
        {
            cout<<i<<":"<<v1[i]<<" "<<v2[i]<<endl;
        }

}

void WriteQuery2File(vector<QueryPoint> Qset,const char* filename)
{
    char tmpfilename[255];
    sprintf(tmpfilename,"%s.txt",filename);
    remove(tmpfilename);
    FILE *f=fopen(tmpfilename,"w");
    for(QueryPoint Q:Qset)
    {
        char stringline[255];
        sprintf(stringline,"%d %d %f %llu %d\n",Q.Ni,Q.Nj,Q.dist_Ni,Q.keywords,Q.k);
        fputs(stringline,f);
    }
    fclose(f);
}

bool ReadQueryFromFile(const char* filename,vector<QueryPoint> &Qset)
{
    char tmpfilename[255];
    sprintf(tmpfilename,"%s.txt",filename);
    FILE *f=fopen(tmpfilename,"r");
    if(!f)
    {
        cout<<"QueryPoint file read error!"<<endl;
        return false;
    }
    QueryPoint Q;
    while(fscanf(f,"%d %d %f %llu %d",&Q.Ni,&Q.Nj,&Q.dist_Ni,&Q.keywords,&Q.k)!=EOF)
    {
        Qset.push_back(Q);
    }
    fclose(f);
    return true;
}


int main(int argc,char** argv)
{
    //    Query File default to map_queryPoints_querykeywordsnumber_k_cachepages
    string configFileName = "config.prop";
    ConfigType cr(configFileName,argc, argv);
    cr.ListConfig();
    InitClock();
    OpenDiskComm(cr.getDataFileName().c_str(),cr.getParameterCachePages());

    int cntnotmached=0;
    int implcnt=0;
    vector<QueryPoint> Qset;
    vector<QueryPoint> QnotMatched;
    vector<QueryPoint> QresultNotZero;
    std::string queryfile = cr.getQueryFileName();
    
    if(!ReadQueryFromFile(queryfile.c_str(),Qset))
    {
        return -1;
    }


    for(auto Q:Qset)
    {
        
        implcnt++;//Should write to file
        //Baselinemethod1
        cout<<"///////////////////////////////////////////////////////////////////////"<<endl;
        cout<<"Baseline method 1:"<<endl;
        unsigned long time0=clock();
        Baselinemethod1(Q);
        unsigned long time0e=clock();
        if(time0e<time0)
            time0=LONG_MAX*2+1-time0+time0e;
        else
            time0=time0e-time0;
        totalbase1+=time0;
        int resultnum0=rS.size();
        int candidatenum0=cS.size();
        base1poiverified+=cS.size();
        //cout<<"Result #:"<<resultnum0<<endl;
        //printQueryResult();
        base1datapageaccessed+=printPageAccess();
        cout<<"QueryTime:"<<time0<<endl;
        //Baselinemethod2
        cout<<"///////////////////////////////////////////////////////////////////////"<<endl;
        cout<<"Baseline method 2:"<<endl;
        unsigned long time1=clock();
        Baselinemethod2(Q);
        unsigned long time1e=clock();

        if(time1e<time1)
            time1=LONG_MAX*2+1-time1+time1e;
        else
            time1=time1e-time1;
        totalbase2+=time1;
        int resultnum1=rS.size();
        int candidatenum1=cS.size();
        base2poiverified+=cS.size();
        l1prunned+=candidatenum0-candidatenum1;
        //cout<<"Result #:"<<resultnum1<<endl;

        base2datapageaccessed+=printPageAccess();
        cout<<"QueryTime:"<<time1<<endl;
        //Enhanced method
        cout<<"///////////////////////////////////////////////////////////////////////"<<endl;
        cout<<"Enhanced method :"<<endl;
        unsigned long time2=clock();
        EnhancedMethod(Q);
        unsigned long time2e=clock();

        if(time2e<time2)
            time2=LONG_MAX*2+1-time2+time2e;
        else
            time2=time2e-time2;

        totalenhanced+=time2;
        int resultnum2=rS.size();
        int enhancedcandidatenum=cS.size();
        //cout<<"Result #:"<<resultnum2<<endl;
        enhancedpoiverified+=cS.size()-l4prunnedperquery;
        l2prunned+=(candidatenum1-l3prunnedperquery-l4prunnedperquery);
        l3prunned+=l3prunnedperquery;
        l4prunned+=l4prunnedperquery;
        enhanceddatapageaccessed+=printPageAccess();

        if(rS.size()>0)
            QresultNotZero.push_back(Q);
        cout<<"QueryTime:"<<time2<<endl;
         cout<<"///////////////////////////////////////////////////////////////////////"<<endl;
        cout<<"Baseline1 Clock Elapsed:"<<time0<<endl;
        cout<<"Baseline2 Clock Elapsed:"<<time1<<endl;
        cout<<"Enhanced  Clock Elapsed:"<<time2<<endl;
        //printQueryResult();
        cout<<"///////////////////////////////////////////////////////////////////////"<<endl;

        //cout<<"***********************************************************************"<<endl;
        //cout<<endl;
    }

    cout<<"Total Query Implemented:"<<implcnt<<endl;
    
    if(QresultNotZero.size()||QnotMatched.size())
    {
        string querynotzero=queryfile+"NotZero";
        string querynotmatched=queryfile+"NotMatched";
        if(QresultNotZero.size())WriteQuery2File(QresultNotZero,querynotzero.c_str());
        if(QnotMatched.size())WriteQuery2File(QnotMatched,querynotmatched.c_str());
    }

    string queryresult= cr.getQueryResultFileName();
    ofstream fout(queryresult);
    fout<<"QueryName:"<<queryfile<<endl;
    fout<<"Query Number Implemented:"<<implcnt<<endl;
    fout<<"Result not Zero:"<<QresultNotZero.size()<<endl;
    fout<<"Result not Matched:"<<QnotMatched.size()<<endl;
    fout<<"AvgQueryTime(msec):"<<endl;
    fout<<"     Baseline1:"<<totalbase1/float(implcnt)<<endl;
    fout<<"     Baseline2:"<<totalbase2/float(implcnt)<<endl;
    fout<<"     EnhancedM:"<<totalenhanced/float(implcnt)<<endl;
    fout<<"AvgPOINumVerified:"<<endl;
    fout<<"     Baseline1:"<<base1poiverified/float(implcnt)<<endl;
    fout<<"     Baseline2:"<<base2poiverified/float(implcnt)<<endl;
    fout<<"     EnhancedM:"<<enhancedpoiverified/float(implcnt)<<endl;
    fout<<"AvgPageAccessed:"<<endl;
    fout<<"     Baseline1:"<<float(base1datapageaccessed)/implcnt<<endl;
    fout<<"     Baseline2:"<<float(base2datapageaccessed)/implcnt<<endl;
    fout<<"     EnhancedM:"<<float(enhanceddatapageaccessed)/implcnt<<endl;
    fout<<"AvgNumPrunedbyEachLemma:"<<endl;
    fout<<"     Lemma1:"<<float(l1prunned)/implcnt<<endl;
    fout<<"     Lemma2:"<<float(l2prunned)/implcnt<<endl;
    fout<<"     Lemma3:"<<float(l3prunned)/implcnt<<endl;
    fout<<"     Lemma4:"<<float(l4prunned)/implcnt<<endl;
    fout<<"AvgNumEachLemmaApplied:"<<endl;
    fout<<"     Lemma3:"<<float(l3appliednum)/implcnt<<endl;
    fout<<"     Lemma4:"<<float(l4appliednum)/implcnt<<endl;
    fout<<"AvgNumNodeExpanded:"<<endl;
    fout<<"     Baseline2:"<<float(base2numnodeexpand)/implcnt<<endl;
    fout<<"     EnhancedM:"<<float(enhancednumnodeexpand)/implcnt<<endl;
    fout<<"AvgNumEdgeExpanded:"<<endl;
    fout<<"     Baseline2:"<<float(base2numedgeexpand)/implcnt<<endl;
    fout<<"     EnhancedM:"<<float(enhancednumedgeexpand)/implcnt<<endl;
    fout<<endl;
    fout.close();
    CloseDiskComm();
    PrintElapsed();

    return 0;
}


