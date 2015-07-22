#include <iostream>
#include <unordered_map>
#include <queue>
#include <bitset>
#include <algorithm>
#include "diskbased.h"
#include "ConfigType.h"
#include "KeywordsGenerator.h"
#include "egtree.h"
#include <fstream>

using namespace std;
// all paremeters
#define CandidateSet vector<POI>
#define ResultSet vector<int>
#define EA 1
#define EGBU 2
#define EGTD 3
#define INFINITE_MAX 100000000.0
#define visited 4
#define hVisited 5
#define unVisited 6
//#define halfVisited vector<int>
int nOfDomiTest;
int nOfEdgeExpd;

CandidateSet cS;
ResultSet rS;

// for dijkstra
//halfVQueue hvQ;
map<edgePair, edgeState, eSComparison> visitedState ;
map<int,float> distTQ;
dVQueue dvq;

// for EGBU
struct PartAddr {
	int part;
	int addr;
};
map<int,PartAddr> paID;
vector<TreeNode> EGT;
set<int> visitedSet; //record the  id of visited treeNodeID

// to verify whether R contains by L
bool LcontianRIKwd(const set<int> L,const set<int> R)
{
    //High bit set to zero
    //R&=tmp;
	
	if(L.size() < R.size()) return false;
	set<int>::iterator it;
	for(it = R.begin(); it != R.end(); it++) {
		if (find(L.begin(), L.end(), (*it)) != L.end()) {
			
		} else {
			return false;
		}
	}
    return  true;
	
	//set<int> temp;
	//set_intersection(L.begin(), L.end(), R.begin(), R.end(), temp.begin());
}
// smaller is better
bool isBeDominate(POI poi, CandidateSet cS, const QueryPoint &Q) {
	for(int i=0; i<cS.size(); ) {
		POI temp = cS[i];
		nOfDomiTest++;
		bool pDt = true;
		bool tDp = true;
		int pEt = 0;
		for(int j=0; j<ATTRIBUTE_DIMENSION; j++) {
			if(Q.subSpace[j]==0) {
			} else {
				if(poi.attr[j]<temp.attr[j]) {
					tDp = false;
				} else if(poi.attr[j] == temp.attr[j]) {
					pEt++;
				} else {
					pDt = false;
				}
			}
		}
		// poi dominate temp
		if(pDt&&pEt!=Q.subSpace.count()) {
			cS.erase(cS.begin()+i);
		}
		if(tDp&&pEt!=Q.subSpace.count()) {
			return true;
		}
		if(pEt == Q.subSpace.count()) {
			i++;
		}
		//nOfDomiTest++;
	}
	return false;
}

void initialQuery(int algorithmId, FILE *ptAddr) {
	nOfDomiTest = 0;
	nOfEdgeExpd = 0;
	int nOfNode;
	int nodeID;
	if(algorithmId == EA) { // EA

	} else if (algorithmId == EGBU){ // EGBU
		// load partAddr
		fread( &nOfNode, sizeof(int), 1, ptAddr );
		for(int i=0; i<nOfNode; i++) {
			PartAddr pa;
			fread( &nodeID, sizeof(int), 1, ptAddr );
			fread( &pa.part, sizeof(int), 1, ptAddr );
			fread( &pa.addr, sizeof(int), 1, ptAddr );
			paID[nodeID] = pa;
		}
		// load egtree
		egtree_load(EGT);

	} else { // EGTD
		// load partAddr
		fread( &nOfNode, sizeof(int), 1, ptAddr );
		for(int i=0; i<nOfNode; i++) {
			PartAddr pa;
			fread( &nodeID, sizeof(int), 1, ptAddr );
			fread( &pa.part, sizeof(int), 1, ptAddr );
			fread( &pa.addr, sizeof(int), 1, ptAddr );
			paID[nodeID] = pa;
		}
		// load egtree
		egtree_load(EGT);
	}
}

void printResult() {
	for(int i=0; i<rS.size(); i++) {
		cout<<rS[i]<<" ";
	}
}

void EABasedAlgorithm(const QueryPoint &Q) {
	DStepQueue sQ;
    int poicnt=0;
    int roadnodecnt=0;

    //Data
    int Ni=Q.Ni,Nj=Q.Nj;//Ni must less than Nj
    int AdjGrpAddr,AdjListSize,NewNodeID,PtGrpKey,PtNumOnEdge;
    //unsigned long long keywords;
	set<int> tempKwds;
	int *tempKwd;
    float EdgeDist,PtDist;
    //Precess sub edge of Q.Ni,Q.Nj divided by Q
    AdjGrpAddr=getAdjListGrpAddr(Ni);
    getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);
    for (int i=0; i<AdjListSize; i++){
        getVarE(ADJNODE_A,Ref(NewNodeID),AdjGrpAddr,i);
        //getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
        if(NewNodeID == Nj){
            getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
            getVarE(PTKEY_A,Ref(PtGrpKey),AdjGrpAddr,i);

            nOfEdgeExpd++;

            //cout<<"("<<Ni<<","<<Nj<<") EdgeDist:"<<EdgeDist;
            if(PtGrpKey==-1){
                //cout<<"No POI existed on Edge where Q located."<<endl;
            }
            else {
                getFixedF(SIZE_P,Ref(PtNumOnEdge),PtGrpKey);
                //cout<<" PtNum on Edge where Q located:"<<PtNumOnEdge<<endl;
                //Notice the order POI visited on edge
                for(int j=0; j<PtNumOnEdge; j++){
                    poicnt++;
                    getVarE(PT_DIST,Ref(PtDist),PtGrpKey,j);

                    getVarE(PT_VCT,tempKwd,PtGrpKey,j);
					int num = sizeof(tempKwd)/sizeof(int);
					for(int loop=0; loop<sizeof(tempKwd); loop+sizeof(int)){
						int temp;
						memcpy(&temp,tempKwd+loop,sizeof(int));
						tempKwds.insert(temp);
					}

					if(LcontianRIKwd(tempKwds,Q.kwd)){
                        POI tmp;
						tmp.dist_toquery = fabs(PtDist-Q.dist_Ni);
						getVarE(PT_ATTRIBUTE,tmp.attr,PtGrpKey,j);
						getVarE(PT_P,&tmp.poid,PtGrpKey,j);
						
						if(tmp.dist_toquery <= Q.distCnst) {
							if(isBeDominate(tmp, cS, Q)) {
							} else {
								cS.push_back(tmp);
							}
						}             
                    }// endif
                }
                
            }// endelse
        }
    }
	// record the visited state of edge
    edgePair epi,epj;
	edgeState esi,esj;
	epi.Ni = 0;
	epi.Nj = Q.Ni;
	esi.vState = visited;
	esi.iDisToQuery = 0;
	esi.jDisToQuery = Q.dist_Ni;

	epj.Ni = 0;
	epj.Nj = Q.Nj;
	esj.iDisToQuery = 0;
	esj.jDisToQuery = EdgeDist-Q.dist_Ni;

	visitedState[epi] = esi;
	visitedState[epj] = esj;


	// record the visited node information
	dijkVisit dvi,dvj;
	dvi.N = Q.Ni;
	dvi.disTQ = Q.dist_Ni;

	dvj.N = Q.Nj;
	dvj.disTQ = EdgeDist-Q.dist_Ni;
	
	dvq.push(dvi);
	dvq.push(dvj);

	while(!dvq.empty()) {
		dijkVisit tmpDV = dvq.top();
		dvq.pop();
		if(tmpDV.disTQ <= Q.distCnst) {
			AdjGrpAddr=getAdjListGrpAddr(tmpDV.N);
			getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);
			for (int i=0; i<AdjListSize; i++){
				getVarE(ADJNODE_A,Ref(NewNodeID),AdjGrpAddr,i);
				getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
				getVarE(PTKEY_A,Ref(PtGrpKey),AdjGrpAddr,i);
				edgePair ep;
				ep.Ni = NewNodeID<tmpDV.N?NewNodeID:tmpDV.N;
				ep.Nj = NewNodeID>tmpDV.N?NewNodeID:tmpDV.N;
				map<edgePair, edgeState, eSComparison>::iterator iter;
				iter = visitedState.find(ep);
				if(iter!=visitedState.end()) {
					if(visitedState[ep].vState == hVisited) {// 计算两端的POI，修改状态
						visitedState[ep].vState == visited;
						//---添加距离visitedState
						if(tmpDV.N < NewNodeID) {
							visitedState[ep].iDisToQuery = tmpDV.disTQ;
						} else {
							visitedState[ep].jDisToQuery = tmpDV.disTQ;
						}
						//从两端加入POIs
						if(PtGrpKey==-1){
						//cout<<"No POI existed on Edge where Q located."<<endl;
						} else {

							getFixedF(SIZE_P,Ref(PtNumOnEdge),PtGrpKey);
							//cout<<" PtNum on Edge where Q located:"<<PtNumOnEdge<<endl;
							//Notice the order POI visited on edge
							for(int j=0; j<PtNumOnEdge; j++){
								poicnt++;
								getVarE(PT_DIST,Ref(PtDist),PtGrpKey,j);

								getVarE(PT_VCT,tempKwd,PtGrpKey,j);
								int num = sizeof(tempKwd)/sizeof(int);
								for(int loop=0; loop<sizeof(tempKwd); loop+sizeof(int)){
									int temp;
									memcpy(&temp,tempKwd+loop,sizeof(int));
									tempKwds.insert(temp);
								}

								if(LcontianRIKwd(tempKwds,Q.kwd)){
									POI tmp;
									//tmp.dist_toquery = fabs(PtDist-Q.dist_Ni);
									float d1,d2;
									d1 = visitedState[ep].iDisToQuery + PtDist;
									d2 = visitedState[ep].jDisToQuery + EdgeDist - PtDist;
									tmp.dist_toquery = d1<d2 ? d1:d2;
									
									getVarE(PT_ATTRIBUTE,tmp.attr,PtGrpKey,j);
									getVarE(PT_P,&tmp.poid,PtGrpKey,j);
						
									if(tmp.dist_toquery <= Q.distCnst) {
										if(isBeDominate(tmp, cS, Q)) {
										} else {
											cS.push_back(tmp);
										}
									}             
								}// endif
							}
						}
					} else if(visitedState[ep].vState == visited) { //修改距离
						//---添加距离visitedState
						if(tmpDV.N < NewNodeID) {
							visitedState[ep].iDisToQuery = tmpDV.disTQ;
						} else {
							visitedState[ep].jDisToQuery = tmpDV.disTQ;
						}
					} else {

					}
				} else {
					edgePair ep;
					edgeState es;
					if(tmpDV.N < NewNodeID) {
						ep.Ni = tmpDV.N;
						ep.Nj = NewNodeID;
						es.iDisToQuery = tmpDV.disTQ;
						es.jDisToQuery = INFINITE_MAX;
					} else {
						ep.Ni = NewNodeID;
						ep.Nj = tmpDV.N;
						es.iDisToQuery = INFINITE_MAX;
						es.jDisToQuery = tmpDV.disTQ;
					}
					if(distTQ.find(NewNodeID)!=distTQ.end()){
						if(distTQ[NewNodeID]>(tmpDV.disTQ+EdgeDist)) {
							distTQ[NewNodeID] = tmpDV.disTQ+EdgeDist;
						}
					} else {
						distTQ[NewNodeID] = tmpDV.disTQ+EdgeDist;
					}
					if(distTQ[NewNodeID]<=Q.distCnst) {
						//将整条边加入
						es.vState = visited;
						visitedState[ep] = es;

						if(PtGrpKey==-1){
						//cout<<"No POI existed on Edge where Q located."<<endl;
						} else {
							getFixedF(SIZE_P,Ref(PtNumOnEdge),PtGrpKey);
							//cout<<" PtNum on Edge where Q located:"<<PtNumOnEdge<<endl;
							//Notice the order POI visited on edge
							for(int j=0; j<PtNumOnEdge; j++){
								poicnt++;
								getVarE(PT_DIST,Ref(PtDist),PtGrpKey,j);

								getVarE(PT_VCT,tempKwd,PtGrpKey,j);
								int num = sizeof(tempKwd)/sizeof(int);
								for(int loop=0; loop<sizeof(tempKwd); loop+sizeof(int)){
									int temp;
									memcpy(&temp,tempKwd+loop,sizeof(int));
									tempKwds.insert(temp);
								}

								if(LcontianRIKwd(tempKwds,Q.kwd)){
									POI tmp; //这时候可能不是真实距离,也不会影响结果
									if(tmpDV.N<NewNodeID) {
										tmp.dist_toquery = tmpDV.disTQ+PtDist;
									} else {
										tmp.dist_toquery = tmpDV.disTQ+EdgeDist-PtDist;
									}
									
								
									getVarE(PT_ATTRIBUTE,tmp.attr,PtGrpKey,j);
									getVarE(PT_P,&tmp.poid,PtGrpKey,j);
						
									if(tmp.dist_toquery <= Q.distCnst) {
										if(isBeDominate(tmp, cS, Q)) {
										} else {
											cS.push_back(tmp);
										}
									}             
								}// endif
							}
						}

					} else {// 如果只有部分在边上，hVisited
						es.vState = hVisited;
						visitedState[ep] = es;
					}
				}

				//getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
				//getVarE(PTKEY_A,Ref(PtGrpKey),AdjGrpAddr,i);
				if(PtGrpKey==-1){
						//cout<<"No POI existed on Edge where Q located."<<endl;
				} else {

					getFixedF(SIZE_P,Ref(PtNumOnEdge),PtGrpKey);
					//cout<<" PtNum on Edge where Q located:"<<PtNumOnEdge<<endl;
					//Notice the order POI visited on edge
					for(int j=0; j<PtNumOnEdge; j++){
						poicnt++;
						getVarE(PT_DIST,Ref(PtDist),PtGrpKey,j);

						getVarE(PT_VCT,tempKwd,PtGrpKey,j);
						int num = sizeof(tempKwd)/sizeof(int);
						for(int loop=0; loop<sizeof(tempKwd); loop+sizeof(int)){
							int temp;
							memcpy(&temp,tempKwd+loop,sizeof(int));
							tempKwds.insert(temp);
						}

						if(LcontianRIKwd(tempKwds,Q.kwd)){
							POI tmp;
							tmp.dist_toquery = fabs(PtDist-Q.dist_Ni);
							getVarE(PT_ATTRIBUTE,tmp.attr,PtGrpKey,j);
							getVarE(PT_P,&tmp.poid,PtGrpKey,j);
						
							if(tmp.dist_toquery <= Q.distCnst) {
								if(isBeDominate(tmp, cS, Q)) {
								} else {
									cS.push_back(tmp);
								}
							}             
						}// endif
					}
				}

			
			}// end for
		}// endif
	}// end while
    
	// 处理hVisited

	map<edgePair, edgeState, eSComparison>::iterator iterhV;
	for(iterhV=visitedState.begin(); iterhV!=visitedState.end(); ++iterhV) {
		if(iterhV->second.vState == hVisited) {//处理只有距离当前节点满足距离约束的POI
			int nodei,nodej;
			if(iterhV->second.iDisToQuery == INFINITE_MAX) {
				nodei = iterhV->first.Nj;
				nodej = iterhV->first.Ni;
			} else {
				nodei = iterhV->first.Nj;
				nodej = iterhV->first.Ni;
			}
			//获取边上的POI
			AdjGrpAddr=getAdjListGrpAddr(nodei);
			getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);
			for (int i=0; i<AdjListSize; i++){
				getVarE(ADJNODE_A,Ref(NewNodeID),AdjGrpAddr,i);
				//getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
				if(NewNodeID == nodej){
					getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
					getVarE(PTKEY_A,Ref(PtGrpKey),AdjGrpAddr,i);

					nOfEdgeExpd++;

					//cout<<"("<<Ni<<","<<Nj<<") EdgeDist:"<<EdgeDist;
					if(PtGrpKey==-1){
						//cout<<"No POI existed on Edge where Q located."<<endl;
					}
					else {
						getFixedF(SIZE_P,Ref(PtNumOnEdge),PtGrpKey);
						//cout<<" PtNum on Edge where Q located:"<<PtNumOnEdge<<endl;
						//Notice the order POI visited on edge
						for(int j=0; j<PtNumOnEdge; j++){
							poicnt++;
							getVarE(PT_DIST,Ref(PtDist),PtGrpKey,j);

							getVarE(PT_VCT,tempKwd,PtGrpKey,j);
							int num = sizeof(tempKwd)/sizeof(int);
							for(int loop=0; loop<sizeof(tempKwd); loop+sizeof(int)){
								int temp;
								memcpy(&temp,tempKwd+loop,sizeof(int));
								tempKwds.push_back(temp);
							}

							if(LcontianRIKwd(tempKwds,Q.kwd)){
								POI tmp;
								if(nodei<nodej){
									tmp.dist_toquery = iterhV->second.iDisToQuery+PtDist;
								} else {
									tmp.dist_toquery = iterhV->second.jDisToQuery+EdgeDist-PtDist;
								}
								getVarE(PT_ATTRIBUTE,tmp.attr,PtGrpKey,j);
								getVarE(PT_P,&tmp.poid,PtGrpKey,j);
						
								if(tmp.dist_toquery <= Q.distCnst) {
									if(isBeDominate(tmp, cS, Q)) {
									} else {
										cS.push_back(tmp);
									}
								}             
							}// endif
						}
                
					}// endelse
				}
			}
		}
	}
	

	for(int i=0; i<cS.size(); i++) {
		rS.push_back(cS[i].poid);
		//cout<<cS[i].poid<<" ";
	}
	//cout<<endl;
    //cout<<"POI Read           #:"<<poicnt<<endl;  
    //cout<<"Road Dominate Tested   #:"<<nOfDomiTest<<endl;
    //cout<<"Road Edge Expanded #:"<<nOfEdgeExpd<<endl;
    //base2numedgeexpand+=edgeexpanded;
    //base2numnodeexpand+=roadnodecnt;
}

void EGBUAlgorithm(const QueryPoint &Q) {
	//process leadnode, two situation:1, not in same part 2, in the same part
	int Ni=Q.Ni,Nj=Q.Nj;//Ni must less than Nj
    int AdjGrpAddr,AdjListSize,NewNodeID,PtGrpKey,PtNumOnEdge;
    //unsigned long long keywords;
	vector<int> tempKwds;
	int *tempKwd;
    float EdgeDist,PtDist;
	// locate the treeNodeID
	int trNodeIDi,trNodeIDj;
	int currTNID;
	float minDist;

	trNodeIDi = paID[Ni].part;
	trNodeIDj = paID[Nj].part;
	if(trNodeIDi == trNodeIDj) { //两个在一个划分块里面
		currTNID = trNodeIDi;

	} else { //查询点在边界边上面

	}


    //Precess sub edge of Q.Ni,Q.Nj divided by Q
    AdjGrpAddr=getAdjListGrpAddr(Ni);
    getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);
}

void EGTDAlgorithm(const QueryPoint &Q) {

}

void queryAlgorithm(int algorithmId, const QueryPoint &Q) {
	initialQuery(algorithmId);
	if(algorithmId == EA) { // EA
		EABasedAlgorithm(Q);
	} else if (algorithmId == EGBU){ // EGBU
		EGBUAlgorithm(Q);
	} else { // EGTD
		EGTDAlgorithm(Q);
	}
}

int main() {

	return 0;

}