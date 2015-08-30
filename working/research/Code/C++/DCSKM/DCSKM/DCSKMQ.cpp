#include <iostream>
#include <unordered_map>
#include <queue>
#include <bitset>
#include <algorithm>
#include <time.h>
#include "diskbased.h"
#include "ConfigType.h"
#include "KeywordsGenerator.h"
#include "egtree.h"
#include <fstream>
#include "DCSKMQ.h"

using namespace std;

// global paremeters
#define CandidateSet vector<POI>
#define ResultSet vector<int> // record the id of poi

#define EA 1
#define EGBU 2
#define EGTD 3
#define INFINITE_MAX 100000000.0
#define visited 4
#define hVisited 5
#define unVisited 6

// -----------performance parameters----------------
int nOfDominateTest;
int nOfEdgeExpended; // whether the pois of this edge are visited
int nOfPOIVisited; 
clock_t start, end;

// -----------algorihtms execute parameters------------
int algorithmId;
float queryEdgeDist;
CandidateSet cS;
ResultSet rS;

QueryPoint Q;

// for dijkstra
// halfVQueue hvQ;
//map<edgePair, edgeState, eSComparison> edgeStates;
map<int, edgeState> edgeStates; // <edge, edgeState>
map<int, float> distTQ; // <vertex, distance>
dVQueue dvq; 
set<int> visitedVtx;

// for EGBU
struct PartAddr {
	int part;
	int addr;
};

map<int, PartAddr> paID; // <vertex, address&partID>
vector<TreeNode> EGT; // load the egtree index
//set<int> visitedTreeNode; //record the  id of visited treeNodeID

struct TreeNodeC
{
	bool operator () (const TreeNode& left, const TreeNode& right) const
	{
		return left.minDistTQ > right.minDistTQ;
	}
};

typedef	priority_queue<TreeNode, vector<TreeNode>, TreeNodeC> TreeNodeQ;
TreeNodeQ tnq;
vector<int> pathID; // record the path from query to upmost treenode id

// for EGTD
set<int> visitedTNSet; // record the visited tree nodes

struct CandTDValue {
	locSkyline ls;
	set<POI> poi;
};

//自定义排序函数，要求是递增函数
bool SortByQuerySum(const CandTDValue &v1, const CandTDValue &v2)//注意：本函数的参数的类型一定要与vector中元素的类型一致  
{
	int sum1 = 0, sum2 = 0;
	for (int i = 0; i < Q.subSpace.size(); i++) {
		if (Q.subSpace[i] == 1) {
			sum1 += v1.ls.blockID[i];
			sum2 += v2.ls.blockID[i];
		}
	}
	return sum1 > sum2;
	//return v1.kwdCombination.size() > v2.kwdCombination.size();//降序排列  
}
map<int, CandTDValue> cSTD; // record the candidate POI, which is organized by the block id <blockID, ls&pois>

// to verify whether R contains by L
bool LcontianRIKwd(const set<int> L, const set<int> R)
{
	//High bit set to zero
	//R&=tmp;
	/*
	if (L.size() < R.size()) return false;
	set<int>::iterator it;
	for (it = R.begin(); it != R.end(); it++) {
		if (find(L.begin(), L.end(), (*it)) != L.end()) {

		}
		else {
			return false;
		}
	}
	return  true;
	*/
	set<int> temp;
	set_intersection(L.begin(), L.end(), R.begin(), R.end(), temp.begin());
	if (temp.size() == R.size()) return true;
	return false;
}
// smaller is better
void updatecSDominate(POI poi, CandidateSet &cS) {
	vector<POI>::iterator it = cS.begin();
	for (; it != cS.end(); ) {
		POI temp = *it;
		nOfDominateTest++;
		bool pDt = true;
		bool tDp = true;
		int pEt = 0;
		for (int j = 0; j<ATTRIBUTE_DIMENSION; j++) {
			if (Q.subSpace[j] == 0) {
			}
			else {
				if (poi.attr[j]<temp.attr[j]) {
					tDp = false;
				}
				else if (poi.attr[j] == temp.attr[j]) {
					pEt++;
				}
				else {
					pDt = false;
				}
			}
		}
		// poi dominate temp
		if (pDt&&pEt != Q.subSpace.count()) {
			it  = cS.erase(it);
			continue;
		}
		if (tDp&&pEt != Q.subSpace.count()) {
			return ;
		}
		if (!pDt&&!tDp) { // 互不支配
			it++;
			continue;
		}
		if (pEt == Q.subSpace.count()) {
			cS.push_back(poi);
			return ;
		}
		//nOfDominateTest++;
	}
	if (it == cS.end()) {
		cS.push_back(poi);
	}
}

bool isBeBoundDominate(float attrBound[ATTRIBUTE_DIMENSION][2], CandidateSet cS, const QueryPoint &Q) {
	for (int i = 0; i<cS.size(); ) {
		POI temp = cS[i];
		nOfDominateTest++;
		bool bDt = true;
		bool tDb = true;
		int pEt = 0;
		for (int j = 0; j<ATTRIBUTE_DIMENSION; j++) {
			if (Q.subSpace[j] == 0) {
			}
			else {
				if (attrBound[j][0]<temp.attr[j]) {
					tDb = false;
				}
				else if (attrBound[j][0] == temp.attr[j]) {
					pEt++;
				}
				else {
					bDt = false;
				}
			}
		}
		// poi dominate temp
		if (tDb&&pEt != Q.subSpace.count()) {
			return true;
		}
		if (pEt == Q.subSpace.count()) {
			i++;
		}
		//nOfDominateTest++;
	}
	return false;
}


bool isBeLocSkyDominate(float attrBound[ATTRIBUTE_DIMENSION][2], CandidateSet cS, const QueryPoint &Q) {
	for (int i = 0; i<cS.size(); ) {
		POI temp = cS[i];
		nOfDominateTest++;
		bool bDt = true;
		bool tDb = true;
		int pEt = 0;
		for (int j = 0; j<ATTRIBUTE_DIMENSION; j++) {
			if (Q.subSpace[j] == 0) {
			}
			else {
				if (attrBound[j][0]<temp.attr[j]) {
					tDb = false;
				}
				else if (attrBound[j][0] == temp.attr[j]) {
					pEt++;
				}
				else {
					bDt = false;
				}
			}
		}
		// poi dominate temp
		if (tDb&&pEt != Q.subSpace.count()) {
			return true;
		}
		if (pEt == Q.subSpace.count()) {
			i++;
		}
		//nOfDominateTest++;
	}
	return false;
}

void partAddrLoad(const char* filename, map<int, PartAddr> &partID) {

	printf("loading partAddrFile\n");
	FILE *paFile;
	paFile = fopen(filename, "r+");
	CheckFile(paFile, filename);
	int nOfNode, nodeID;
	fread(&nOfNode, sizeof(int), 1, paFile);
	for (int i = 0; i<nOfNode; i++) {
		PartAddr pa;
		fread(&nodeID, sizeof(int), 1, paFile);
		fread(&pa.part, sizeof(int), 1, paFile);
		fread(&pa.addr, sizeof(int), 1, paFile);
		partID[nodeID] = pa;
	}
	fclose(paFile);
}

void initialQuery(const char* fileprefix) {
	nOfDominateTest = 0;
	nOfEdgeExpended = 0;
	nOfPOIVisited = 0;
	cS.clear();
	rS.clear();
	edgeStates.clear();
	distTQ.clear();
	
	char tmpFileName[255];

	if (algorithmId == EA) { // EA

	} //---------------M-- load egtree and partAddr
	else {
		pathID.clear();
		// load egtree
		sprintf(tmpFileName, "%s\\egindex.eg_inx", fileprefix);
		egtree_load(tmpFileName, EGT);
		// load partAddr
		sprintf(tmpFileName, "%s\\part.inf", fileprefix);
		partAddrLoad(tmpFileName, paID);
		// compute the distance from bottom to up
		bottomTUpDist(Q, EGT);
		if (algorithmId == EGBU) {

		}
		else {
			visitedTNSet.clear(); 
			cSTD.clear();		
		}

	}
}

void printResult() {
	for (int i = 0; i<rS.size(); i++) {
		cout << rS[i] << " ";
	}
}

void EABasedAlgorithm(const QueryPoint &Q) {
	//Data
	int Ni = Q.Ni, Nj = Q.Nj;//Ni must less than Nj
	int AdjGrpAddr, AdjListSize, NewNodeID, PtGrpKey, PtNumOnEdge;
	//unsigned long long keywords;
	set<int> tempKwds;
	int *tempKwd;
	int nKwd;
	float EdgeDist, PtDist;
	//Precess sub edge of Q.Ni,Q.Nj divided by Q
	AdjGrpAddr = getAdjListGrpAddr(Ni);
	getFixedF(SIZE_A, Ref(AdjListSize), AdjGrpAddr);
	for (int i = 0; i < AdjListSize; i++) {
		getVarE(ADJNODE_A, Ref(NewNodeID), AdjGrpAddr, i);
		//getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
		if (NewNodeID == Nj) {
			getVarE(DIST_A, Ref(EdgeDist), AdjGrpAddr, i);
			getVarE(PTKEY_A, Ref(PtGrpKey), AdjGrpAddr, i);

			//nOfEdgeExpended++;
			//cout<<"("<<Ni<<","<<Nj<<") EdgeDist:"<<EdgeDist;

			if (PtGrpKey == -1) {
				//cout<<"No POI existed on Edge where Q located."<<endl;
			}
			else {
				getFixedF(SIZE_P, Ref(PtNumOnEdge), PtGrpKey);
				nOfEdgeExpended++;
				//cout<<" PtNum on Edge where Q located:"<<PtNumOnEdge<<endl;
				//Notice the order POI visited on edge
				for (int j = 0; j<PtNumOnEdge; j++) {
					nOfPOIVisited++;
					getVarE(PT_DIST, Ref(PtDist), PtGrpKey, j);
					getVarE(PT_NKWD, Ref(nKwd), PtGrpKey, j);
					getVarE(PT_KWD, tempKwd, PtGrpKey, j);
					//int num = sizeof(tempKwd) / sizeof(int);
					for (int loop = 0; loop<nKwd; loop++) {
						int temp;
						memcpy(&temp, tempKwd + loop*sizeof(int), sizeof(int));
						tempKwds.insert(temp);
					}

					if (LcontianRIKwd(tempKwds, Q.kwd)) {
						POI tmp;
						tmp.dist_toquery = fabs(PtDist - Q.dist_Ni);
						//getVarE(PT_ATTRIBUTE, tmp.attr, PtGrpKey, j);
						//getVarE(PT_P, &tmp.poid, PtGrpKey, j);
						if (tmp.dist_toquery <= Q.distCnst) {
							getVarE(PT_ATTRIBUTE, tmp.attr, PtGrpKey, j);
							getVarE(PT_P, &tmp.poid, PtGrpKey, j);
							updatecSDominate(tmp, cS);
						}
					}// endif
				}

			}// endelse

			break; // exit the loop
		}
	}
	// record the visited state of edge
	//edgePair epi, epj;
	int edgeKey; // key 
	edgeState esi; //value
	//epi.Ni = 0;
	//epi.Nj = Q.Ni;
	edgeKey = getKey(Ni, Nj);
	esi.vState = visited;
	esi.iDisToQuery = Q.dist_Ni;
	esi.jDisToQuery = EdgeDist - Q.dist_Ni;
	edgeStates[edgeKey] = esi;
	//epj.Ni = 0;
	//epj.Nj = Q.Nj;
	//esj.iDisToQuery = 0;
	//esj.jDisToQuery = EdgeDist - Q.dist_Ni;
	//edgeStates[epi] = esi;

	// record the visited node information
	distTQ[Ni] = Q.dist_Ni;
	distTQ[Nj] = EdgeDist - Q.dist_Ni;
	/*
	dijkVisit dvi, dvj;
	dvi.N = Ni;
	dvi.disTQ = Q.dist_Ni;

	dvj.N = Nj;
	dvj.disTQ = EdgeDist - Q.dist_Ni;
	if(dvi.disTQ < Q.distCnst) dvq.push(dvi);
	if (dvj.disTQ < Q.distCnst) dvq.push(dvj);
	*/
	// traverse the vertex in the dijkstra order
	float minDist;
	int minVertexID;

	while (!distTQ.empty()) {
		//dijkVisit tmpDV = dvq.top();
		//dvq.pop();
		minDist = -1;
		map<int, float>::iterator itd= distTQ.begin();
		for (; itd != distTQ.end(); itd++) {
			if (minDist == -1) {
				minVertexID = itd->first;
				minDist = itd->second;
			}
			else {
				if (itd->second < minDist) {
					minDist = itd->second;
					minVertexID = itd->first;
				}
			}
		}
		if (minDist >= Q.distCnst) break;
		//visited.insert(minpos);
		visitedVtx.insert(minVertexID);
		distTQ.erase(minVertexID);

		// handle the adjacent edge of this vertex
		AdjGrpAddr = getAdjListGrpAddr(minVertexID);
		getFixedF(SIZE_A, Ref(AdjListSize), AdjGrpAddr);
		for (int i = 0; i<AdjListSize; i++) {
			getVarE(ADJNODE_A, Ref(NewNodeID), AdjGrpAddr, i);
			getVarE(DIST_A, Ref(EdgeDist), AdjGrpAddr, i);
			getVarE(PTKEY_A, Ref(PtGrpKey), AdjGrpAddr, i);
			
			//edgePair ep;
			int ek;
			ek = getKey(minVertexID, NewNodeID);
			if (edgeStates.find(ek) != edgeStates.end()) { // ek has been visited
				if (edgeStates[ek].vState == hVisited) {// 计算两端的POI，修改状态
					edgeStates[ek].vState == visited;
					//---添加距离edgeStates
					if (minVertexID < NewNodeID) {
						edgeStates[ek].iDisToQuery = minDist;
					}
					else {
						edgeStates[ek].jDisToQuery = minDist;
					}
					//从两端加入POIs
					if (PtGrpKey == -1) {
						//cout<<"No POI existed on Edge where Q located."<<endl;
					}
					else {
						nOfEdgeExpended++;
						getFixedF(SIZE_P, Ref(PtNumOnEdge), PtGrpKey);
						//cout<<" PtNum on Edge where Q located:"<<PtNumOnEdge<<endl;
						//Notice the order POI visited on edge
						for (int j = 0; j<PtNumOnEdge; j++) {
							nOfPOIVisited++;
							getVarE(PT_DIST, Ref(PtDist), PtGrpKey, j);
							getVarE(PT_NKWD, Ref(nKwd), PtGrpKey, j);
							getVarE(PT_KWD, tempKwd, PtGrpKey, j);
							//int num = sizeof(tempKwd) / sizeof(int);
							for (int loop = 0; loop<nKwd; loop++) {
								int temp;
								memcpy(&temp, tempKwd + loop*sizeof(int), sizeof(int));
								tempKwds.insert(temp);
							}

							if (LcontianRIKwd(tempKwds, Q.kwd)) {
								POI tmp;
								//tmp.dist_toquery = fabs(PtDist-Q.dist_Ni);
								float d1, d2;
								d1 = edgeStates[ek].iDisToQuery + PtDist;
								d2 = edgeStates[ek].jDisToQuery + EdgeDist - PtDist;
								tmp.dist_toquery = d1<d2 ? d1 : d2;

								if (tmp.dist_toquery <= Q.distCnst) {
									getVarE(PT_ATTRIBUTE, tmp.attr, PtGrpKey, j);
									getVarE(PT_P, &tmp.poid, PtGrpKey, j);
									updatecSDominate(tmp, cS);
								}
							}// endif
						}
					}
				}
				else { //修改距离
					   //---添加距离edgeStates
					if (minVertexID < NewNodeID) {
						edgeStates[ek].iDisToQuery = minDist;
					}
					else {
						edgeStates[ek].jDisToQuery = minDist;
					}
				}
			}
			else { // ek has not been visited
				//edgePair ep;
				edgeState es;
				if (minVertexID < NewNodeID) {
					//ep.Ni = tmpDV.N;
					//ep.Nj = NewNodeID;
					es.iDisToQuery = minDist;
					es.jDisToQuery = INFINITE_MAX;
				}
				else {
					//ep.Ni = NewNodeID;
					//ep.Nj = tmpDV.N;
					es.iDisToQuery = INFINITE_MAX;
					es.jDisToQuery = minDist;
				}

				// 如何处理dijks问题？？？

				if (visitedVtx.find(NewNodeID) != visitedVtx.end()) {
					if (distTQ[NewNodeID]>(minDist + EdgeDist)) {
						distTQ[NewNodeID] = minDist + EdgeDist;
					}
				}
				else {
					distTQ[NewNodeID] = minDist + EdgeDist;
				}
				// label this edge state
				if (distTQ[NewNodeID] <= Q.distCnst) {
					//将整条边加入
					es.vState = visited;
					edgeStates[ek] = es;
					
					if (PtGrpKey == -1) {
						//cout<<"No POI existed on Edge where Q located."<<endl;
					}
					else {
						nOfEdgeExpended++;
						getFixedF(SIZE_P, Ref(PtNumOnEdge), PtGrpKey);
						//cout<<" PtNum on Edge where Q located:"<<PtNumOnEdge<<endl;
						//Notice the order POI visited on edge
						for (int j = 0; j<PtNumOnEdge; j++) {
							nOfPOIVisited++;
							getVarE(PT_DIST, Ref(PtDist), PtGrpKey, j);
							getVarE(PT_NKWD, Ref(nKwd), PtGrpKey, j);
							getVarE(PT_KWD, tempKwd, PtGrpKey, j);
							//int num = sizeof(tempKwd) / sizeof(int);
							for (int loop = 0; loop<nKwd; loop++) {
								int temp;
								memcpy(&temp, tempKwd + loop*sizeof(int), sizeof(int));
								tempKwds.insert(temp);
							}

							if (LcontianRIKwd(tempKwds, Q.kwd)) {
								POI tmp;
								//tmp.dist_toquery = fabs(PtDist-Q.dist_Ni);
								float d1, d2;
								d1 = edgeStates[ek].iDisToQuery + PtDist;
								d2 = edgeStates[ek].jDisToQuery + EdgeDist - PtDist;
								tmp.dist_toquery = d1<d2 ? d1 : d2;

								if (tmp.dist_toquery <= Q.distCnst) {
									getVarE(PT_ATTRIBUTE, tmp.attr, PtGrpKey, j);
									getVarE(PT_P, &tmp.poid, PtGrpKey, j);
									updatecSDominate(tmp, cS);
								}
							}// endif
						}
					}
				}
				else {// 如果只有部分在边上，hVisited
					es.vState = hVisited;
					edgeStates[ek] = es;
				}
			}
		}// end for
		
	}// end while

	 // 处理hVisited

	map<int, edgeState>::iterator iterhV = edgeStates.begin();
	for (; iterhV != edgeStates.end(); ++iterhV) {
		if (iterhV->second.vState == hVisited) {//处理只有距离当前节点满足距离约束的POI
			int nodei, nodej;
			int brknodei, brknodej;
			breakKey(iterhV->first, brknodei, brknodej);
			if (iterhV->second.iDisToQuery == INFINITE_MAX) {
				nodei = brknodej;
				nodej = brknodei;
			}
			else {
				nodei = brknodei;
				nodej = brknodej;
			}
			//获取边上的POI
			AdjGrpAddr = getAdjListGrpAddr(nodei);
			getFixedF(SIZE_A, Ref(AdjListSize), AdjGrpAddr);
			for (int i = 0; i<AdjListSize; i++) {
				getVarE(ADJNODE_A, Ref(NewNodeID), AdjGrpAddr, i);
				//getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
				if (NewNodeID == nodej) {
					getVarE(DIST_A, Ref(EdgeDist), AdjGrpAddr, i);
					getVarE(PTKEY_A, Ref(PtGrpKey), AdjGrpAddr, i);

					//cout<<"("<<Ni<<","<<Nj<<") EdgeDist:"<<EdgeDist;
					if (PtGrpKey == -1) {
						//cout<<"No POI existed on Edge where Q located."<<endl;
					}
					else {
						nOfEdgeExpended++;
						getFixedF(SIZE_P, Ref(PtNumOnEdge), PtGrpKey);
						//cout<<" PtNum on Edge where Q located:"<<PtNumOnEdge<<endl;
						//Notice the order POI visited on edge
						for (int j = 0; j<PtNumOnEdge; j++) {
							nOfPOIVisited++;
							getVarE(PT_DIST, Ref(PtDist), PtGrpKey, j);
							getVarE(PT_NKWD, Ref(nKwd), PtGrpKey, j);
							getVarE(PT_KWD, tempKwd, PtGrpKey, j);
							//int num = sizeof(tempKwd) / sizeof(int);
							for (int loop = 0; loop<nKwd; loop++) {
								int temp;
								memcpy(&temp, tempKwd + loop*sizeof(int), sizeof(int));
								tempKwds.insert(temp);
							}

							if (LcontianRIKwd(tempKwds, Q.kwd)) {
								POI tmp;
								if (nodei<nodej) {
									tmp.dist_toquery = iterhV->second.iDisToQuery + PtDist;
								}
								else {
									tmp.dist_toquery = iterhV->second.jDisToQuery + EdgeDist - PtDist;
								}

								if (tmp.dist_toquery <= Q.distCnst) {
									getVarE(PT_ATTRIBUTE, tmp.attr, PtGrpKey, j);
									getVarE(PT_P, &tmp.poid, PtGrpKey, j);
									updatecSDominate(tmp, cS);
								}
							}// endif
						}
					}// endelse
				}
			}
		}
	}

	for (int i = 0; i<cS.size(); i++) {
		rS.push_back(cS[i].poid);
		//cout<<cS[i].poid<<" ";
	}
	//cout<<endl;
	//cout<<"POI Read           #:"<<nOfPOIVisited<<endl;  
	//cout<<"Road Dominate Tested   #:"<<nOfDominateTest<<endl;
	//cout<<"Road Edge Expanded #:"<<nOfEdgeExpended<<endl;
	//base2numedgeexpand+=edgeexpanded;
	//base2numnodeexpand+=roadnodecnt;
}

void bottomTUpDist(QueryPoint Q, vector<TreeNode> &EGT) {
	//bottomToUp 由下至上构建distTQ，保存路径

	int pid = paID[Q.Ni].part; //假设查询点位于treenode 内部
	pathID.push_back(pid);
	vector<int>::iterator it;
	//计算两端点到边界的距离，去最小
	EGT[pid].minDistTQ = INFINITE_MAX;
	for (int i = 0; i<EGT[pid].borders.size(); i++) {
		distTQ[i] = INFINITE_MAX;
	}
	int posi, posj, size = EGT[pid].leafnodes.size();
	it = find(EGT[pid].leafnodes.begin(), EGT[pid].leafnodes.end(), Q.Ni);
	posi = (*it);
	it = find(EGT[pid].leafnodes.begin(), EGT[pid].leafnodes.end(), Q.Nj);
	posj = (*it);
	for (int i = 0; i<EGT[pid].borders.size(); i++) {
		float distI = Q.dist_Ni + EGT[pid].mind[i*size + posi];
		float distJ = queryEdgeDist - Q.dist_Ni + EGT[pid].mind[i*size + posj];
		if (distI<distJ) {
			EGT[pid].distTQ[i] = distI;
			if (distI<EGT[pid].minDistTQ) EGT[pid].minDistTQ = distI;
		}
		else {
			EGT[pid].distTQ[i] = distJ;
			if (distJ<EGT[pid].minDistTQ) EGT[pid].minDistTQ = distJ;
		}
	}
	int preID;
	while (pid != 0 && EGT[pid].minDistTQ<Q.distCnst) {
		preID = pid;
		pid = EGT[pid].father;
		pathID.push_back(pid);
		//计算父节点的距离
		//vector<int>::iterator it;
		//计算两端点到边界的距离，去最小
		EGT[pid].minDistTQ = INFINITE_MAX;
		for (int i = 0; i<EGT[pid].borders.size(); i++) {
			EGT[pid].distTQ[i] = INFINITE_MAX;
		}
		size = EGT[pid].union_borders.size();
		int posb, posu;
		//对于该层的每个border，计算从下一层到该层的最短距离
		for (int i = 0; i<EGT[pid].borders.size(); i++) {
			for (int j = 0; j<EGT[preID].borders.size(); j++) {
				//定位两个node在union_border中位置		
				it = find(EGT[pid].union_borders.begin(), EGT[pid].union_borders.end(), EGT[preID].borders[j]);
				posb = (*it);
				it = find(EGT[pid].union_borders.begin(), EGT[pid].union_borders.end(), EGT[pid].borders[i]);
				posu = (*it);
				float mindPos;
				if (posb<posu) {
					mindPos = (posu*(posu + 1)) / 2 + posb;
				}
				else {
					mindPos = (posb*(posb + 1)) / 2 + posu;
				}
				if ((EGT[preID].distTQ[j] + EGT[pid].mind[mindPos])<EGT[pid].distTQ[i]) {
					EGT[pid].distTQ[i] = EGT[preID].distTQ[j] + EGT[pid].mind[mindPos];
				}
			}//endfor
			if (EGT[pid].distTQ[i]<EGT[pid].minDistTQ) EGT[pid].minDistTQ = EGT[pid].distTQ[i];
		}//endfor

	}//endwhile

}
//计算s到cans的距离
vector<float> dijkstra_candidate(const QueryPoint &Q, int tid, vector<TreeNode> &EGT) {

	//找最小共同父节点,构建最短路径
	int commonID = tid;
	vector<int> ctPathID;
	ctPathID.push_back(commonID);
	while (true) {
		if (find(pathID.begin(), pathID.end(), commonID) == pathID.end()) {
			commonID = EGT[commonID].father;
			ctPathID.push_back(commonID);
		}
		else {
			break;
		}
	}
	//寻找pathID上commonID的下一层
	int pathid = paID[Q.Ni].part;
	while (EGT[pathid].father != commonID) {
		pathid = EGT[pathid].father;
	}
	//开始右转去连接其他边，这时候需要从上到下处理了
	//单独处理右边最上面的节点
	int upMostID = ctPathID.size() - 1;
	if (EGT[upMostID].distTQ.size() != EGT[upMostID].borders.size()) {
		EGT[upMostID].minDistTQ = INFINITE_MAX;
		for (int i = 0; i<EGT[upMostID].borders.size(); i++) {
			EGT[upMostID].distTQ[i] = INFINITE_MAX;
		}
		int posutm, pospth, sizeC = EGT[commonID].union_borders.size();
		for (int i = 0; i<EGT[upMostID].borders.size(); i++) {
			for (int j = 0; j<EGT[pathid].borders.size(); j++) {
				//定位
				vector<int>::iterator itr;
				itr = find(EGT[commonID].union_borders.begin(), EGT[commonID].union_borders.end(), EGT[upMostID].borders[i]);
				posutm = (*itr);
				itr = find(EGT[commonID].union_borders.begin(), EGT[commonID].union_borders.end(), EGT[pathid].borders[j]);
				pospth = (*itr);
				float mindPos;
				if (posutm<pospth) {
					mindPos = (pospth*(pospth + 1)) / 2 + posutm;
				}
				else {
					mindPos = (posutm*(posutm + 1)) / 2 + pospth;
				}
				if ((EGT[pathid].distTQ[j] + EGT[upMostID].mind[mindPos])<EGT[upMostID].distTQ[i]) {
					EGT[upMostID].distTQ[i] = EGT[pathid].distTQ[j] + EGT[upMostID].mind[mindPos];
				}
			}
			if (EGT[upMostID].distTQ[i]<EGT[upMostID].minDistTQ) EGT[upMostID].minDistTQ = EGT[upMostID].distTQ[i];
		}
	}
	//从导数第二个开始
	int uBID, fID;
	for (int k = ctPathID.size() - 2; k >= 0; k--) {
		uBID = ctPathID[k];
		fID = EGT[uBID].father;
		//preID = pid;
		//pid=EGT[pid].father;
		//pathID.push_back(pid);
		//计算父节点的距离
		//vector<int>::iterator it;
		//计算两端点到边界的距离，去最小
		if (EGT[uBID].distTQ.size() == EGT[uBID].borders.size()) continue;
		EGT[uBID].minDistTQ = INFINITE_MAX;
		for (int i = 0; i<EGT[uBID].borders.size(); i++) {
			EGT[uBID].distTQ[i] = INFINITE_MAX;
		}
		int sizef = EGT[fID].union_borders.size();
		int posdc, posdf;
		//对于该层的每个border，计算从下一层到该层的最短距离
		for (int i = 0; i<EGT[uBID].borders.size(); i++) {
			for (int j = 0; j<EGT[fID].borders.size(); j++) {
				//定位两个node在union_border中位置		
				vector<int>::iterator it;
				it = find(EGT[fID].union_borders.begin(), EGT[fID].union_borders.end(), EGT[uBID].borders[i]);
				posdc = (*it);
				it = find(EGT[fID].union_borders.begin(), EGT[fID].union_borders.end(), EGT[fID].borders[j]);
				posdf = (*it);
				float mindPos;
				if (posdc<posdf) {
					mindPos = (posdf*(posdf + 1)) / 2 + posdc;
				}
				else {
					mindPos = (posdc*(posdc + 1)) / 2 + posdf;
				}
				if ((EGT[fID].distTQ[j] + EGT[fID].mind[mindPos])<EGT[uBID].distTQ[i]) {
					EGT[uBID].distTQ[i] = EGT[fID].distTQ[j] + EGT[fID].mind[mindPos];
				}
			}
			if (EGT[uBID].distTQ[i]<EGT[uBID].minDistTQ) EGT[uBID].minDistTQ = EGT[uBID].distTQ[i];
		}//endfor

	}//endfor
	return EGT[tid].distTQ;
	//求最短路径
}

//------------------for EGTD------------------------
bool LEditKwdContainR(set<set<int>> editKwd, set<int> qkwd) {
	set<set<int>> ::iterator it = editKwd.begin();
	for (; it != editKwd.end(); it++) {
		set<int> editK = *it;
		if (editK.size() < qkwd.size()) return false;
		if (editDistanceRTL(editK, qkwd) < editK.size()*edDis) { // maybe exist
			return true;
		}
	}
	return false;
}

bool canLSDomiNode(CandTDValue ct, locSkyline ls) {
	bool cDl = true;
	bool lDc = true;
	int equal = 0;

	for (int i = 0; i < ATTRIBUTE_DIMENSION; i++) {
		if (Q.subSpace[i] == 1) {
			//if (ct.ls.blockID[i] < ls.blockID[i]) lDc = false;
			if (ct.ls.blockID[i] < ls.blockID[i]) cDl = false;
			if (ct.ls.blockID[i] == ls.blockID[i]) equal++;
		}
	}
	if (!cDl) return false;
	if (cDl&&equal != Q.subSpace.count()) return true;
	if (equal == Q.subSpace.count()) { //两个相等
		int less = 0;
		int eq = 0;
		for (int i = 0; i < ATTRIBUTE_DIMENSION; i++) {
			if (Q.subSpace[i] == 1) {
				if (ct.ls.skylineUpBound[i] < ls.skylineBtBound[i]) less++;
				if (ct.ls.skylineUpBound[i] == ls.skylineBtBound[i]) eq++;
			}
		}
		if ((less + eq) == Q.subSpace.count() && eq>0) return true;
		return false;
	}
}
// test whether candidate set dominates the locskyline of node
bool canLocSkylineNode(map<int, locSkyline> locSky) {
	// 是否在cSTD中存在一个ID支配locSky
	map<int, CandTDValue> ::iterator itCTD = cSTD.begin();

	for (; itCTD != cSTD.begin(); itCTD++) {
		bool dominateAll = true;
		map<int, locSkyline> ::iterator itSky = locSky.begin();
		for (; itSky != locSky.end(); itSky++) {
			if (!canLSDomiNode(itCTD->second, itSky->second)) {//不支配，停止
				dominateAll = false;
				break;
			}
		}
		if (dominateAll) return true;
	}
	return false;
}

void EGTDAlgorithm(const QueryPoint &Q) {
	//process leadnode, two situation:1, not in same part 2, in the same part
	int Ni = Q.Ni, Nj = Q.Nj;//Ni must less than Nj
	int AdjGrpAddr, AdjListSize, NewNodeID, PtGrpKey, PtNumOnEdge;


	set<int> tempKwds;
	set<int> visNodes;
	int *tempKwd;
	set<int> edgeSumKwds;
	float edgeSumAttr[ATTRIBUTE_DIMENSION][2];
	float EdgeDist, PtDist;
	// locate the treeNodeID
	int trNodeIDi, trNodeIDj;
	int currTNID;
	float minDist;
	dijkVisit dvi, dvj;

	// get the upmost nodeID
	int sizeT = pathID.size();
	int upid = pathID[sizeT - 1];
	int cVistPID = paID[Ni].part;
	int currID = upid;
	tnq.push(EGTree[upid]);
	// while loop
	while (!tnq.empty()) {
		TreeNode tn = tnq.top();
		tnq.pop();
		//***************************
		if (tn.isleaf) { // if is leaf node, there are two situations
			if (currID == cVistPID) { // if in the query leaf node
				AdjGrpAddr = getAdjListGrpAddr(Ni);
				getFixedF(SIZE_A, Ref(AdjListSize), AdjGrpAddr);
				AdjGrpAddr = getAdjListGrpAddr(Ni);
				getFixedF(SIZE_A, Ref(AdjListSize), AdjGrpAddr);
				for (int i = 0; i<AdjListSize; i++) {
					getVarE(ADJNODE_A, Ref(NewNodeID), AdjGrpAddr, i);
					//getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
					if (NewNodeID == Nj) {
						getVarE(DIST_A, Ref(EdgeDist), AdjGrpAddr, i);
						//记录查询边距离
						queryEdgeDist = EdgeDist;
						getVarE(PTKEY_A, Ref(PtGrpKey), AdjGrpAddr, i);
						getVarE(SUMKWD_A, &edgeSumKwds.begin(), AdjGrpAddr, i);
						getVarE(SUMATTR_A, edgeSumAttr, AdjGrpAddr, i);

						edgePair epi, epj;
						edgeState esi, esj;
						epi.Ni = 0;
						epi.Nj = Q.Ni;
						esi.vState = visited;
						esi.iDisToQuery = 0;
						esi.jDisToQuery = Q.dist_Ni;

						epj.Ni = 0;
						epj.Nj = Q.Nj;
						esj.iDisToQuery = 0;
						esj.jDisToQuery = EdgeDist - Q.dist_Ni;

						edgeStates[epi] = esi;
						edgeStates[epj] = esj;

						if (!LcontianRIKwd(edgeSumKwds, Q.kwd)) {

							break;
						}
						if (isBeBoundDominate(edgeSumAttr, cS, Q)) {

							break;
						}
						nOfEdgeExpended++;

						//cout<<"("<<Ni<<","<<Nj<<") EdgeDist:"<<EdgeDist;
						if (PtGrpKey == -1) {
							//cout<<"No POI existed on Edge where Q located."<<endl;
						}
						else {
							getFixedF(SIZE_P, Ref(PtNumOnEdge), PtGrpKey);
							//cout<<" PtNum on Edge where Q located:"<<PtNumOnEdge<<endl;
							//Notice the order POI visited on edge
							for (int j = 0; j<PtNumOnEdge; j++) {
								nOfPOIVisited++;
								getVarE(PT_DIST, Ref(PtDist), PtGrpKey, j);
								getVarE(PT_KWD, tempKwd, PtGrpKey, j);

								int num = sizeof(tempKwd) / sizeof(int);
								for (int loop = 0; loop<sizeof(tempKwd); loop + sizeof(int)) {
									int temp;
									memcpy(&temp, tempKwd + loop, sizeof(int));
									tempKwds.insert(temp);
								}

								if (LcontianRIKwd(tempKwds, Q.kwd)) {
									POI tmp;
									tmp.dist_toquery = fabs(PtDist - Q.dist_Ni);
									getVarE(PT_ATTRIBUTE, tmp.attr, PtGrpKey, j);
									getVarE(PT_P, &tmp.poid, PtGrpKey, j);

									if (tmp.dist_toquery <= Q.distCnst) {
										updatecSDominate(tmp, cS);
									}
								}// endif
							}

						}// endelse
					}
				}

				trNodeIDi = paID[Ni].part;
				trNodeIDj = paID[Nj].part;
				if (trNodeIDi == trNodeIDj) { //两个在一个划分块里面
					currTNID = trNodeIDi;
					dvi.N = Q.Ni;
					dvi.disTQ = Q.dist_Ni;

					dvj.N = Q.Nj;
					dvj.disTQ = EdgeDist - Q.dist_Ni;

					dvq.push(dvi);
					dvq.push(dvj);
				} /*else { //查询点在边界边上面
				  if(trNodeIDi < trNodeIDj) {
				  currTNID = trNodeIDi;
				  dvi.N = Q.Ni;
				  dvi.disTQ = Q.dist_Ni;
				  dvq.push(dvi);
				  } else {
				  currTNID = trNodeIDj;
				  dvj.N = Q.Nj;
				  dvj.disTQ = EdgeDist-Q.dist_Ni;
				  dvq.push(dvj);
				  }
				  }*/
				  //------------&&&&&&&&&&&&&&内部节点处理
				while (!dvq.empty()) {
					dijkVisit tmpDV = dvq.top();
					dvq.pop();
					if (paID[tmpDV.N].part != currTNID) continue;
					if (tmpDV.disTQ <= Q.distCnst) {
						AdjGrpAddr = getAdjListGrpAddr(tmpDV.N);
						getFixedF(SIZE_A, Ref(AdjListSize), AdjGrpAddr);
						for (int i = 0; i<AdjListSize; i++) {
							getVarE(ADJNODE_A, Ref(NewNodeID), AdjGrpAddr, i);
							getVarE(DIST_A, Ref(EdgeDist), AdjGrpAddr, i);
							getVarE(PTKEY_A, Ref(PtGrpKey), AdjGrpAddr, i);
							edgePair ep;
							ep.Ni = NewNodeID<tmpDV.N ? NewNodeID : tmpDV.N;
							ep.Nj = NewNodeID>tmpDV.N ? NewNodeID : tmpDV.N;
							map<edgePair, edgeState, eSComparison>::iterator iter;
							iter = edgeStates.find(ep);
							if (iter != edgeStates.end()) {
								if (edgeStates[ep].vState == hVisited) {// 计算两端的POI，修改状态
									edgeStates[ep].vState == visited;
									//---添加距离edgeStates
									if (tmpDV.N < NewNodeID) {
										edgeStates[ep].iDisToQuery = tmpDV.disTQ;
									}
									else {
										edgeStates[ep].jDisToQuery = tmpDV.disTQ;
									}
									//从两端加入POIs
									if (PtGrpKey == -1) {
										//cout<<"No POI existed on Edge where Q located."<<endl;
									}
									else {

										getFixedF(SIZE_P, Ref(PtNumOnEdge), PtGrpKey);
										//cout<<" PtNum on Edge where Q located:"<<PtNumOnEdge<<endl;
										//Notice the order POI visited on edge
										for (int j = 0; j<PtNumOnEdge; j++) {
											nOfPOIVisited++;
											getVarE(PT_DIST, Ref(PtDist), PtGrpKey, j);

											getVarE(PT_KWD, tempKwd, PtGrpKey, j);
											int num = sizeof(tempKwd) / sizeof(int);
											for (int loop = 0; loop<sizeof(tempKwd); loop + sizeof(int)) {
												int temp;
												memcpy(&temp, tempKwd + loop, sizeof(int));
												tempKwds.insert(temp);
											}

											if (LcontianRIKwd(tempKwds, Q.kwd)) {
												POI tmp;
												//tmp.dist_toquery = fabs(PtDist-Q.dist_Ni);
												float d1, d2;
												d1 = edgeStates[ep].iDisToQuery + PtDist;
												d2 = edgeStates[ep].jDisToQuery + EdgeDist - PtDist;
												tmp.dist_toquery = d1<d2 ? d1 : d2;

												getVarE(PT_ATTRIBUTE, tmp.attr, PtGrpKey, j);
												getVarE(PT_P, &tmp.poid, PtGrpKey, j);

												if (tmp.dist_toquery <= Q.distCnst) {
													updatecSDominate(tmp, cS);
												}
											}// endif
										}
									}
								}
								else if (edgeStates[ep].vState == visited) { //修改距离
																			   //---添加距离edgeStates
									if (tmpDV.N < NewNodeID) {
										edgeStates[ep].iDisToQuery = tmpDV.disTQ;
									}
									else {
										edgeStates[ep].jDisToQuery = tmpDV.disTQ;
									}
								}
								else {

								}
							}
							else {
								edgePair ep;
								edgeState es;
								if (tmpDV.N < NewNodeID) {
									ep.Ni = tmpDV.N;
									ep.Nj = NewNodeID;
									es.iDisToQuery = tmpDV.disTQ;
									es.jDisToQuery = INFINITE_MAX;
								}
								else {
									ep.Ni = NewNodeID;
									ep.Nj = tmpDV.N;
									es.iDisToQuery = INFINITE_MAX;
									es.jDisToQuery = tmpDV.disTQ;
								}
								if (distTQ.find(NewNodeID) != distTQ.end()) {
									if (distTQ[NewNodeID]>(tmpDV.disTQ + EdgeDist)) {
										distTQ[NewNodeID] = tmpDV.disTQ + EdgeDist;
									}
								}
								else {
									distTQ[NewNodeID] = tmpDV.disTQ + EdgeDist;
								}
								if (distTQ[NewNodeID] <= Q.distCnst) {
									//将整条边加入
									es.vState = visited;
									edgeStates[ep] = es;



									if (PtGrpKey == -1) {
										//cout<<"No POI existed on Edge where Q located."<<endl;
									}
									else {
										getFixedF(SIZE_P, Ref(PtNumOnEdge), PtGrpKey);
										//获得关键字信息和属性信息，直接过滤
										if (!LcontianRIKwd(edgeSumKwds, Q.kwd)) {
											//这条边不加
											continue;
										}
										if (isBeBoundDominate(edgeSumAttr, cS, Q)) {
											//这条边不加
											continue;
										}
										//cout<<" PtNum on Edge where Q located:"<<PtNumOnEdge<<endl;
										//Notice the order POI visited on edge
										for (int j = 0; j<PtNumOnEdge; j++) {
											nOfPOIVisited++;
											getVarE(PT_DIST, Ref(PtDist), PtGrpKey, j);

											getVarE(PT_KWD, tempKwd, PtGrpKey, j);
											int num = sizeof(tempKwd) / sizeof(int);
											for (int loop = 0; loop<sizeof(tempKwd); loop + sizeof(int)) {
												int temp;
												memcpy(&temp, tempKwd + loop, sizeof(int));
												tempKwds.insert(temp);
											}

											if (LcontianRIKwd(tempKwds, Q.kwd)) {
												POI tmp; //这时候可能不是真实距离,也不会影响结果
												if (tmpDV.N<NewNodeID) {
													tmp.dist_toquery = tmpDV.disTQ + PtDist;
												}
												else {
													tmp.dist_toquery = tmpDV.disTQ + EdgeDist - PtDist;
												}


												getVarE(PT_ATTRIBUTE, tmp.attr, PtGrpKey, j);
												getVarE(PT_P, &tmp.poid, PtGrpKey, j);

												if (tmp.dist_toquery <= Q.distCnst) {
													updatecSDominate(tmp, cS);
												}
											}// endif
										}
									}

								}
								else {// 如果只有部分在边上，hVisited
									if (!LcontianRIKwd(edgeSumKwds, Q.kwd)) {
										//这条边不加
										es.vState = visited;
										edgeStates[ep] = es;
										continue;
									}
									if (isBeBoundDominate(edgeSumAttr, cS, Q)) {
										//这条边不加
										es.vState = visited;
										edgeStates[ep] = es;
										continue;
									}
									es.vState = hVisited;
									edgeStates[ep] = es;
								}
							}

							//getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
							//getVarE(PTKEY_A,Ref(PtGrpKey),AdjGrpAddr,i);
							if (PtGrpKey == -1) {
								//cout<<"No POI existed on Edge where Q located."<<endl;
							}
							else {

								getFixedF(SIZE_P, Ref(PtNumOnEdge), PtGrpKey);
								//cout<<" PtNum on Edge where Q located:"<<PtNumOnEdge<<endl;
								//Notice the order POI visited on edge
								for (int j = 0; j<PtNumOnEdge; j++) {
									nOfPOIVisited++;
									getVarE(PT_DIST, Ref(PtDist), PtGrpKey, j);

									getVarE(PT_KWD, tempKwd, PtGrpKey, j);
									int num = sizeof(tempKwd) / sizeof(int);
									for (int loop = 0; loop<sizeof(tempKwd); loop + sizeof(int)) {
										int temp;
										memcpy(&temp, tempKwd + loop, sizeof(int));
										tempKwds.insert(temp);
									}

									if (LcontianRIKwd(tempKwds, Q.kwd)) {
										POI tmp;
										tmp.dist_toquery = fabs(PtDist - Q.dist_Ni);
										getVarE(PT_ATTRIBUTE, tmp.attr, PtGrpKey, j);
										getVarE(PT_P, &tmp.poid, PtGrpKey, j);

										if (tmp.dist_toquery <= Q.distCnst) {
											updatecSDominate(tmp, cS);
										}
									}// endif
								}
							}


						}// end for
					}// endif
				}// end while

			}
			else { // not in the query leaf node
				   // utilize the mind to retrieve this node meet distance constraint
				   //---------&&&&&&&&&&&&计算满足距离约束的node集合，并
				vector<float> tempDist;
				for (int k = 0; k<tn.leafnodes.size(); k++) {
					tempDist.push_back(INFINITE_MAX);
				}
				for (int k = 0; k<tn.refDistTQ.size(); k++) {
					//float mindist = INFINITE_MAX;
					if (tn.refDistTQ[k]<0) continue;
					int pos, size = tn.leafnodes.size();
					float tempD;
					for (int t = 0; t<tn.leafnodes.size(); t++) {
						pos = k*size + t;
						tempD = Q.distCnst - tn.refDistTQ[k] + tn.mind[pos];
						if (tempD <= Q.distCnst) {
							if (tempD<tempDist[t]) tempDist[t] = tempD;
						}
					}
				}
				//构建dvq
				//--------------------xiugai----------
				dijkVisit tmpDV = dvq.top();
				for (int k = 0; k<tempDist.size(); k++) {
					if (tempDist[k]>Q.distCnst) continue;
					//为每个符合条件的点展开
					/*dijkVisit dv;
					dv.N = tn.leafnodes[k];
					dv.disTQ = tempDis[k];
					dvq.push(dv);*/
					AdjGrpAddr = getAdjListGrpAddr(tn.leafnodes[k]);
					getFixedF(SIZE_A, Ref(AdjListSize), AdjGrpAddr);
					for (int i = 0; i<AdjListSize; i++) {
						getVarE(ADJNODE_A, Ref(NewNodeID), AdjGrpAddr, i);
						getVarE(DIST_A, Ref(EdgeDist), AdjGrpAddr, i);
						getVarE(PTKEY_A, Ref(PtGrpKey), AdjGrpAddr, i);
						getVarE(SUMKWD_A, &edgeSumKwds.begin(), AdjGrpAddr, i);
						getVarE(SUMATTR_A, edgeSumAttr, AdjGrpAddr, i);
						edgePair ep;
						ep.Ni = NewNodeID<tn.leafnodes[k] ? NewNodeID : tn.leafnodes[k];
						ep.Nj = NewNodeID>tn.leafnodes[k] ? NewNodeID : tn.leafnodes[k];
						map<edgePair, edgeState, eSComparison>::iterator iter;
						iter = edgeStates.find(ep);
						if (iter != edgeStates.end()) {
							if (edgeStates[ep].vState == hVisited) {// 计算两端的POI，修改状态
								edgeStates[ep].vState == visited;
								//---添加距离edgeStates
								if (tn.leafnodes[k] < NewNodeID) {
									edgeStates[ep].iDisToQuery = tempDist[k];
								}
								else {
									edgeStates[ep].jDisToQuery = tempDist[k];
								}
								//从两端加入POIs
								if (PtGrpKey == -1) {
									//cout<<"No POI existed on Edge where Q located."<<endl;
								}
								else {

									getFixedF(SIZE_P, Ref(PtNumOnEdge), PtGrpKey);
									//cout<<" PtNum on Edge where Q located:"<<PtNumOnEdge<<endl;
									//Notice the order POI visited on edge
									for (int j = 0; j<PtNumOnEdge; j++) {
										nOfPOIVisited++;
										getVarE(PT_DIST, Ref(PtDist), PtGrpKey, j);

										getVarE(PT_KWD, tempKwd, PtGrpKey, j);
										int num = sizeof(tempKwd) / sizeof(int);
										for (int loop = 0; loop<sizeof(tempKwd); loop + sizeof(int)) {
											int temp;
											memcpy(&temp, tempKwd + loop, sizeof(int));
											tempKwds.insert(temp);
										}

										if (LcontianRIKwd(tempKwds, Q.kwd)) {
											POI tmp;
											//tmp.dist_toquery = fabs(PtDist-Q.dist_Ni);
											float d1, d2;
											d1 = edgeStates[ep].iDisToQuery + PtDist;
											d2 = edgeStates[ep].jDisToQuery + EdgeDist - PtDist;
											tmp.dist_toquery = d1<d2 ? d1 : d2;

											getVarE(PT_ATTRIBUTE, tmp.attr, PtGrpKey, j);
											getVarE(PT_P, &tmp.poid, PtGrpKey, j);

											if (tmp.dist_toquery <= Q.distCnst) {
												updatecSDominate(tmp, cS);
											}
										}// endif
									}
								}
							}
							else if (edgeStates[ep].vState == visited) { //修改距离
																		   //---添加距离edgeStates
								if (tmpDV.N < NewNodeID) {
									edgeStates[ep].iDisToQuery = tempDist[k];
								}
								else {
									edgeStates[ep].jDisToQuery = tempDist[k];
								}
							}
							else {


							}
						}
						else {
							//判断是否要加入到ep中			
							edgeState es;
							if (tn.leafnodes[k] < NewNodeID) {
								ep.Ni = tn.leafnodes[k];
								ep.Nj = NewNodeID;
								es.iDisToQuery = tempDist[k];
								es.jDisToQuery = INFINITE_MAX;
							}
							else {
								ep.Ni = NewNodeID;
								ep.Nj = tn.leafnodes[k];
								es.iDisToQuery = INFINITE_MAX;
								es.jDisToQuery = tempDist[k];
							}
							if (distTQ.find(NewNodeID) != distTQ.end()) {
								if (distTQ[NewNodeID]>(tempDist[k] + EdgeDist)) {
									distTQ[NewNodeID] = tempDist[k] + EdgeDist;
								}
							}
							else {
								distTQ[NewNodeID] = tempDist[k] + EdgeDist;
							}
							if (distTQ[NewNodeID] <= Q.distCnst) {
								//将整条边加入
								es.vState = visited;
								edgeStates[ep] = es;

								if (PtGrpKey == -1) {
									//cout<<"No POI existed on Edge where Q located."<<endl;
								}
								else {
									//获得关键字信息和属性信息，直接过滤
									if (!LcontianRIKwd(edgeSumKwds, Q.kwd)) {
										//这条边不加
										continue;
									}
									if (isBeBoundDominate(edgeSumAttr, cS, Q)) {
										//这条边不加
										continue;
									}

									//cout<<" PtNum on Edge where Q located:"<<PtNumOnEdge<<endl;
									//Notice the order POI visited on edge
									for (int j = 0; j<PtNumOnEdge; j++) {
										nOfPOIVisited++;
										getVarE(PT_DIST, Ref(PtDist), PtGrpKey, j);

										getVarE(PT_KWD, tempKwd, PtGrpKey, j);
										int num = sizeof(tempKwd) / sizeof(int);
										for (int loop = 0; loop<sizeof(tempKwd); loop + sizeof(int)) {
											int temp;
											memcpy(&temp, tempKwd + loop, sizeof(int));
											tempKwds.insert(temp);
										}

										if (LcontianRIKwd(tempKwds, Q.kwd)) {
											POI tmp; //这时候可能不是真实距离,也不会影响结果
											if (tmpDV.N<NewNodeID) {
												tmp.dist_toquery = tmpDV.disTQ + PtDist;
											}
											else {
												tmp.dist_toquery = tmpDV.disTQ + EdgeDist - PtDist;
											}


											getVarE(PT_ATTRIBUTE, tmp.attr, PtGrpKey, j);
											getVarE(PT_P, &tmp.poid, PtGrpKey, j);

											if (tmp.dist_toquery <= Q.distCnst) {
												updatecSDominate(tmp, cS);
											}
										}// endif
									}
								}

							}
							else {// 如果只有部分在边上，hVisited
								  //if被关键字等过滤掉，则直接算是visited
								  //获得关键字信息和属性信息，直接过滤
								if (!LcontianRIKwd(edgeSumKwds, Q.kwd)) {
									//这条边不加
									es.vState = visited;
									edgeStates[ep] = es;
									continue;
								}
								if (isBeBoundDominate(edgeSumAttr, cS, Q)) {
									//这条边不加
									es.vState = visited;
									edgeStates[ep] = es;
									continue;
								}
								es.vState = hVisited;
								edgeStates[ep] = es;
							}
						}


					}//endif
				}//endfor
			}
		}
		else { // if is inter node
			for (int i = 0; i<tn.children.size();) {
				int cid = tn.children[i];
				currID = cid;
				if (visitedTNSet.find(cid) == visitedTNSet.end()) {
					//if is not filter by summary information
					if (!LcontianRIKwd(EGT[cid].union_kwd, Q.kwd)) {
						if (visitedTNSet.find(cid) == visitedTNSet.end()) visitedTNSet.insert(cid);
						continue;
					}
					if (isBeBoundDominate(EGT[cid].attrBound, cS, Q)) {
						if (visitedTNSet.find(cid) == visitedTNSet.end()) visitedTNSet.insert(cid);
						continue;
					}
					if (EGT[cid].minDistTQ > Q.distCnst) {
						if (visitedTNSet.find(cid) == visitedTNSet.end()) visitedTNSet.insert(cid);
						continue;
					}
					// if is not filter by editDistance and locSkyline ***********
					if (!LEditKwdContainR(EGT[cid].editKwd, Q.kwd)) {
						if (visitedTNSet.find(cid) == visitedTNSet.end()) visitedTNSet.insert(cid);
						continue;
					}
					//**********************locSky 
					if (!canLocSkylineNode(EGT[cid].locSky)) {
						if (visitedTNSet.find(cid) == visitedTNSet.end()) visitedTNSet.insert(cid);
						continue;
					}

					//compute the minDistTQ
					vector<float> distTQ;
					distTQ = dijkstra_candidate(Q, cid, EGT);
					float minDist = INFINITE_MAX;
					for (int k = 0; k<distTQ.size(); k++) {
						if (distTQ[k] < minDist) minDist = distTQ[k];
						EGT[cid].refDistTQ.push_back(Q.distCnst - distTQ[k]);
					}
					EGT[cid].minDistTQ = minDist;
					//----------------&&&&&&&
					tnq.push(EGT[cid]);
				}
			}
		}
	}

	// handel hvisited
	// 处理hVisited

	map<edgePair, edgeState, eSComparison>::iterator iterhV;
	for (iterhV = edgeStates.begin(); iterhV != edgeStates.end(); ++iterhV) {
		if (iterhV->second.vState == hVisited) {//处理只有距离当前节点满足距离约束的POI
			int nodei, nodej;
			if (iterhV->second.iDisToQuery == INFINITE_MAX) {
				nodei = iterhV->first.Nj;
				nodej = iterhV->first.Ni;
			}
			else {
				nodei = iterhV->first.Nj;
				nodej = iterhV->first.Ni;
			}
			//获取边上的POI
			AdjGrpAddr = getAdjListGrpAddr(nodei);
			getFixedF(SIZE_A, Ref(AdjListSize), AdjGrpAddr);
			for (int i = 0; i<AdjListSize; i++) {
				getVarE(ADJNODE_A, Ref(NewNodeID), AdjGrpAddr, i);
				//getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
				if (NewNodeID == nodej) {
					getVarE(DIST_A, Ref(EdgeDist), AdjGrpAddr, i);
					getVarE(PTKEY_A, Ref(PtGrpKey), AdjGrpAddr, i);

					nOfEdgeExpended++;

					//cout<<"("<<Ni<<","<<Nj<<") EdgeDist:"<<EdgeDist;
					if (PtGrpKey == -1) {
						//cout<<"No POI existed on Edge where Q located."<<endl;
					}
					else {
						getFixedF(SIZE_P, Ref(PtNumOnEdge), PtGrpKey);
						//cout<<" PtNum on Edge where Q located:"<<PtNumOnEdge<<endl;
						//Notice the order POI visited on edge
						for (int j = 0; j<PtNumOnEdge; j++) {
							nOfPOIVisited++;
							getVarE(PT_DIST, Ref(PtDist), PtGrpKey, j);

							getVarE(PT_KWD, tempKwd, PtGrpKey, j);
							int num = sizeof(tempKwd) / sizeof(int);
							for (int loop = 0; loop<sizeof(tempKwd); loop + sizeof(int)) {
								int temp;
								memcpy(&temp, tempKwd + loop, sizeof(int));
								tempKwds.insert(temp);
							}

							if (LcontianRIKwd(tempKwds, Q.kwd)) {
								POI tmp;
								if (nodei<nodej) {
									tmp.dist_toquery = iterhV->second.iDisToQuery + PtDist;
								}
								else {
									tmp.dist_toquery = iterhV->second.jDisToQuery + EdgeDist - PtDist;
								}
								getVarE(PT_ATTRIBUTE, tmp.attr, PtGrpKey, j);
								getVarE(PT_P, &tmp.poid, PtGrpKey, j);

								if (tmp.dist_toquery <= Q.distCnst) {
									updatecSDominate(tmp, cS);
								}
							}// endif
						}

					}// endelse
				}
			}
		}
	}

	for (int i = 0; i<cS.size(); i++) {
		rS.push_back(cS[i].poid);
		//cout<<cS[i].poid<<" ";
	}


}

void queryAlgorithm(const char* fileprefix) {
	initialQuery(fileprefix); //---------------------**********
	if (algorithmId == EA) { // EA
		EABasedAlgorithm(Q);
	}
	else if (algorithmId == EGBU) { // EGBU	
		EGBUAlgorithm(Q);
	}
	else { // EGTD
		EGTDAlgorithm(Q);
	}
}

int main(int argc, char *argv[]) {
	string configFileName = "config.prop";
	ConfigType cr(configFileName, argc, argv);

	//const char *indexFile = cr.getIndexFileName().c_str();
	const char *indexFile = "F:\\experiment\\index";
	algorithmId = EA;
	OpenDiskComm(indexFile, 128, algorithmId);
	queryAlgorithm(indexFile);
	CloseDiskComm();
	return 0;

}