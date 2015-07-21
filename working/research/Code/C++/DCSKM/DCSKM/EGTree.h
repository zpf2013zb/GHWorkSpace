#ifndef _EGTREE_H_
#define _EGTREE_H_

#include<stdio.h>
#include<metis.h>
#include<vector>
#include<stdlib.h>
#include<memory.h>
#include<unordered_map>
#include<map>
#include<set>
#include<deque>
#include<stack>
#include<algorithm>
//#include<netshare.h>

#include "btree.h"
#include "utility.h"
#include "netshare.h"
#include "ConfigType.h"
#include "KeywordsGenerator.h"
#include <random>
#include <bitset>

#include<sys/time.h>
using namespace std;

/********************************PreDefinition************************************/
// MACRO for timing
struct timeval tv;
long long ts, te;
#define TIME_TICK_START gettimeofday( &tv, NULL ); ts = tv.tv_sec * 100000 + tv.tv_usec / 10;
#define TIME_TICK_END gettimeofday( &tv, NULL ); te = tv.tv_sec * 100000 + tv.tv_usec / 10;
#define TIME_TICK_PRINT(T) printf("%s RESULT: %lld (0.01MS)\r\n", (#T), te - ts );
// offset 
#define _FILE_OFFSET_BITS 64
// set all edge weight to 1(unweighted graph)
#define ADJWEIGHT_SET_TO_ALL_ONE true
// we assume edge weight is integer, thus (input edge) * WEIGHT_INFLATE_FACTOR = (our edge weight)
#define WEIGHT_INFLATE_FACTOR 100000
// egtree fanout
#define PARTITION_PART 4
// egtree leaf node capacity = tau(in paper)
#define LEAF_CAP 32
#define ATTRIBUTE_DIMENSION 6
// egtree I/O file
#define FILE_NODE "road.node" //input nodes and edges
#define FILE_EDGE "road.edge"
#define FILE_NODES_GTREE_PATH "egt.paths" //output egtree information
#define FILE_GTREE 			  "egt.egtree"
#define FILE_ONTREE_MIND	  "egt.minds"
#define FILE_PART "egt.partition"


/********************************DataStructure************************************/
//-----basic infromation of gtree----------------------
typedef struct{
	double x,y;
	vector<int> adjnodes;
	vector<int> adjweight;
	bool isborder;
	vector<int> egtreepath; // this is used to do sub-graph locating
}Node; //##########思考Node是否需要记录那么多数据并且思考Edge结构体

typedef struct{
//--------------basic structure of gtree
	vector<int> borders; // the borders of this node
	vector<int> children;
	bool isleaf;
	vector<int> leafnodes;
	int father;
// ----- min dis -----
	vector<int> union_borders; // for non leaf node, merging borders of children
	vector<int> mind; // min dis, row by row of union_borders #########压缩处理，下面这些可能不需要，而且需要思考添加Skyline以及关键字信息
// ----- for pre query init, OCCURENCE LIST in paper -----
	//vector<int> nonleafinvlist;
	//vector<int> leafinvlist;
	//vector<int> up_pos;
	//vector<int> current_pos;
// --------------extend structure of egtree for EGBU
	vector<int> union_kwd;
	float attrBound[ATTRIBUTE_DIMENSION][2];
	int pterToPF; //pointer to point file
// --------------extend structure of egtree for EGTD



}TreeNode;

int noe; // number of edges
int nLeafNode;
vector<Node> Nodes;
extern vector<TreeNode> EGTree;


// init status struct
typedef struct{
	int tnid; // tree node id
	set<int> nset; // node set
}Status;

//---------------extend structure of egtree-------------------
typedef struct {
	int nid;
	int ptid;
	float edgeLth;
	vector<int> kwds;
	float attrBound[ATTRIBUTE_DIMENSION][2];
	int pointer;
}EAdjFEntry;

typedef struct {
	int poid;
    float distTV;
	int attr[ATTRIBUTE_DIMENSION];
    //unsigned long long vct;//Vector of keywords denoted by 64-bit
	vector<int> kwd;
} EPotFEntry;


// use for metis
// idx_t = int64_t / real_t = double
idx_t nvtxs; // |vertices|
idx_t ncon; // number of weight per vertex #####描述的有问题？
idx_t* xadj; // array of adjacency of indices
idx_t* adjncy; // array of adjacency nodes
idx_t* vwgt; // array of weight of nodes
idx_t* adjwgt; // array of weight of edges in adjncy
idx_t nparts; // number of parts to partition
idx_t objval; // edge cut for partitioning solution
idx_t* part; // array of partition vector

/********************************Function************************************/
//-------------basic function of gtree
void options_setting(); // METIS setting options
void init_input(int nOfNode, EdgeMapType EdgeMap); // load nodes and edges 
void data_transform_init( set<int> &nset ); // transform original data format to that suitable for METIS
void init(int nOfNode, EdgeMapType EdgeMap); // combining two steps above
void finalize(); // free space
unordered_map<int,int> graph_partition( set<int> &nset ); // graph partition
void build(); // egtree construction
void egtree_save(); // dump gtree index to file
void egtree_load();// load gtree index from file
vector<int> dijkstra_candidate( int s, vector<int> &cands, vector<Node> &graph ); // dijkstra search, used for single-source shortest path search WITHIN one gtree leaf node!
void hierarchy_shortest_path_calculation(); // calculate the distance matrix
void hierarchy_shortest_path_save(); // dump distance matrix into file
void hierarchy_shortest_path_load(); // load distance matrix from file
int mainFunction(int nOfNode, EdgeMapType EdgeMap); // main function

//---------------extend function of egtree
void makeEPtFiles(FILE *ptFile,char* treefile); // construct the extend point file
void makeEAdjListFiles(FILE *alFile); // construct the extend adjacentList file
void BuildEBinaryStorage(const char* fileprefix); // construct the extend binary storage
// 
#endif