#include "StdAfx.h"
#include "EGTree.h"

// METIS setting options
void options_setting() {
	METIS_SetDefaultOptions(options);
	options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY; // _RB
	options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT; // _VOL
	options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM; // _RM
	options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_RANDOM; // _GROW _EDGE _NODE
	options[METIS_OPTION_RTYPE] = METIS_RTYPE_FM; // _GREEDY _SEP2SIDED _SEP1SIDED
												  // options[METIS_OPTION_NCUTS] = 1;
												  // options[METIS_OPTION_NITER] = 10;
												  /* balance factor, used to be 500 */
	options[METIS_OPTION_UFACTOR] = 500;
	// options[METIS_OPTION_MINCONN];
	options[METIS_OPTION_CONTIG] = 1;
	// options[METIS_OPTION_SEED];
	options[METIS_OPTION_NUMBERING] = 0;
	// options[METIS_OPTION_DBGLVL] = 0;
}
// read the data from EdgeMap instead of file
void init_input(int nOfNode, EdgeMapType EdgeMap) {
	nLeafNode = 0;
	// process node information
	printf("PROCESSING NODE...");
	// notable that vertex id in node is from 0
	for (int i = 0; i<nOfNode; i++) {
		Node node;
		node.isborder = false;
		Nodes.push_back(node);
	}
	printf("COMPLETE. NODE_COUNT=%d\n", (int)Nodes.size());

	// load edge
	printf("PROCESSING EDGE...");
	//fin = fopen(FILE_EDGE, "r");
	int eid;
	int snid, enid;
	double weight;
	int iweight;
	noe = 0;
	EdgeMapType::iterator iter = EdgeMap.begin();
	for (; iter != EdgeMap.end(); iter++) {
		edge* e = iter->second;
		noe++;
		snid = e->Ni;
		enid = e->Nj;
		weight = e->dist;

		iweight = (int)(weight * WEIGHT_INFLATE_FACTOR);
		Nodes[snid].adjnodes.push_back(enid);
		Nodes[snid].adjweight.push_back(iweight);
		Nodes[enid].adjnodes.push_back(snid);
		Nodes[enid].adjweight.push_back(iweight);
	}
	//fclose(fin);
	printf("COMPLETE.\n");
}

/*
void init_input(int nOfNode, const EdgeMapType EdgeMap){
FILE *fin;

// load node
printf("LOADING NODE...");
fin = fopen(FILE_NODE,"r");
int nid;
double x,y;
while( fscanf(fin, "%d %lf %lf", &nid, &x, &y) == 3 ){
Node node = { x, y };
Nodes.push_back(node);
}
fclose(fin);
printf("COMPLETE. NODE_COUNT=%d\n", (int)Nodes.size());

// load edge
printf("LOADING EDGE...");
fin = fopen(FILE_EDGE, "r");
int eid;
int snid, enid;
double weight;
int iweight;
noe = 0;
while( fscanf(fin,"%d %d %d %lf", &eid, &snid, &enid, &weight ) == 4 ){
noe ++;
iweight = (int) (weight * WEIGHT_INFLATE_FACTOR);
Nodes[snid].adjnodes.push_back( enid );
Nodes[snid].adjweight.push_back( iweight );
Nodes[enid].adjnodes.push_back( snid );
Nodes[enid].adjweight.push_back( iweight );
}
fclose(fin);
printf("COMPLETE.\n");
}
*/

// transform original data format to that suitable for METIS
void data_transform_init(set<int> &nset) {
	// nvtxs, ncon
	nvtxs = nset.size();
	ncon = 1;

	xadj = new idx_t[nset.size() + 1];
	adjncy = new idx_t[noe * 2];
	adjwgt = new idx_t[noe * 2];


	int xadj_pos = 1;
	int xadj_accum = 0;
	int adjncy_pos = 0;

	// xadj, adjncy, adjwgt
	unordered_map<int, int> nodemap;
	nodemap.clear();

	xadj[0] = 0;
	int i = 0;
	for (set<int>::iterator it = nset.begin(); it != nset.end(); it++, i++) {
		// init node map
		nodemap[*it] = i;

		int nid = *it;
		int fanout = Nodes[nid].adjnodes.size();
		for (int j = 0; j < fanout; j++) {
			int enid = Nodes[nid].adjnodes[j];
			// ensure edges within
			if (nset.find(enid) != nset.end()) {
				// xadj_accum used to record the start and end position of adjacent nodes for current visited nodes
				xadj_accum++;

				adjncy[adjncy_pos] = enid;
				adjwgt[adjncy_pos] = Nodes[nid].adjweight[j];
				adjncy_pos++;
			}
		}
		xadj[xadj_pos++] = xadj_accum;
	}

	// adjust nodes number started by 0  ###########这部分以及下面部分要修改，为何要做映射，为什么权重改为1？(可以继续使用)
	for (int i = 0; i < adjncy_pos; i++) {
		adjncy[i] = nodemap[adjncy[i]];
	}

	// adjwgt -> 1
	if (ADJWEIGHT_SET_TO_ALL_ONE) {
		for (int i = 0; i < adjncy_pos; i++) {
			adjwgt[i] = 1;
		}
	}

	// nparts
	nparts = PARTITION_PART;

	// part
	part = new idx_t[nset.size()];
}

void init(int nOfNode, const EdgeMapType EdgeMap) {
	init_input(nOfNode, EdgeMap);
	options_setting();
}

void finalize() {
	delete xadj;
	delete adjncy;
	delete adjwgt;
	delete part;
}

// graph partition
// input: nset = a set of node id
// output: <node, node belong to partition id>
unordered_map<int, int> graph_partition(set<int> &nset) {
	unordered_map<int, int> result;

	// transform data to metis
	data_transform_init(nset);

	// partition, result -> part
	// k way partition
	int returnstate;
	returnstate = METIS_PartGraphKway(
		&nvtxs,
		&ncon,
		xadj,
		adjncy,
		NULL,
		NULL,
		adjwgt,
		&nparts,
		NULL,
		NULL,
		options,
		&objval,
		part
		);
	// return result state
	if (METIS_ERROR_INPUT == returnstate) {
		printf("Input Error!");
		return result;
	}
	else if (METIS_ERROR_MEMORY == returnstate) {
		printf("Memory Error!");
		return result;
	}
	else if (METIS_ERROR == returnstate) {
		printf("Other Type Error!");
		return result;
	}
	else {

	}

	// push to result
	result.clear();
	int i = 0;
	// 又将结果变回来了，在上面采用nodemap可能是为了满足一定的条件，part中的结果是什么呢？子部分用什么来表示呢？
	for (set<int>::iterator it = nset.begin(); it != nset.end(); it++, i++) {
		result[*it] = part[i];
	}

	// finalize
	finalize();

	return result;
}

// egtree construction
void build(EdgeMapType EdgeMap) {
	// init root
	TreeNode root;
	root.isleaf = false;
	root.father = -1;
	EGTree.push_back(root);

	// init stack
	stack<Status> buildstack;
	Status rootstatus;
	rootstatus.tnid = 0;
	rootstatus.nset.clear();
	// careful the index of node start from 1
	for (int i = 1; i <= Nodes.size(); i++) {
		rootstatus.nset.insert(i);
	}
	buildstack.push(rootstatus);

	// start to build
	unordered_map<int, int> presult;
	set<int> childset[PARTITION_PART];


	while (buildstack.size() > 0) {
		// pop top
		Status current = buildstack.top();
		buildstack.pop();

		// update egtreepath
		for (set<int>::iterator it = current.nset.begin(); it != current.nset.end(); it++) {
			Nodes[*it].egtreepath.push_back(current.tnid);
		}

		// check cardinality
		if (current.nset.size() <= LEAF_CAP) {
			// build leaf node
			nLeafNode++;
			EGTree[current.tnid].isleaf = true;
			EGTree[current.tnid].leafnodes.clear();
			for (set<int>::iterator it = current.nset.begin(); it != current.nset.end(); it++) {
				EGTree[current.tnid].leafnodes.push_back(*it);
				//---------修改partID
				//partID[*it] = current.tnid;

				//------------------------------计算整个区域的上下界和关键字
				//---------------------M--replace the adjList with Nodes--------
				for (int i = 0; i<Nodes[*it].adjnodes.size(); i++) {
					int Nj = Nodes[*it].adjnodes[i];	// Nk can be smaller or greater than Ni !!!
					edge* e = EdgeMap[getKey(*it, Nj)];

					if (it == current.nset.begin()) { //直接初始化
						for (int j = 0; j<ATTRIBUTE_DIMENSION; j++) {
							EGTree[current.tnid].attrBound[j][0] = e->attrBound[j][0];
							EGTree[current.tnid].attrBound[j][1] = e->attrBound[j][1];
						}
						copy(e->kwds.begin(), e->kwds.end(), EGTree[current.tnid].union_kwd.begin());
					}
					else {
						for (int j = 0; j<ATTRIBUTE_DIMENSION; j++) {
							if (EGTree[current.tnid].attrBound[j][0]>e->attrBound[j][0]) {
								EGTree[current.tnid].attrBound[j][0] = e->attrBound[j][0];
							}
							if (EGTree[current.tnid].attrBound[j][1]<e->attrBound[j][1]) {
								EGTree[current.tnid].attrBound[j][1] = e->attrBound[j][1];
							}
						}
						set_union(e->kwds.begin(), e->kwds.end(), EGTree[current.tnid].union_kwd.begin(), EGTree[current.tnid].union_kwd.end(), EGTree[current.tnid].union_kwd.begin());
					}

				}//endfor i

			}//endfor iteration
			continue;
		}

		// partition
		//		printf("PARTITIONING...NID=%d...SIZE=%d...", current.tnid, (int)current.nset.size() );
		presult = graph_partition(current.nset);
		//		printf("COMPLETE.\n");

		// construct child node set
		for (int i = 0; i < PARTITION_PART; i++) {
			childset[i].clear();
		}
		int slot;
		// put the nodes into corresponding sub-partition(slot)
		for (set<int>::iterator it = current.nset.begin(); it != current.nset.end(); it++) {
			slot = presult[*it];
			childset[slot].insert(*it);
		}

		// generate child tree nodes
		int childpos;
		for (int i = 0; i < PARTITION_PART; i++) {
			TreeNode tnode;
			tnode.isleaf = false;
			tnode.father = current.tnid;

			// insert to EGTree first
			EGTree.push_back(tnode);
			childpos = EGTree.size() - 1;
			EGTree[current.tnid].children.push_back(childpos);

			// calculate border nodes
			EGTree[childpos].borders.clear();
			for (set<int>::iterator it = childset[i].begin(); it != childset[i].end(); it++) {

				bool isborder = false;
				for (int j = 0; j < Nodes[*it].adjnodes.size(); j++) {
					if (childset[i].find(Nodes[*it].adjnodes[j]) == childset[i].end()) {
						isborder = true;
						break;
					}
				}
				if (isborder) {
					EGTree[childpos].borders.push_back(*it);
					// update globally
					Nodes[*it].isborder = true;
				}
			}

			// add to stack
			Status ongoingstatus;
			ongoingstatus.tnid = childpos;
			ongoingstatus.nset = childset[i];
			buildstack.push(ongoingstatus);

		}

	}

}

// dump EGTree index to file
void egtree_save(const char* filename) {
	// FILE_GTREE
	printf("making egtreeFile\n");
	FILE *fout = fopen(filename, "wb");
	int *buf = new int[Nodes.size()];
	float *buff = new float[Nodes.size()];
	for (int i = 0; i < EGTree.size(); i++) {
		// borders
		int count_borders = EGTree[i].borders.size();
		fwrite(&count_borders, sizeof(int), 1, fout);
		copy(EGTree[i].borders.begin(), EGTree[i].borders.end(), buf);
		fwrite(buf, sizeof(int), count_borders, fout);
		// children
		int count_children = EGTree[i].children.size();
		fwrite(&count_children, sizeof(int), 1, fout);
		copy(EGTree[i].children.begin(), EGTree[i].children.end(), buf);
		fwrite(buf, sizeof(int), count_children, fout);
		// isleaf
		fwrite(&EGTree[i].isleaf, sizeof(bool), 1, fout);
		// leafnodes
		int count_leafnodes = EGTree[i].leafnodes.size();
		fwrite(&count_leafnodes, sizeof(int), 1, fout);
		copy(EGTree[i].leafnodes.begin(), EGTree[i].leafnodes.end(), buf);
		fwrite(buf, sizeof(int), count_leafnodes, fout);
		// father
		fwrite(&EGTree[i].father, sizeof(int), 1, fout);

		//**************************************************
		//保存EGBU和EGTD信息

		// union_border
		int count_unionBorder = EGTree[i].union_borders.size();
		fwrite(&count_unionBorder, sizeof(int), 1, fout);
		copy(EGTree[i].union_borders.begin(), EGTree[i].union_borders.end(), buf);
		fwrite(buf, sizeof(int), count_unionBorder, fout);
		// mind
		int count_mind = EGTree[i].mind.size();
		fwrite(&count_mind, sizeof(float), 1, fout);
		copy(EGTree[i].mind.begin(), EGTree[i].mind.end(), buff);
		fwrite(buff, sizeof(float), count_mind, fout);
		// union_kwd
		int count_unionKwd = EGTree[i].union_kwd.size();
		fwrite(&count_unionKwd, sizeof(int), 1, fout);
		copy(EGTree[i].union_kwd.begin(), EGTree[i].union_kwd.end(), buf);
		fwrite(buf, sizeof(int), count_unionKwd, fout);
		// attrBound
		//int count_leafnodes = EGTree[i].leafnodes.size();
		//fwrite(&count_leafnodes, sizeof(int), 1, fout);
		//copy(EGTree[i].leafnodes.begin(), EGTree[i].leafnodes.end(), buf);
		fwrite(EGTree[i].attrBound, sizeof(float), ATTRIBUTE_DIMENSION * 2, fout);
		//perToPF
		//int count_leafnodes = EGTree[i].leafnodes.size();
		//fwrite(&count_leafnodes, sizeof(int), 1, fout);
		//copy(EGTree[i].leafnodes.begin(), EGTree[i].leafnodes.end(), buf);
		fwrite(&EGTree[i].pterToPF, sizeof(int), 1, fout);
		// editKwd;
		int count_editKwd = EGTree[i].editKwd.size();
		fwrite(&count_editKwd, sizeof(int), 1, fout);
		int tempKwd;
		set<set<int>> ::iterator itEdit = EGTree[i].editKwd.begin();
		for (; itEdit != EGTree[i].editKwd.end(); itEdit++) {
			set<int> temp = *itEdit;
			tempKwd = temp.size();
			fwrite(&tempKwd, sizeof(int), 1, fout);
			set<int> ::iterator it = temp.begin();
			int kwd;
			for (; it != temp.end(); it++) {
				kwd = *it;
				fwrite(&kwd, sizeof(int), 1, fout);
			}
		}
		copy(EGTree[i].editKwd.begin(), EGTree[i].editKwd.end(), buf);
		fwrite(buf, sizeof(int), count_editKwd, fout);
		// locSky
		int count_locSky = EGTree[i].locSky.size();
		fwrite(&count_locSky, sizeof(int), 1, fout);
		map<int, locSkyline> ::iterator itLoc = EGTree[i].locSky.begin();
		for (; itLoc != EGTree[i].locSky.end(); itLoc++) {
			int key = itLoc->first;
			locSkyline value = itLoc->second;
			fwrite(&key, sizeof(int), 1, fout);
			fwrite(value.blockID, sizeof(int), ATTRIBUTE_DIMENSION, fout);
			fwrite(value.skylineBtBound, sizeof(float), ATTRIBUTE_DIMENSION, fout);
			fwrite(value.skylineUpBound, sizeof(float), ATTRIBUTE_DIMENSION, fout);
		}
		copy(EGTree[i].locSky.begin(), EGTree[i].locSky.end(), buf);
		fwrite(buf, sizeof(int), count_locSky, fout);
	}
	fclose(fout);

	//**************************************************
	//保存EGBU和EGTD信息




	// FILE_NODES_GTREE_PATH
	fout = fopen(FILE_NODES_GTREE_PATH, "wb");
	for (int i = 0; i < Nodes.size(); i++) {
		int count = Nodes[i].egtreepath.size();
		fwrite(&count, sizeof(int), 1, fout);
		copy(Nodes[i].egtreepath.begin(), Nodes[i].egtreepath.end(), buf);
		fwrite(buf, sizeof(int), count, fout);
	}
	fclose(fout);
	delete[] buf;
}

// load EGTree index from file
void egtree_load(const char* filename, vector<TreeNode>& EGTree) {
	// FILE_GTREE
	printf("loading egtreeFile\n");
	FILE *fin = fopen(filename, "rb");
	int *buf = new int[Nodes.size()];
	float *buff = new float[Nodes.size()];
	int count_borders, count_children, count_leafnodes;
	bool isleaf;
	int father;

	// clear EGTree
	EGTree.clear();

	while (fread(&count_borders, sizeof(int), 1, fin)) {
		TreeNode tn;
		// borders
		tn.borders.clear();
		fread(buf, sizeof(int), count_borders, fin);
		for (int i = 0; i < count_borders; i++) {
			tn.borders.push_back(buf[i]);
		}
		// children
		fread(&count_children, sizeof(int), 1, fin);
		fread(buf, sizeof(int), count_children, fin);
		for (int i = 0; i < count_children; i++) {
			tn.children.push_back(buf[i]);
		}
		// isleaf
		fread(&isleaf, sizeof(bool), 1, fin);
		tn.isleaf = isleaf;
		// leafnodes
		fread(&count_leafnodes, sizeof(int), 1, fin);
		fread(buf, sizeof(int), count_leafnodes, fin);
		for (int i = 0; i < count_leafnodes; i++) {
			tn.leafnodes.push_back(buf[i]);
		}
		// father
		fread(&father, sizeof(int), 1, fin);
		tn.father = father;

		//**************************************************
		//加载EGBU和EGTD信息

		// union_border
		int count_unionBorder;
		fread(&count_unionBorder, sizeof(int), 1, fin);
		fread(buf, sizeof(int), count_unionBorder, fin);
		for (int i = 0; i < count_unionBorder; i++) {
			tn.union_borders.push_back(buf[i]);
		}
		// mind
		int count_mind;
		fread(&count_mind, sizeof(float), 1, fin);
		//copy(EGTree[i].leafnodes.begin(), EGTree[i].leafnodes.end(), buf);
		fread(buff, sizeof(float), count_mind, fin);
		for (int i = 0; i < count_mind; i++) {
			tn.mind.push_back(buff[i]);
		}
		// union_kwd
		int count_unionKwd;
		fread(&count_unionKwd, sizeof(int), 1, fin);
		//copy(EGTree[i].leafnodes.begin(), EGTree[i].leafnodes.end(), buf);
		fread(buf, sizeof(int), count_unionKwd, fin);
		for (int i = 0; i < count_unionKwd; i++) {
			tn.union_kwd.insert(buf[i]);
		}
		// attrBound
		//int count_leafnodes = EGTree[i].leafnodes.size();
		fread(tn.attrBound, sizeof(float), ATTRIBUTE_DIMENSION * 2, fin);
		//copy(EGTree[i].leafnodes.begin(), EGTree[i].leafnodes.end(), buf);
		//fwrite(buf, sizeof(int), count_leafnodes, fin);

		//perToPF
		//int count_leafnodes = EGTree[i].leafnodes.size();
		fread(&tn.pterToPF, sizeof(int), 1, fin);
		//copy(EGTree[i].leafnodes.begin(), EGTree[i].leafnodes.end(), buf);
		//fread(buf, sizeof(int), count_leafnodes, fin);

		// editKwd;
		int count_editKwd;
		fread(&count_editKwd, sizeof(int), 1, fin);
		int tempKwd;
		for (int i = 0; i < count_editKwd; i++) {
			fread(&tempKwd, sizeof(int), 1, fin);
			set<int> tmp;
			int kwd;
			for (int t = 0; t < tempKwd; t++) {
				fread(&kwd, sizeof(int), 1, fin);
				tmp.insert(kwd);
			}
			tn.editKwd.insert(tmp);
		}

		// locSky
		int count_locSky;
		fread(&count_locSky, sizeof(int), 1, fin);
		for (int k = 0; k < count_locSky; k++) {
			int key;
			locSkyline ls;
			fread(&key, sizeof(int), 1, fin);
			fread(ls.blockID, sizeof(int), ATTRIBUTE_DIMENSION, fin);
			fread(&ls.skylineBtBound, sizeof(float), ATTRIBUTE_DIMENSION, fin);
			fread(&ls.skylineUpBound, sizeof(float), ATTRIBUTE_DIMENSION, fin);
		}
		EGTree.push_back(tn);
	}
	fclose(fin);
	delete[] buf;
	delete[] buff;
	//--------------------------M-- no use---------------------
	/*
	//**************************************************
	//加载EGBU和EGTD信息


	// FILE_NODES_GTREE_PATH
	int count;
	fin = fopen( FILE_NODES_GTREE_PATH, "rb" );
	int pos = 0;
	while( fread( &count, sizeof(int), 1, fin ) ){
	fread( buf, sizeof(int), count, fin );
	// clear egtreepath
	Nodes[pos].egtreepath.clear();
	for ( int i = 0; i < count; i++ ){
	Nodes[pos].egtreepath.push_back( buf[i] );
	}
	// pos increase
	pos ++;
	}
	fclose(fin);
	delete[] buf;
	*/
}

// dijkstra search, used for single-source shortest path search WITHIN one EGTree leaf node!
// input: s = source node
//        cands = candidate node list
//        graph = search graph(this can be set to subgraph)
vector<int> dijkstra_candidate(int s, vector<int> &cands, vector<Node> &graph) {
	// init
	set<int> todo;
	todo.clear();
	todo.insert(cands.begin(), cands.end());

	unordered_map<int, int> result;
	result.clear();
	set<int> visited;
	visited.clear();
	unordered_map<int, int> q;
	q.clear();
	q[s] = 0;

	// start
	int min, minpos, adjnode, weight;
	while (!todo.empty() && !q.empty()) {
		min = -1;
		for (unordered_map<int, int>::iterator it = q.begin(); it != q.end(); it++) {
			if (min == -1) {
				minpos = it->first;
				min = it->second;
			}
			else {
				if (it->second < min) {
					min = it->second;
					minpos = it->first;
				}
			}
		}
		
		// put min to result, add to visited
		result[minpos] = min;
		visited.insert(minpos);
		q.erase(minpos);

		if (todo.find(minpos) != todo.end()) {
			todo.erase(minpos);
		}

		// expand
		for (int i = 0; i < graph[minpos].adjnodes.size(); i++) {
			adjnode = graph[minpos].adjnodes[i];
			if (visited.find(adjnode) != visited.end()) {
				continue;
			}
			weight = graph[minpos].adjweight[i];

			if (q.find(adjnode) != q.end()) {
				if (min + weight < q[adjnode]) {
					q[adjnode] = min + weight;
				}
			}
			else {
				q[adjnode] = min + weight;
			}

		}
	}

	// output
	vector<int> output;
	for (int i = 0; i < cands.size(); i++) {
		output.push_back(result[cands[i]]);
	}

	// return
	return output;
}

//test is be dominate
// smaller is better
int rdominatel(InerNode left, InerNode right) {
	int size = ATTRIBUTE_DIMENSION;
	bool lDr = true;
	bool rDl = true;
	int equal = 0;
	for (int i = 0; i < size; i++) {
		if (left.attr[i] < right.attr[i]) rDl = false;
		if (left.attr[i] > right.attr[i]) lDr = false;
		if (left.attr[i] == right.attr[i]) equal++;
	}
	if (lDr && (equal != size)) { //dominate local
		return -1;
	}

	if (rDl) { //
		return 1;
	}

	if (!rDl && !lDr) {
		return 0;
	}

}
// must in descending order of size
bool sortBySize(InerNode left, InerNode right) {
	if (left.kwd.size()>right.kwd.size()) return true;
	else return false;
}

bool sortByKSize(set<int> left, set<int> right) {
	if (left.size()>right.size()) return true;
	else return false;
}

int editDistanceRTL(set<int> kwd1, set<int> kwd2) {
	int n = kwd1.size();
	int m = kwd2.size();
	set<int> difference;
	// contain in kwd1 not in kwd2, n>=m
	if (n > m) {
		set_difference(kwd1.begin(), kwd1.end(), kwd2.begin(), kwd2.end(), difference.begin());
	}
	else {
		set_difference(kwd2.begin(), kwd2.end(), kwd1.begin(), kwd1.end(), difference.begin());
	}
	return difference.size();
}


void handleKwdAttr(int tn, vector<InerNode> nodeInerNode) {
	//first handle kwd
	vector<InerNode> temp(nodeInerNode);
	sort(temp.begin(), temp.end(), sortBySize);
	vector<InerNode>::iterator itiT = temp.begin();
	vector<InerNode>::iterator itjT;

	//------------------------注意，这里面很可能会有问题，访问同时删除---------------------
	for (; itiT != temp.end();) {
		InerNode iNi = *itiT;
		for (itjT = itiT + 1; itjT != temp.end();) {
			InerNode iNj = *itjT;
			if (editDistanceRTL(iNi.kwd, iNj.kwd)<iNi.kwd.size()*edDis) {
				itjT = temp.erase(itjT);
			}
			else { // l can represent r
				   // remove r 
				itjT++;
			}
		}
		itiT++;
	}
	// 保存kwd信息

	vector<InerNode>::iterator itKwd = temp.begin();
	for (; itKwd != temp.end(); itKwd++) {
		set<int> tmpKwd = itKwd->kwd;
		//kwdC kc;
		//kc.kwdCombination = tmpKwd;
		EGTree[tn].editKwd.insert(tmpKwd);
	}

	//then handle attr

	vector<InerNode>::iterator iti = nodeInerNode.begin();
	vector<InerNode>::iterator itj;
	for (; iti != nodeInerNode.end(); ) {
		InerNode iNi = *iti;
		for (itj = iti + 1; itj != nodeInerNode.end(); ) {
			InerNode iNj = *itj;
			if (rdominatel(iNi, iNj) == -1) {
				itj = nodeInerNode.erase(itj);
			}
			else if (rdominatel(iNi, iNj) == 0) { // not dominate each other
												  // remove r 
				itj++;
			}
			else {
				iti = nodeInerNode.erase(iti);
			}
		}
	}

	// 计算block
	for (int i = 0; i < nodeInerNode.size(); i++) {
		InerNode inode = nodeInerNode[i];
		int num[ATTRIBUTE_DIMENSION];
		int mapKey = 0;
		for (int j = 0; j < ATTRIBUTE_DIMENSION; j++) {
			for (int t = 1; t <= splitBlock; t++) {
				float thre = 1.0*t / splitBlock;
				if (inode.attr[j] <= thre) {
					num[j] = t;
					break;
				}
			}
			mapKey = mapKey * 10 + num[j];
		}
		if (EGTree[tn].locSky.find(mapKey) != EGTree[tn].locSky.end()) {//已经有了
			locSkyline ls = EGTree[tn].locSky[mapKey];
			for (int t = 0; t < ATTRIBUTE_DIMENSION; t++) {
				if (ls.skylineBtBound[t] > num[t]) ls.skylineBtBound[t] = inode.attr[t];
				if (ls.skylineUpBound[t] < num[t]) ls.skylineUpBound[t] = inode.attr[t];
			}
		}
		else { //不存在,新建
			locSkyline ls;
			for (int t = 0; t < ATTRIBUTE_DIMENSION; t++) {
				ls.skylineBtBound[t] = inode.attr[t];
				ls.skylineUpBound[t] = inode.attr[t];
				ls.blockID[t] = num[t];
			}
			EGTree[tn].locSky[mapKey] = ls;
		}
	}

}

// combine the information of child node to tn
void handleINKwdAddr(int tn) {
	set<set<int>> unionEditKwd; //used to record the 
	map<int, locSkyline> unionLocSky;
	// if there is problem when combine????
	for (int k = 0; k < EGTree[tn].children.size(); k++) {
		int cid = EGTree[tn].children[k];
		set_union(EGTree[cid].editKwd.begin(), EGTree[cid].editKwd.end(), unionEditKwd.begin(), unionEditKwd.end(), unionEditKwd.begin());
		//set_union(EGTree[cid].locSky.begin(), EGTree[cid].locSky.end(), unionLocSky.begin(), unionLocSky.end(), unionLocSky.begin());
	}
	// handle kwd information 
	sort(unionEditKwd.begin(), unionEditKwd.end(), sortByKSize);
	set<set<int>>::iterator itiT = unionEditKwd.begin();
	set<set<int>>::iterator itjT;

	//------------------------注意，这里面很可能会有问题，访问同时删除---------------------
	for (; itiT != unionEditKwd.end();) {
		set<int> iNi = *itiT;
		for (itjT = ++itiT; itjT != unionEditKwd.end();) {
			itiT--; //--------------note this change-------------
			set<int> iNj = *itjT;
			if (editDistanceRTL(iNi, iNj)<iNi.size()*edDis) {
				itjT = unionEditKwd.erase(itjT);
			}
			else { // l can represent r
				   // remove r 
				itjT++;
			}
		}
		itiT++;
	}
	// 保存kwd信息

	set<set<int>>::iterator itKwd = unionEditKwd.begin();
	for (; itKwd != unionEditKwd.end(); itKwd++) {
		set<int> tmpKwd = *itKwd;
		//kwdC kc;
		//kc.kwdCombination = tmpKwd;
		EGTree[tn].editKwd.insert(tmpKwd);
	}

	// handle attr information
	for (int k = 0; k < EGTree[tn].children.size(); k++) {
		int cid = EGTree[tn].children[k];
		//set_union(EGTree[cid].editKwd.begin(), EGTree[cid].editKwd.end(), unionEditKwd.begin(), unionEditKwd.end(), unionEditKwd.begin());
		if (k == 0) {
			set_union(EGTree[cid].locSky.begin(), EGTree[cid].locSky.end(), unionLocSky.begin(), unionLocSky.end(), unionLocSky.begin());
		}
		else {
			map<int, locSkyline>::iterator iterLocSkyline = EGTree[cid].locSky.begin();
			for (; iterLocSkyline != EGTree[cid].locSky.end(); iterLocSkyline++) {
				int key = iterLocSkyline->first;
				locSkyline value = iterLocSkyline->second;
				if (EGTree[tn].locSky.find(key) != EGTree[tn].locSky.end()) {//已经存在
					for (int t = 0; t < ATTRIBUTE_DIMENSION; t++) {
						if (EGTree[tn].locSky[key].skylineBtBound[t] > value.skylineBtBound[t]) EGTree[tn].locSky[key].skylineBtBound[t] = value.skylineBtBound[t];
						if (EGTree[tn].locSky[key].skylineUpBound[t] < value.skylineUpBound[t]) EGTree[tn].locSky[key].skylineUpBound[t] = value.skylineUpBound[t];
					}
				}
				else { //不存在直接加入进去
					EGTree[tn].locSky[key] = value;
				}
			}
		}

	}

}


// calculate the distance matrix
void hierarchy_shortest_path_calculation() {
	// level traversal
	vector< vector<int> > treenodelevel;

	vector<int> current; // record all the treenodes in current level 
	current.clear();
	current.push_back(0);
	treenodelevel.push_back(current);
	// put all the nodes into treenodelevel according to their levels
	vector<int> mid;
	while (current.size() != 0) {
		mid = current;
		current.clear();
		for (int i = 0; i < mid.size(); i++) {
			for (int j = 0; j < EGTree[mid[i]].children.size(); j++) {
				current.push_back(EGTree[mid[i]].children[j]);
			}
		}
		if (current.size() == 0) break;
		treenodelevel.push_back(current);
	}

	// bottom up calculation
	// temp graph
	vector<Node> graph;
	graph = Nodes;
	vector<int> cands;
	vector<int> result;
	unordered_map<int, unordered_map<int, int> > vertex_pairs;

	// do dijkstra
	int s, t, tn, nid, cid, weight;
	vector<int> tnodes, tweight;
	set<int> nset;
	//set<set<int>> nodeKwd;// 关键字信息等到该节点处理完成后统一处理，而关键字信息一样
	vector<InerNode> nodeInerNode;

	for (int i = treenodelevel.size() - 1; i >= 0; i--) {
		for (int j = 0; j < treenodelevel[i].size(); j++) {

			tn = treenodelevel[i][j];

			nodeInerNode.clear();
			cands.clear();

			if (EGTree[tn].isleaf) {
				//sort lefenodes and union_borders
				sort(EGTree[tn].leafnodes.begin(), EGTree[tn].leafnodes.end());
				sort(EGTree[tn].borders.begin(), EGTree[tn].borders.end());
				// cands = leafnodes
				cands = EGTree[tn].leafnodes;
				// union borders = borders;
				EGTree[tn].union_borders = EGTree[tn].borders;

				//--------------record the mind information
				vertex_pairs.clear();

				// for each border, do min dis
				//int cc = 0;

				for (int k = 0; k < EGTree[tn].union_borders.size(); k++) {
					//printf("DIJKSTRA...LEAF=%d BORDER=%d\n", tn, EGTree[tn].union_borders[k] );
					result = dijkstra_candidate(EGTree[tn].union_borders[k], cands, graph);
					//printf("DIJKSTRA...END\n");

					// save to map
					for (int p = 0; p < result.size(); p++) {
						EGTree[tn].mind.push_back(result[p]);
						vertex_pairs[EGTree[tn].union_borders[k]][cands[p]] = result[p];
					}
				}

				//---------------extend for EGTD Algorithm-----------------
				int keyID;
				for (int k = 0; k < EGTree[tn].leafnodes.size(); k++) {
					// for each leafnode we handle each adjacent edge of it
					int nodei = EGTree[tn].leafnodes[k];
					int adjSize = Nodes[nodei].adjnodes.size();
					for (int p = 0; p < adjSize; p++) {
						int nodej = Nodes[nodei].adjnodes[p];
						keyID = getKey(nodei, nodej);
						if (visitedEdgeKey.count(keyID) == 0) {// this edge has not been visited
							visitedEdgeKey.insert(keyID);
							edge *e = EdgeMap[keyID];
							vector<InerNode> inode = e->pts;
							for (int b = 0; b < inode.size(); b++) {
								InerNode tempN = inode[b];
								//set<int> poiKwd = tempN.kwd;
								//nodeKwd.insert(tempN.kwd);
								nodeInerNode.push_back(tempN);
							}

						}
					}

					//统一处理所有的InterNode信息
					handleKwdAttr(tn, nodeInerNode);
				}
			}
			else {

				nset.clear();
				for (int k = 0; k < EGTree[tn].children.size(); k++) {
					cid = EGTree[tn].children[k];
					nset.insert(EGTree[cid].borders.begin(), EGTree[cid].borders.end());


					//------------------------处理中间节点的信息,将多个孩子节点的值赋给它


					if (k == 0) { //直接初始化
						for (int l = 0; l<ATTRIBUTE_DIMENSION; l++) {
							EGTree[tn].attrBound[l][0] = EGTree[cid].attrBound[l][0];
							EGTree[tn].attrBound[l][1] = EGTree[cid].attrBound[l][1];
							//memcpy(EGTree[tn].attrBound,EGTree[cid].attrBound,sizeof(EGTree[tn].attrBound));
						}
						copy(EGTree[cid].union_kwd.begin(), EGTree[cid].union_kwd.end(), EGTree[tn].union_kwd.begin());

						//------------------extend for EGTD----------------------
						// handle kwd
						//EGTree[tn].editKwd = EGTree[cid].editKwd;

						// handle attr
						//EGTree[tn].locSky = EGTree[cid].locSky;
					}
					else {
						for (int l = 0; l<ATTRIBUTE_DIMENSION; l++) {
							if (EGTree[tn].attrBound[l][0]>EGTree[cid].attrBound[l][0]) {
								EGTree[tn].attrBound[l][0] = EGTree[cid].attrBound[l][0];
							}
							if (EGTree[tn].attrBound[l][1]<EGTree[cid].attrBound[l][1]) {
								EGTree[tn].attrBound[l][1] = EGTree[cid].attrBound[l][1];
							}
						}
						set_union(EGTree[cid].union_kwd.begin(), EGTree[cid].union_kwd.end(), EGTree[tn].union_kwd.begin(), EGTree[tn].union_kwd.end(), EGTree[tn].union_kwd.begin());
						// handle kwd 
						//handleINKwdAddr(tn);
						// handle attr
					}

				}
				//---------------extend for EGTD---------------
				handleINKwdAddr(tn);


				// union borders = cands;

				cands.clear();
				for (set<int>::iterator it = nset.begin(); it != nset.end(); it++) {
					cands.push_back(*it);
				}
				EGTree[tn].union_borders = cands;
				//sort lefenodes and union_borders
				sort(EGTree[tn].union_borders.begin(), EGTree[tn].union_borders.end());
				sort(cands.begin(), cands.end());

				//------------------------record the mind information
				vertex_pairs.clear();

				// for each border, do min dis
				//int cc = 0;

				for (int k = 0; k < EGTree[tn].union_borders.size(); k++) {
					//printf("DIJKSTRA...LEAF=%d BORDER=%d\n", tn, EGTree[tn].union_borders[k] );
					result = dijkstra_candidate(EGTree[tn].union_borders[k], cands, graph);
					//printf("DIJKSTRA...END\n");

					// save to map
					for (int p = 0; p < result.size(); p++) {
						if (k <= p) {
							EGTree[tn].mind.push_back(result[p]);
						}
						vertex_pairs[EGTree[tn].union_borders[k]][cands[p]] = result[p];
					}
				}



			}
			//------------------------------以前代码-----------------------------------
			// start to do min dis
			/*
			vertex_pairs.clear();

			// for each border, do min dis
			//int cc = 0;

			for ( int k = 0; k < EGTree[tn].union_borders.size(); k++ ){
			//printf("DIJKSTRA...LEAF=%d BORDER=%d\n", tn, EGTree[tn].union_borders[k] );
			result = dijkstra_candidate( EGTree[tn].union_borders[k], cands, graph );
			//printf("DIJKSTRA...END\n");

			// save to map
			for ( int p = 0; p < result.size(); p ++ ){
			EGTree[tn].mind.push_back( result[p] );
			vertex_pairs[EGTree[tn].union_borders[k]][cands[p]] = result[p];
			}
			}
			*/

			// IMPORTANT! after all border finished, degenerate graph,###用于简化Dijkstra算法距离计算
			// first, remove inward edges
			for (int k = 0; k < EGTree[tn].borders.size(); k++) {
				s = EGTree[tn].borders[k];
				tnodes.clear();
				tweight.clear();
				for (int p = 0; p < graph[s].adjnodes.size(); p++) {
					nid = graph[s].adjnodes[p];
					weight = graph[s].adjweight[p];
					// if adj node in same tree node
					// ????
					if (graph[nid].egtreepath.size() <= i || graph[nid].egtreepath[i] != tn) {
						// only leave those useful
						tnodes.push_back(nid);
						tweight.push_back(weight);

					}
				}
				// cut it
				graph[s].adjnodes = tnodes;
				graph[s].adjweight = tweight;
			}
			// second, add inter connected edges
			for (int k = 0; k < EGTree[tn].borders.size(); k++) {
				for (int p = 0; p < EGTree[tn].borders.size(); p++) {
					if (k == p) continue;
					s = EGTree[tn].borders[k];
					t = EGTree[tn].borders[p];
					graph[s].adjnodes.push_back(t);
					graph[s].adjweight.push_back(vertex_pairs[s][t]);
				}
			}
		}
	}
}

// --------------------------M-- no use------------
/*
// dump distance matrix into file
void hierarchy_shortest_path_save(){
FILE* fout = fopen( FILE_ONTREE_MIND, "wb" );
int* buf;
int count;
for ( int i = 0; i < EGTree.size(); i++ ){
// union borders
count = EGTree[i].union_borders.size();
fwrite( &count, sizeof(int), 1, fout );
buf = new int[count];
copy( EGTree[i].union_borders.begin(), EGTree[i].union_borders.end(), buf );
fwrite( buf, sizeof(int), count, fout );
delete[] buf;
// mind
count = EGTree[i].mind.size();
fwrite( &count, sizeof(int), 1, fout );
buf = new int[count];
copy( EGTree[i].mind.begin(), EGTree[i].mind.end(), buf );
fwrite( buf, sizeof(int), count, fout );
delete[] buf;
}
fclose(fout);
}

// --------------------------M-- no use------------
// load distance matrix from file
void hierarchy_shortest_path_load(){
FILE* fin = fopen( FILE_ONTREE_MIND, "rb" );
int* buf;
int count, pos = 0;
while( fread( &count, sizeof(int), 1, fin ) ){
// union borders
buf = new int[count];
fread( buf, sizeof(int), count, fin );
EGTree[pos].union_borders.clear();
for ( int i = 0; i < count; i++ ){
EGTree[pos].union_borders.push_back(buf[i]);
}
delete[] buf;
// mind
fread( &count, sizeof(int), 1, fin );
buf = new int[count];
fread( buf, sizeof(int), count, fin );
EGTree[pos].mind.clear();
for ( int i = 0; i < count; i++ ){
EGTree[pos].mind.push_back(buf[i]);
}
pos++;
delete[] buf;
}
fclose(fin);
}
*/

int mainFunction(int nOfNode, EdgeMapType EdgeMap) { // main function{
													 // init
	TIME_TICK_START
		init(nOfNode, EdgeMap);
	TIME_TICK_END
		TIME_TICK_PRINT("INIT")

		// gtree_build
		TIME_TICK_START
		build(EdgeMap);
	TIME_TICK_END
		TIME_TICK_PRINT("BUILD")

		// calculate distance matrix
		TIME_TICK_START
		hierarchy_shortest_path_calculation();
	TIME_TICK_END
		TIME_TICK_PRINT("MIND")

		// dump EGTree
		//-------------------------------M-- save in genData----
		//egtree_save();
		// dump distance matrix
		//hierarchy_shortest_path_save();

		return 0;
}