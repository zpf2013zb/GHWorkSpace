#include "StdAfx.h"
#include "EGTree.h"

// METIS setting options
void options_setting(){
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
void init_input(int nOfNode, EdgeMapType EdgeMap) {
	nLeafNode = 0;
	// process node information
	printf("PROCESSING NODE...");
	//fin = fopen(FILE_NODE,"r");
	for(int i=0; i<nOfNode; i++) {
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
	EdgeMapType::iterator iter=EdgeMap.begin();
	while (iter!=EdgeMap.end()) {
		edge* e=iter->second;
		noe ++;
		snid = e->Ni;
		enid = e->Nj;
		weight = e->dist;

		iweight = (int) (weight * WEIGHT_INFLATE_FACTOR);
		Nodes[snid].adjnodes.push_back( enid );
		Nodes[snid].adjweight.push_back( iweight );
		Nodes[enid].adjnodes.push_back( snid );
		Nodes[enid].adjweight.push_back( iweight );		
		iter++;
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
void data_transform_init( set<int> &nset ){
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
	unordered_map<int,int> nodemap;
	nodemap.clear();

	xadj[0] = 0;
	int i = 0;
	for ( set<int>::iterator it = nset.begin(); it != nset.end(); it++, i++ ){
		// init node map
		nodemap[*it] = i;

		int nid = *it;
		int fanout = Nodes[nid].adjnodes.size();
		for ( int j = 0; j < fanout; j++ ){
			int enid = Nodes[nid].adjnodes[j];
			// ensure edges within
			if ( nset.find( enid ) != nset.end() ){
				// xadj_accum used to record the start and end position of adjacent nodes for current visited nodes
				xadj_accum ++;

				adjncy[adjncy_pos] = enid;
				adjwgt[adjncy_pos] = Nodes[nid].adjweight[j];
				adjncy_pos ++;
			}
		}
		xadj[xadj_pos++] = xadj_accum;
	}

	// adjust nodes number started by 0  ###########这部分以及下面部分要修改，为何要做映射，为什么权重改为1？(可以继续使用)
	for ( int i = 0; i < adjncy_pos; i++ ){
		adjncy[i] = nodemap[adjncy[i]];
	}

	// adjwgt -> 1
	if (ADJWEIGHT_SET_TO_ALL_ONE){
		for ( int i = 0; i < adjncy_pos; i++ ){
			adjwgt[i] = 1;
		}
	}

	// nparts
	nparts = PARTITION_PART;

	// part
	part = new idx_t[nset.size()];
}

void init(int nOfNode, const EdgeMapType EdgeMap){
	init_input(nOfNode, EdgeMap);
	options_setting();
}

void finalize(){
    delete xadj;
    delete adjncy;
    delete adjwgt;
    delete part; 
}

// graph partition
// input: nset = a set of node id
// output: <node, node belong to partition id>
unordered_map<int,int> graph_partition( set<int> &nset ){
	unordered_map<int,int> result;

	// transform data to metis
	data_transform_init( nset );		

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
	if ( METISERRORINPUT == returnstate  ) {
		printf("Input Error!");
		return result;
	} else if(METISERRORMEMORY == returnstate) {
		printf("Memory Error!");
		return result;
	} else if(METISERROR == returnstate) {
		printf("Other Type Error!");
		return result;
	} else {

	}

	// push to result
	result.clear(); 
	int i = 0;
	// 又将结果变回来了，在上面采用nodemap可能是为了满足一定的条件，part中的结果是什么呢？子部分用什么来表示呢？
	for ( set<int>::iterator it = nset.begin(); it != nset.end(); it++, i++ ){
		result[*it] = part[i];
	}

	// finalize
	finalize();

	return result;
}

// egtree construction
void build(){
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
	for ( int i = 0; i < Nodes.size(); i++ ){
		rootstatus.nset.insert(i);
	}
	buildstack.push( rootstatus );  

	// start to build
	unordered_map<int,int> presult;
	set<int> childset[PARTITION_PART];


	while( buildstack.size() > 0 ){
		// pop top
		Status current = buildstack.top();
		buildstack.pop();

		// update egtreepath
		for ( set<int>::iterator it = current.nset.begin(); it != current.nset.end(); it++ ){
			Nodes[*it].egtreepath.push_back( current.tnid );
		}

		// check cardinality
		if ( current.nset.size() <= LEAF_CAP ){
			// build leaf node
			nLeafNode++;
			EGTree[current.tnid].isleaf = true;
			EGTree[current.tnid].leafnodes.clear();
			for ( set<int>::iterator it = current.nset.begin(); it != current.nset.end(); it++ ){
				EGTree[current.tnid].leafnodes.push_back( *it );
				//---------修改partID
				//partID[*it] = current.tnid;
			}
			continue;
		}

		// partition
//		printf("PARTITIONING...NID=%d...SIZE=%d...", current.tnid, (int)current.nset.size() );
		presult = graph_partition( current.nset );
//		printf("COMPLETE.\n");

		// construct child node set
		for ( int i = 0; i < PARTITION_PART; i++ ){
			childset[i].clear();
		}
		int slot;
		// put the nodes into corresponding sub-partition(slot)
		for ( set<int>::iterator it = current.nset.begin(); it != current.nset.end(); it++ ){
			slot = presult[*it];
			childset[slot].insert(*it);
		}

		// generate child tree nodes
		int childpos;
		for ( int i = 0; i < PARTITION_PART; i++ ){
			TreeNode tnode;
			tnode.isleaf = false;
			tnode.father = current.tnid;
			
			// insert to EGTree first
            EGTree.push_back(tnode);
			childpos = EGTree.size() - 1;
			EGTree[current.tnid].children.push_back( childpos );

			// calculate border nodes
			EGTree[childpos].borders.clear();
			for ( set<int>::iterator it = childset[i].begin(); it != childset[i].end(); it++ ){

				bool isborder = false;
				for ( int j = 0; j < Nodes[*it].adjnodes.size(); j++ ){
					if ( childset[i].find( Nodes[*it].adjnodes[j] ) == childset[i].end() ){
						isborder = true;
						break;
					}
				}
				if ( isborder ){
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
void egtree_save(){
	// FILE_GTREE
	FILE *fout = fopen( FILE_GTREE, "wb" );
	int *buf = new int[ Nodes.size() ];
	for ( int i = 0; i < EGTree.size(); i++ ){
		// borders
		int count_borders = EGTree[i].borders.size();
		fwrite( &count_borders, sizeof(int), 1, fout );
		copy( EGTree[i].borders.begin(), EGTree[i].borders.end(), buf );
		fwrite( buf, sizeof(int), count_borders, fout );
		// children
		int count_children = EGTree[i].children.size();
		fwrite( &count_children, sizeof(int), 1, fout );
		copy( EGTree[i].children.begin(), EGTree[i].children.end(), buf );
		fwrite( buf, sizeof(int), count_children, fout );
		// isleaf
		fwrite( &EGTree[i].isleaf, sizeof(bool), 1, fout );
		// leafnodes
		int count_leafnodes = EGTree[i].leafnodes.size();
		fwrite( &count_leafnodes, sizeof(int), 1, fout );
		copy( EGTree[i].leafnodes.begin(), EGTree[i].leafnodes.end(), buf );
		fwrite( buf, sizeof(int), count_leafnodes, fout );
		// father
		fwrite( &EGTree[i].father, sizeof(int), 1, fout );
	}
	fclose(fout);

	// FILE_NODES_GTREE_PATH
	fout = fopen( FILE_NODES_GTREE_PATH, "wb" );
	for ( int i = 0; i < Nodes.size(); i++ ){
		int count = Nodes[i].egtreepath.size();
		fwrite( &count, sizeof(int), 1, fout );
		copy( Nodes[i].egtreepath.begin(), Nodes[i].egtreepath.end(), buf );
		fwrite( buf, sizeof(int), count, fout );
	}
	fclose(fout);
	delete[] buf;
}

// load EGTree index from file
void egtree_load(){
	// FILE_GTREE
	FILE *fin = fopen( FILE_GTREE, "rb" );
	int *buf = new int[ Nodes.size() ];
	int count_borders, count_children, count_leafnodes;
	bool isleaf;
	int father;

	// clear EGTree
	EGTree.clear();

	while( fread( &count_borders, sizeof(int), 1, fin ) ){
		TreeNode tn;
		// borders
		tn.borders.clear();
		fread( buf, sizeof(int), count_borders, fin );
		for ( int i = 0; i < count_borders; i++ ){
			tn.borders.push_back(buf[i]);
		}
		// children
		fread( &count_children, sizeof(int), 1, fin );
		fread( buf, sizeof(int), count_children, fin );
		for ( int i = 0; i < count_children; i++ ){
			tn.children.push_back(buf[i]);
		}
		// isleaf
		fread( &isleaf, sizeof(bool), 1, fin );
		tn.isleaf = isleaf;
		// leafnodes
		fread( &count_leafnodes, sizeof(int), 1, fin );
		fread( buf, sizeof(int), count_leafnodes, fin );
		for ( int i = 0; i < count_leafnodes; i++ ){
			tn.leafnodes.push_back(buf[i]);
		}
		// father
		fread( &father, sizeof(int), 1, fin );
		tn.father = father;

		EGTree.push_back(tn);
	}
	fclose(fin);
	
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
}

// dijkstra search, used for single-source shortest path search WITHIN one EGTree leaf node!
// input: s = source node
//        cands = candidate node list
//        graph = search graph(this can be set to subgraph)
vector<int> dijkstra_candidate( int s, vector<int> &cands, vector<Node> &graph ){
	// init
	set<int> todo;
	todo.clear();
	todo.insert(cands.begin(), cands.end());
	
	unordered_map<int,int> result;
	result.clear();
	set<int> visited;
	visited.clear();
	unordered_map<int,int> q;
	q.clear();
	q[s] = 0;

	// start
	int min, minpos, adjnode, weight;
	while( ! todo.empty() && ! q.empty() ){
		min = -1;
		for ( unordered_map<int,int>::iterator it = q.begin(); it != q.end(); it ++ ){
			if ( min == -1 ){
				minpos = it -> first;
				min = it -> second;
			}
			else{
				if ( it -> second < min ){
					min = it -> second;
					minpos = it -> first;
				}
			}
		}

		// put min to result, add to visited
		result[minpos] = min;
		visited.insert( minpos );
		q.erase(minpos);

		if ( todo.find( minpos ) != todo.end() ){
			todo.erase( minpos );
		}

		// expand
		for ( int i = 0; i < graph[minpos].adjnodes.size(); i++ ){
			adjnode = graph[minpos].adjnodes[i];
			if ( visited.find( adjnode ) != visited.end() ){
				continue;
			}
			weight = graph[minpos].adjweight[i];

			if ( q.find(adjnode) != q.end() ){
				if ( min + weight < q[adjnode] ){
					q[adjnode] = min + weight;
				}
			}
			else{
				q[adjnode] = min + weight;
			}
		
		}
	}

	// output
	vector<int> output;
	for ( int i = 0; i < cands.size(); i++ ){
		output.push_back( result[cands[i]] );
	}

	// return
	return output;
}

// calculate the distance matrix, algorithm shown in section 5.2 of paper
void hierarchy_shortest_path_calculation(){
	// level traversal
	vector< vector<int> > treenodelevel;
	
	vector<int> current;
	current.clear();
	current.push_back(0);
	treenodelevel.push_back(current);
	// put all the nodes into treenodelevel according to their levels
	vector<int> mid;
	while( current.size() != 0 ){
		mid = current;
		current.clear();
		for ( int i = 0; i < mid.size(); i++ ){
			for ( int j = 0; j < EGTree[mid[i]].children.size(); j++ ){
				current.push_back( EGTree[mid[i]].children[j] );
			}
		}
		if ( current.size() == 0 ) break;
		treenodelevel.push_back( current );
	}
	
	// bottom up calculation
	// temp graph
	vector<Node> graph;
	graph = Nodes;
	vector<int> cands;
	vector<int> result;
	unordered_map<int, unordered_map<int,int> > vertex_pairs;

	// do dijkstra
	int s, t, tn, nid, cid, weight;
	vector<int> tnodes, tweight;
	set<int> nset;

	for ( int i = treenodelevel.size() - 1; i >= 0; i-- ){
		for ( int j = 0; j < treenodelevel[i].size(); j++ ){
			tn = treenodelevel[i][j];

			cands.clear();
			if ( EGTree[tn].isleaf ){
				// cands = leafnodes
				cands = EGTree[tn].leafnodes;
				// union borders = borders;
				EGTree[tn].union_borders = EGTree[tn].borders;
			}
			else{
				nset.clear();
				for ( int k = 0; k < EGTree[tn].children.size(); k++ ){
					cid = EGTree[tn].children[k];
					nset.insert( EGTree[cid].borders.begin(), EGTree[cid].borders.end() );
				}
				// union borders = cands;
				
				cands.clear();
				for ( set<int>::iterator it = nset.begin(); it != nset.end(); it ++ ){
					cands.push_back( *it );
				}
				EGTree[tn].union_borders = cands;
			}
				
			// start to do min dis
			vertex_pairs.clear();
				
			// for each border, do min dis
			int cc = 0;

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

			// IMPORTANT! after all border finished, degenerate graph,###用于简化Dijkstra算法距离计算
			// first, remove inward edges
			for ( int k = 0; k < EGTree[tn].borders.size(); k++ ){
				s = EGTree[tn].borders[k];
				tnodes.clear();
				tweight.clear();
				for ( int p = 0; p < graph[s].adjnodes.size(); p++ ){
					nid = graph[s].adjnodes[p];
					weight = graph[s].adjweight[p];
					// if adj node in same tree node
// ????
					if ( graph[nid].egtreepath.size() <= i || graph[nid].egtreepath[i] != tn ){
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
			for ( int k = 0; k < EGTree[tn].borders.size(); k++ ){
				for ( int p = 0; p < EGTree[tn].borders.size(); p++ ){
					if ( k == p ) continue;
					s = EGTree[tn].borders[k];
					t = EGTree[tn].borders[p];
					graph[s].adjnodes.push_back( t );
					graph[s].adjweight.push_back( vertex_pairs[s][t] );
				}
			}
		}
	}
}

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
//---------------extend function of egtree
void makePtFiles(FILE *ptFile,char* treefile)
{
    PtMaxUsing=true;
    BTree* bt=initialize(treefile);
    printf("making PtFiles\n");

    int RawAddr=0,key=0,size;	// changed
    EdgeMapType::iterator iter=EdgeMap.begin();
    while (iter!=EdgeMap.end())
    {
        edge* e=iter->second;
        if (e->pts.size()>0)  	// do not index empty groups
        {
            sort(e->pts.begin(),e->pts.end(),ComparInerNode());

            RawAddr=ftell(ptFile);	// set addr to correct amt.
            size=e->pts.size();
            fwrite(&(e->Ni),1,sizeof(int),ptFile);
            fwrite(&(e->Nj),1,sizeof(int),ptFile);
            fwrite(&(e->dist),1,sizeof(float),ptFile);
            fwrite(&(size),1,sizeof(int),ptFile);
            fwrite(&(e->pts[0]),e->pts.size(),sizeof(InerNode),ptFile);
            e->FirstRow=key;
            PtMaxKey=key+e->pts.size()-1;	// useful for our special ordering !

            //printf("(key,value)=(%d,%d)\n",key,RawAddr);
            addentry(bt,&top_level,i_capacity,1,key,&num_written_blocks,RawAddr);
            key+=sizeof(int)*3+sizeof(float);
            key+=e->pts.size()*sizeof(InerNode); //Modified by Qin Xu
        }
        else
            e->FirstRow=-1;		// also later used by AdjFile

        iter++;
    }
    finalize(bt);
    bt->UserField=num_D;
    delete bt;
    PtMaxUsing=false;
}

// Adj FlatFile Field:
//		Header:	size(int)
//		Entry:	Nk(int), eDist(float), PtGrpKey(int), PtSize(int)		changed
void makeAdjListFiles(FILE *alFile)
{
    printf("making alFiles, dependency on makePtFiles\n");

    int key=0,size,PtSize;
    fwrite(&NodeNum,1,sizeof(int),alFile);

    // slotted header info.
    int addr=sizeof(int)+sizeof(int)*NodeNum;
    for (int Ni=1; Ni<=NodeNum; Ni++)
    {
        fwrite(&addr,1,sizeof(int),alFile);
        addr+=sizeof(int)+AdjList[Ni].size()*(2*sizeof(int)+sizeof(float));
    }

    float distsum=0;
    for (int Ni=1; Ni<=NodeNum; Ni++)
    {
        size=AdjList[Ni].size();
        fwrite(&(size),1,sizeof(int),alFile);

        for (int k=0; k<AdjList[Ni].size(); k++)
        {
            int Nk=AdjList[Ni][k];	// Nk can be smaller or greater than Ni !!!
            edge* e=EdgeMap[getKey(Ni,Nk)];
            PtSize=e->pts.size();
            fwrite(&Nk,1,sizeof(int),alFile);
            fwrite(&(e->dist),1,sizeof(float),alFile);
            fwrite(&(e->FirstRow),1,sizeof(int),alFile); // use FirstRow for other purpose ...
            //printf("(Ni,Nj,dataAddr)=(%d,%d,%d)\n",Ni,Nk,e->FirstRow);

            distsum+=e->dist;
        }
        key=Ni;
    }
    distsum=distsum/2;
    printf("total edge dist: %f\n",distsum);
    printf("total keywords num:%d\n",num_K);
}

void BuildBinaryStorage(const char* fileprefix)
{
    BlkLen=getBlockLength();
    char tmpFileName[255];

    FILE *ptFile,*edgeFile;
    sprintf(tmpFileName,"%s.p_d",fileprefix);
    remove(tmpFileName); // remove existing file
    ptFile=fopen(tmpFileName,"w+");
    sprintf(tmpFileName,"%s.p_bt",fileprefix);
    remove(tmpFileName); // remove existing file
    makePtFiles(ptFile,tmpFileName);

    sprintf(tmpFileName,"%s.al_d",fileprefix);
    remove(tmpFileName); // remove existing file
    edgeFile=fopen(tmpFileName,"w+");
    makeAdjListFiles(edgeFile);

    fclose(ptFile);
    fclose(edgeFile);
}

int mainFunction(int nOfNode, const EdgeMapType EdgeMap){ // main function{
	// init
	TIME_TICK_START
	init(int nOfNode, const EdgeMapType EdgeMap);
	TIME_TICK_END
	TIME_TICK_PRINT("INIT")

	// gtree_build
	TIME_TICK_START
	build();
	TIME_TICK_END
	TIME_TICK_PRINT("BUILD")

	// dump EGTree
	gtree_save();
	
	// calculate distance matrix
	TIME_TICK_START
	hierarchy_shortest_path_calculation();
	TIME_TICK_END
	TIME_TICK_PRINT("MIND")

	// dump distance matrix
	hierarchy_shortest_path_save();

	return 0;
}