#include "btree.h"

#include "utility.h"
#include "netshare.h"

float alpha;
int PtNum;
int PtFileSize,AdjFileSize;
FILE *PtFile,*AdjFile;
BTree *PtTree;
int i_capacity;
char *gBuffer;	// temp, global buffer for atomic use
FreqCache FC_A,FC_P; // FC for ensureCache exclusive use !
int BlkLen;

#define PtTid 	1000
#define PtFid 	2000
#define AdjTid 	3000
#define AdjFid 	4000
#define Ref(obj) ((void*)&obj)

// no use
struct NodeStruct
{
    char level;
    int num_entries;
    int* key;
    int* son;
};
// no use
NodeStruct* createNodeStruct()
{
    NodeStruct* node=new NodeStruct();
    node->level=-1;
    node->num_entries=-1;
    node->key=new int[i_capacity];
    node->son=new int[i_capacity];
    return node;
}
// no use
void ReadIndexBlock(BTree* bt,int block,NodeStruct* node)
{
    int UserId=(bt==PtTree)?PtTid:AdjTid;
    if (!getCacheBlock(gBuffer,UserId,block))
    {
        CachedBlockFile* cf=bt->file;
        cf->read_block(gBuffer,block);
        storeCacheBlock(gBuffer,UserId,block,isIndex);
    }

    int j=0;	// read node header
    memcpy(&(node->level), &gBuffer[j], sizeof(char));
    j += sizeof(char);
    memcpy(&(node->num_entries), &gBuffer[j], sizeof(int));
    j += sizeof(int);

    // read node content
    for (int i=0; i<node->num_entries; i++)
    {
        memcpy(&(node->key[i]),&gBuffer[j],sizeof(int));
        memcpy(&(node->son[i]),&gBuffer[j+sizeof(int)],sizeof(int));
        j += sizeof(int)+sizeof(int);
    }
}
// no use
int inline pointQuery(BTree* bt,int key,int& TreeKey)
{
    static NodeStruct* node=createNodeStruct();
    TreeKey=-2;
    if (key<0) return -3;	// assume non-negative keys

    ReadIndexBlock(bt,bt->root,node);
    while (true)
    {
        int curLevel=(int)(node->level);
        int Son=-1;

        if (curLevel>1)
        {
            for (int i=(node->num_entries-1); i>=0; i--)
                if (key<=node->key[i]) Son=i;
            if (Son<0) return -4;
            ReadIndexBlock(bt,node->son[Son],node);
        }
        else      // curLevel is 1
        {
            //printf("%d %d %d\n",key,node->key[0],node->key[node->num_entries-1]);
            for (int i=(node->num_entries-1); i>=0; i--)
            {
                if (key>=node->key[i])   // use a different test than above
                {
                    TreeKey=node->key[i];
                    return node->son[i];
                }
            }
            return -5;
        }
    }
    return -6;
}
// read the blockID from FC-cache-myfile
char* getFlatBlock(FILE* myfile,int BlockId)
{
    int UserId,FileSize;
    FreqCache* FC;

    if (myfile==PtFile)
    {
        FC=&FC_P;
        UserId=PtFid;
        FileSize=PtFileSize;
    }
    else
    {
        FC=&FC_A;
        UserId=AdjFid;
        FileSize=AdjFileSize;
    }
    if (inFreqCache(*FC,UserId,BlockId)) return (FC->buffer);

    if (!getCacheBlock(FC->buffer,UserId,BlockId))
    {
        int readsize=BlkLen<(FileSize-BlockId*BlkLen)?BlkLen:(FileSize-BlockId*BlkLen);//Modified By Qin Xu
        fseek(myfile,BlockId*BlkLen,SEEK_SET);
        fread(FC->buffer,readsize,sizeof(char),myfile);
        storeCacheBlock(FC->buffer,UserId,BlockId,isData);
    }
    FC->UserId=UserId;
    FC->BlockId=BlockId;
    return (FC->buffer);
}

//Modified by Qin Xu
int PTGRP_HEADSIZE=3*sizeof(int)+sizeof(float);
int PTGRP_ITEMSIZE=sizeof(InerNode);
int ADJGRP_HEADSIZE=sizeof(int);
int ADJGRP_ITEMSIZE=2*sizeof(int)+sizeof(float);	// changed

enum FixedF {SIZE_A,NI_P,NJ_P,DIST_P,SIZE_P};
enum VarE {ADJNODE_A,DIST_A,PTKEY_A,PT_P,PT_DIST,PT_ATTRIBUTE,PT_NOK,PT_VCT};
// get the variant from the adjFile or PtFile, put the pos item from addr BaseAddr to buffer
void getVarE(VarE type,void* buf,int BaseAddr,int pos)
{
    // note default values !
    FILE* f=AdjFile;
    int size=sizeof(int),addr=-1;
    int VarBase=BaseAddr+ADJGRP_HEADSIZE+pos*ADJGRP_ITEMSIZE;

    // for VarE in AdjFile
    if (type==ADJNODE_A) addr=VarBase;
    if (type==DIST_A)
    {
        addr=VarBase+sizeof(int);
        size=sizeof(float);
    }
    if (type==PTKEY_A) addr=VarBase+sizeof(int)+sizeof(float);

    // for VarE in PtFile is distance
    if (type==PT_P)
    {
        addr=BaseAddr+PTGRP_HEADSIZE+pos*PTGRP_ITEMSIZE;
        f=PtFile;
        //size=sizeof(float);
    }
	 if (type==PT_DIST)
    {
        addr=BaseAddr+PTGRP_HEADSIZE+pos*PTGRP_ITEMSIZE+sizeof(int);
        f=PtFile;
        size=sizeof(float);
    }

    //Add by Qin xu
    //for Keywords of  POI
	
	 if (type==PT_ATTRIBUTE)
    {
        addr=BaseAddr+PTGRP_HEADSIZE+pos*PTGRP_ITEMSIZE+sizeof(int)+sizeof(float);
        f=PtFile;
        size=ATTRIBUTE_DIMENSION*sizeof(float);
    }
    if (type==PT_NOK)
    {
        addr=BaseAddr+PTGRP_HEADSIZE+pos*PTGRP_ITEMSIZE+sizeof(int)+sizeof(float)+ATTRIBUTE_DIMENSION*sizeof(float);
        f=PtFile;
        size=sizeof(int);
    }
	 if (type==PT_VCT) //????是不是只存放了一个指针
    {
        addr=BaseAddr+PTGRP_HEADSIZE+pos*PTGRP_ITEMSIZE+sizeof(int)+sizeof(float)+ATTRIBUTE_DIMENSION*sizeof(float);
		int addrNOK = addr;
		addr= sizeof(int);
        f=PtFile;
		int NOK;
		char* BlockAddr=getFlatBlock(f,addrNOK/BlkLen);
		memcpy(&NOK,BlockAddr+(addrNOK%BlkLen),sizeof(int));
        size=sizeof(int)*NOK;
    }

    char* BlockAddr=getFlatBlock(f,addr/BlkLen);
    memcpy(buf,BlockAddr+(addr%BlkLen),size);
}
// put the head information from addr BaseAddr to buf
void getFixedF(FixedF type,void* buf,int BaseAddr)
{
    // note default values !
    FILE* f=PtFile;
    int size=sizeof(int),addr=-1;

    // for FixedF in PtFile
    if (type==NI_P) addr=BaseAddr;
    if (type==NJ_P) addr=BaseAddr+sizeof(int);
    if (type==DIST_P)
    {
        addr=BaseAddr+2*sizeof(int);
        size=sizeof(float);
    }
    if (type==SIZE_P)
    {
        addr=BaseAddr+2*sizeof(int)+sizeof(float);
    }

    // for FixedF in AdjFile
    if (type==SIZE_A)
    {
        addr=BaseAddr;
        f=AdjFile;
    }

    char* BlockAddr=getFlatBlock(f,addr/BlkLen);
    memcpy(buf,BlockAddr+(addr%BlkLen),size);
}
// return the NodeID address in AdjList, note the NodeID from 1
int getAdjListGrpAddr(int NodeID)  	// using AdjFile
{
	//----------NodeID start from 1-----------
    //int addr=sizeof(int)+NodeID*sizeof(int),GrpAddr;
	int addr=NodeID*sizeof(int),GrpAddr;
    char* BlockAddr=getFlatBlock(AdjFile,addr/BlkLen);
    memcpy(Ref(GrpAddr),BlockAddr+(addr%BlkLen),sizeof(int));
    return GrpAddr;
}

int getFileSize(FILE* f)  	// side effect, setting f to begin
{
    fseek(f,0,SEEK_END);
    int filesize=ftell(f);
    fseek(f,0,SEEK_SET);
    return filesize;
}
// open the file and index
void OpenDiskComm(const char* fileprefix,int _cachesize)
{
    char tmpFileName[255];

    BlkLen=getBlockLength();
    int header_size = sizeof(char) + sizeof(int);
    gBuffer=new char[BlkLen];
    i_capacity=(BlkLen-header_size)/(sizeof(int)+sizeof(int));
    printf("Blocklength=%d, i_cap=%d\n",BlkLen,i_capacity);

    InitFreqCache(FC_A);
    InitFreqCache(FC_P);
    InitCache(_cachesize);
    sprintf(tmpFileName,"%s.p_d",fileprefix);
    PtFile=fopen(tmpFileName,"r");
    CheckFile(PtFile,tmpFileName);
    PtFileSize=getFileSize(PtFile);
    sprintf(tmpFileName,"%s.p_bt",fileprefix);
	// load btree information from file and construct 128 cache block
    PtTree=new BTree(tmpFileName, 128); // number of cached file entries: 128
    PtNum=PtTree->UserField;

    sprintf(tmpFileName,"%s.al_d",fileprefix);
    AdjFile=fopen(tmpFileName,"r");
    CheckFile(AdjFile,tmpFileName);
    fread(Ref(NodeNum),1,sizeof(int),AdjFile);	//	NodeNum=AdjTree->UserField;
    AdjFileSize=getFileSize(AdjFile);
    printf("PtFileSize: %d, AdjFileSize: %d, NodeNum: %d, PtNum: %d\n",PtFileSize,AdjFileSize,NodeNum,PtNum);
}

void CloseDiskComm()
{
    fclose(PtFile);
    fclose(AdjFile);
    delete PtTree;
    delete[] gBuffer;

    //printPageAccess();
    DestroyCache();
    DestroyFreqCache(FC_A);
    DestroyFreqCache(FC_P);
}
// no use
void printSubTree(BTree* bt,NodeStruct* node,int indent)
{
    char space[indent+1];
    for (int i=0; i<indent; i++) space[i]=' ';
    space[indent]='\0';

    int curLevel=(int)(node->level);
    printf("%sLevel: %d\n",space,curLevel);

    NodeStruct* cNode=createNodeStruct();
    for(int i=0; i<node->num_entries; i++)
    {
        printf("%s(%d,%d)\n",space,node->key[i],node->son[i]);
        if (curLevel>1)
        {
            ReadIndexBlock(bt,node->son[i],cNode);
            printSubTree(bt,cNode,indent+1);
        }
        else if (curLevel==1)
        {
            // fetch data page
        }
    }
    printf("\n");
}
//no use
void printTree(BTree* bt)
{
    NodeStruct* root_node=createNodeStruct();
    ReadIndexBlock(bt,bt->root,root_node);
    printSubTree(bt,root_node,0);
}

// functions for visualizing clusters

FastArray<float> xcrd,ycrd;
FILE** views;	//=fopen(visualf,"w");
int Cnum;
// no use
void openVisualFiles(const char* nodename,const char* outprefix,int _Cnum)
{
    char nodef[255],visualf[255];
    sprintf(nodef,"data/%s",nodename);
    FILE* cnode=fopen(nodef,"r");

    int id;//Ni,Nj;Modified By Qin Xu
    float x,y;//,dist,;

    // read node info with format: NodeId xcrd ycrd
    while (!feof(cnode))  	// assume NodeID in ascending order from 0
    {
        fscanf(cnode,"%d %f %f\n", &id, &x, &y);
        xcrd.push_back(x);
        ycrd.push_back(y);
    }
    fclose(cnode);

    Cnum=_Cnum;
    views=new FILE*[Cnum+1];		// last file for outlier
    for (int i=0; i<Cnum+1; i++)  	// including file for outliers
    {
        if (i==Cnum)
            sprintf(visualf,"visual/%s.odat",outprefix);
        else
            sprintf(visualf,"visual/%s.pdat%d",outprefix,i);
        remove(visualf);
        views[i]=fopen(visualf,"w");
    }
}
// no use
void writeVisualRecord(int Ni,int Nj,float vP,float eDist,int Cid)
{
    float x1,y1,x2,y2,x,y;

    x1=xcrd[Ni];
    y1=ycrd[Ni];
    x2=xcrd[Nj];
    y2=ycrd[Nj];

    x=x1+(vP/eDist)*(x2-x1);
    y=y1+(vP/eDist)*(y2-y1);

    if (Cid<0)
        fprintf(views[Cnum],"%f %f\n",x,y);
    else
        fprintf(views[Cid],"%f %f\n",x,y);
}
// no use
void closeVisualFiles()  	// close files
{
    xcrd.clear();
    ycrd.clear();
    for (int i=0; i<Cnum+1; i++) fclose(views[i]);	// including file for outliers
    delete[] views;
}


