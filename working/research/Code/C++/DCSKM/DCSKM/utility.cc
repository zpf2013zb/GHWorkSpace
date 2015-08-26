#include "utility.h"

int START_TIME;
double ZipfMaxVal;

void InitClock()
{
    START_TIME=clock();
    srand(time(NULL));
    //srand48(time(NULL));
}
void PrintElapsed()
{
    printf("Elapsed: %d msec\n",clock()-START_TIME);//May not denote timeval??Qin Xu
}

void CheckFile(FILE* fp,const char* filename)
{
    if (fp==NULL)
    {
        printf("Invalid file '%s'\n",filename);
        exit(0);
    }
}

void AddAll(FastList<int> &S,FastList<int> &C)
{
    for (int i=0; i<C.size(); i++) S.push_back(C[i]);
}

// remove elements of C from S 
void RemoveAll(FastList<int> &S,FastList<int> &C)
{
    // new version experienced 65 times faster than	old version

    /* [old RemoveAll]
    for (int i=0;i<C.size();i++)
    	S.erase(remove(S.begin(),S.end(),C[i]),S.end());*/

    FastList<int> Rst;

    sort(S.begin(),S.end());	// sorted in ascending order
    sort(C.begin(),C.end());

    int i=0,j=0;
    while (i<S.size()&&j<C.size())
    {
        if (S[i]==C[j])
        {
            i++;
            j++;
        }
        else if (S[i]<C[j])  	// move i
        {
            while (i<S.size()&&S[i]<C[j])
            {
                Rst.push_back(S[i]);
                i++;
            }
        }
        else if (S[i]>C[j])  	// move j
        {
            while (j<C.size()&&S[i]>C[j])
                j++;
        }
    }

    while (i<S.size())  	// all elements of C should have been examined
    {
        Rst.push_back(S[i]);
        i++;
    }

    S.assign(Rst.begin(),Rst.end());
    Rst.clear();
}

void TrimSpace(char* str)
{
    if (str==NULL) return;

    char space[]= {'\t','\n','\f','\r',' '};
    int pos=0;
    for (int i=0; i<strlen(str); i++)
    {
        bool found=false;
        for (int j=0; j<5; j++)
            if (str[i]==space[j]) found=true;

        if (!found)
        {
            str[pos]=str[i];
            pos++;
        }
    }
    str[pos]='\0';
}

//void AddConfigFromFile(ConfigType &cr,const char* filename)
//{
//    const int LINE_LEN=1024;
//    char line[LINE_LEN],key[LINE_LEN],value[LINE_LEN];
//
//    ifstream br(filename);
//    if (! br.is_open())
//    {
//        printf("Error opening file \"%s\"",filename);
//        exit (1);
//    }
//
//    while (br.getline(line,LINE_LEN))
//    {
//        if (strstr(line,"//")!=NULL) continue; // remove comments
//        char* chPos=strchr(line,'=');
//        if (chPos!=NULL)
//        {
//            int pos=((int)(chPos-line))/sizeof(char);
//            int keyLen=pos;
//            int valueLen=strlen(line)-1-keyLen;
//            memcpy(key,&line[0],keyLen);
//            key[keyLen]='\0';
//            memcpy(value,&line[pos+1],valueLen);
//            value[valueLen]='\0';
//            TrimSpace(key);
//            TrimSpace(value);
//            cr[key]=value;
//        }
//    }
//    br.close();
//}
//
//void AddConfigFromCmdLine(ConfigType &cr,int argc,char** argv)
//{
//    int i=0;
//    while (i<argc)
//    {
//        while ((i<argc)&&(argv[i][0]!='-')) i++;	// shortcut condition
//        if (i+1<argc)
//        {
//            char* key=&(argv[i][1]);
//            char* value=argv[i+1];
//            TrimSpace(key);
//            TrimSpace(value);
//            cr[key]=value;
//            i+=2;
//        }
//        else
//            return;
//    }
//}
//
//void ListConfig(ConfigType &cr)
//{
//    ConfigType::iterator p=cr.begin();
//    while (p!=cr.end())
//    {
//        printf("%s=%s\n",p->first.c_str(),p->second.c_str());
//        p++;
//    }
//}
//
//float getConfigFloat(const char* key,ConfigType &cr,bool required,float _default)
//{
//    float value=_default;
//    if (cr.count(key))
//        value=atof(cr[key].c_str());
//    else
//    {
//        if (required)
//        {
//            printf("Config key \"%s\" not found\n",key);
//            exit(1);
//        }
//    }
//    return value;
//}
//
//int getConfigInt(const char* key,ConfigType &cr,bool required,int _default)
//{
//    int value=_default;
//    if (cr.count(key))
//        value=atoi(cr[key].c_str());
//    else
//    {
//        if (required)
//        {
//            printf("Config key \"%s\" not found\n",key);
//            exit(1);
//        }
//    }
//    return value;
//}
//
//const char* getConfigStr(const char* key,ConfigType &cr,bool required,const char* _default)
//{
//    const char* value=_default;
//    if (cr.count(key))
//        value=cr[key].c_str();
//    else
//    {
//        if (required)
//        {
//            printf("Config key \"%s\" not found\n",key);
//            exit(1);
//        }
//    }
//    return value;
//}
//
//int Cardinality(BitStore* bs)
//{
//    int Dcnt=0;
//    for (int i=0; i<bs->size(); i++)
//        if ((* bs)[i]) Dcnt++;
//    return Dcnt;
//}

void printPattern(BitStore* bs) 
{
    printf("{");
    for (int i=0; i<bs->size(); i++)
        if ((* bs)[i]) printf("%d ",i);
    printf("}");
}

int getSetIndex(vector<string> &vec,char* str)
{
    for (int i=0; i<vec.size(); i++)
        if (vec[i].compare(str)==0) return i;

    // otherwise, insert it into vec
    vec.push_back(str);
    return (vec.size()-1);
}

//returns random numbers following a gaussian (normal) distribution
double gaussian(double mean, double sigma)
{
    double v1,v2,s,x;
    do
    {
        //v1 = drand48();
		v1 = (float)rand() / RAND_MAX;
        if (rand()%2) v1 = -v1;
        //v2 = drand48();
		v2 = (float)rand() / RAND_MAX;
        if (rand()%2) v2 = -v2;
        s = v1*v1 + v2*v2;
    }
    while (s >= 1.);
    x = v1 * sqrt ( -2. * log (s) / s);
    //  x is normally distributed with mean 0 and sigma 1.
    x = x * sigma + mean;
    return x;
}

void InitZipfMaxVal(int maxnum,double theta)
{
    ZipfMaxVal = 0.0;
    for(int i=1; i<=maxnum; i++)
        ZipfMaxVal += 1.0/pow( (double)i, theta);
}

int zipf(double theta)
{
    //double r= drand48()*ZipfMaxVal, sum= 1.0;
	double r = ((float)rand() / RAND_MAX) *ZipfMaxVal, sum = 1.0;
    int i=1;
    while( sum<r)
    {
        i++;
        sum += 1.0/pow( (double)i, theta);
    }
    return (i-1);
}

long poisson(long lambda)
{
    long n = 0;
    double c = exp((double)(-lambda));
    double p = 1.0;
    while (p>=c)
    {
        p = p * ((float)rand() / RAND_MAX);
        n++;
    }
    //return (n-1); //this is in original code
    return n; //this is changed by me
}

float AvgAbsMeanDev(vector<float> &S)
{
    int S_len=S.size();
    if (S_len<=1) return 0;

    float sum=0,diff=0,mean=0;
    for (int p=0; p<S_len; p++) sum+=S[p];
    mean=sum/S_len;
    for (int p=0; p<S_len; p++) diff+=ValAbs(mean-S[p]);

    return ((sum==0)?(1):(diff/sum));
}

float variance(vector<float> &S)
{
    int S_len=S.size();
    if (S_len<=1) return 0;

    float sum=0,sqsum=0;
    for (int p=0; p<S_len; p++)
    {
        sum+=S[p];
        sqsum+=S[p]*S[p];
    }

    float sse=sqsum-(sum*sum/S_len);
    if (sse<0) sse=0;// due to precision
    return (sse/(S_len-1));
}

//-----------
// LRU buffer
//-----------

int SEQ_PAGE_ACCESS=0;
int RAN_PAGE_ACCESS=0;
int INDEX_PAGE_ACCESS=0;
int DATA_PAGE_ACCESS=0;
int CACHE_ACCESSED=0;

int blocklength=DEFAULT_BLOCKLENGTH;
int cachesize;	 	// max. number of cache blocks

// old: int *cache_cont; store transBlockID
int *cache_block,*cache_user; 	// block Id, user (file ID) separate !

int *lastTime;		// >=0 means used
int nextTime;
char **cache;		// cache cont

//int LastBlockIdStored;	// for determining SEQ or RAN
int LastBlockId,LastUserId;

int getBlockLength()
{
    return blocklength;
}

void InitCache(int csize)
{
    assert(blocklength>0);
    assert(csize>0);
    LastBlockId=-2;
    LastUserId=-2;
    nextTime=0;
    cachesize=csize;
    cache_block = new int[cachesize];
    cache_user = new int[cachesize];

    cache = new char*[cachesize];
    lastTime=new int[cachesize];
    for (int i=0; i<cachesize; i++)
    {
        lastTime[i]=-1;
        cache_block[i]=-1;
        cache_user[i]=-1;
        cache[i] = new char[blocklength];
    }
}

void RefreshStat()
{
    SEQ_PAGE_ACCESS=0;
    RAN_PAGE_ACCESS=0;
    INDEX_PAGE_ACCESS=0;
    DATA_PAGE_ACCESS=0;
    CACHE_ACCESSED=0;
}

void RefreshCache()  	// call before query execution
{
    LastBlockId=-2;
    LastUserId=-2;
    nextTime=0;
    for (int i=0; i<cachesize; i++)
    {
        lastTime[i]=-1;
        cache_block[i]=-1;
        cache_user[i]=-1;
    }
}

void DestroyCache()
{
    // no need to flush the cache (since the query alg. are read-only )
    delete[] cache_block;
    delete[] cache_user;
    delete[] lastTime;
    for (int i=0; i<cachesize; i++) delete[] cache[i];
    delete[] cache;
}

// light-weight get function, return the reference of cache[i]
char* getCacheBlockRef(int UserId,int BlockId)
{
    CACHE_ACCESSED++;
    for (int i=0; i<cachesize; i++)
        if ((cache_block[i]==BlockId)&&(cache_user[i]==UserId)&&
                (lastTime[i]>=0))
        {
            lastTime[i]=nextTime++;
            return cache[i];
        }
    return NULL;	// NULL if non-existent
}

bool getCacheBlock(char* buffer,int UserId,int BlockId)
{
    CACHE_ACCESSED++;
    for (int i=0; i<cachesize; i++)
        if ((cache_block[i]==BlockId)&&(cache_user[i]==UserId)&&
                (lastTime[i]>=0))
        {
            memcpy(buffer,cache[i],blocklength);
            lastTime[i]=nextTime++;
            return true;
        }
    return false;
}

// user's responsibility
// the place for counting number of page accesses
void storeCacheBlock(char* buffer,int UserId,int BlockId,AccessMode mode)
{
    int index=-1;
    for (int i=0; i<cachesize; i++)	// search for empty block
        if (lastTime[i]<0)
        {
            index=i;
            break;
        }

    if (index<0)
    {
        index=0;	// full, evict LRU block
        for (int i=0; i<cachesize; i++)
            if (lastTime[i]<lastTime[index]) index=i;
    }

    memcpy(cache[index],buffer,blocklength);
    cache_block[index]=BlockId;
    cache_user[index]=UserId;
    lastTime[index]=nextTime++;

    if (mode==isData)
        DATA_PAGE_ACCESS++;
    else if (mode==isIndex)
        INDEX_PAGE_ACCESS++;

    if ((LastUserId==UserId)&&(abs(LastBlockId-BlockId)<=1))
        SEQ_PAGE_ACCESS++;
    else
        RAN_PAGE_ACCESS++;
    LastBlockId=BlockId;
    LastUserId=UserId;
}

int64_t printPageAccess()
{
    //printf("Cache requests:\t%d\n",CACHE_ACCESSED);
    //printf("Index page accesses:\t%d\n",INDEX_PAGE_ACCESS);
    //printf("Data page accesses:\t%d\n",DATA_PAGE_ACCESS);
    //printf("Seq.:\t%d\n",SEQ_PAGE_ACCESS);
    //printf("Ran.:\t%d\n",RAN_PAGE_ACCESS);
    //printf("Total page accesses:\t%d\n",INDEX_PAGE_ACCESS+DATA_PAGE_ACCESS);
    return int64_t(DATA_PAGE_ACCESS);
}


//-----------
// Freq Cache
//-----------

void InitFreqCache(FreqCache& fc)
{
    fc.buffer=new char[blocklength];	// assume BlkLen inited
    fc.UserId=-1;
    fc.BlockId=-1;
}

// assume buf of size BlkLen
void storeFreqCache(FreqCache& fc,char* buf,int uid,int bid)
{
    memcpy(fc.buffer,buf,blocklength);
    fc.UserId=uid;
    fc.BlockId=bid;
}

bool inFreqCache(FreqCache& fc,int uid,int bid)
{
    return (fc.UserId==uid&&fc.BlockId==bid);
}

void DestroyFreqCache(FreqCache& fc)
{
    delete[] fc.buffer;
}

