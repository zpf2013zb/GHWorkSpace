#ifndef _BLKFILE_H_
#define _BLKFILE_H_
// handle
#define BFHEAD_LENGTH (sizeof(int)*2)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef char Block[]; //Bolck identical to char[]
#define TRUE 1
#define FALSE 0

class BlockFile
{
    FILE* fp;
    char* filename;
    int blocklength;
    int act_block; 	        // SEEK_CUR
    int number; // the number of block
    bool new_flag;
	// put the content of bytes length num to file 
    void put_bytes(const char* bytes,int num) 
    {
        fwrite(bytes,num,1,fp);
    }
	// read the content to bytes length num from file 
    void get_bytes(char* bytes,int num)
    {
        fread(bytes,num,1,fp);
    }
	// change the num to char* and wirte it to file from &int with lenght sizeof(int)
    void fwrite_number(int num);
	// read a int from file
    int fread_number();
	// set the file pointer to start postion of block bnum
    void seek_block(int bnum)
    {
        fseek(fp,(bnum-act_block)*blocklength,SEEK_CUR);
    }

public:
	// construct a blockFile with name and blocklength
	// exist then read blocklength and number else
	// creat a new file and initial the 0 block
    BlockFile(char* name, int b_length);
    ~BlockFile();
	// read the content from the head,[blocklength-BFHEAD_LENGTH]
    void read_header(char * header);
	// set the content to the head,[blocklength-BFHEAD_LENGTH]
    void set_header(char* header);
	// read the content of position i to block b, noteful (pos++)
    bool read_block(Block b,int i);
	// write the content of block b to position i
    bool write_block(const Block b,int i);
	// append block b to the end of fp 
    int append_block(const Block b);
    bool file_new()
    {
        return new_flag;
    }
    int get_blocklength()
    {
        return blocklength;
    }
    int get_num_of_blocks()
    {
        return number;
    }
};

class CachedBlockFile : public BlockFile
{
    enum uses {free,used};
    int ptr; // record the loop pointer in cachesize
    int cachesize; // record the size of cache
    int *cache_cont; // record the block id stored in this cache[i]
    uses *fuf_cont; // record if this cache[i] is used or not
    char **cache; // record the content of block in cache[i]
	// return a free space, otherwise write block ptr to file 
	// and return it as free space 
    int next();
	// if block index in cache, in then return the i where cache_cont[i] == index
    int in_cache(int index); 

public:
	// initial *cache_cont; *fuf_cont; **cache;
    CachedBlockFile(char* name,int blength, int csize);
    ~CachedBlockFile();
	// rewrite the read&writh_block of BlockFile
	// first read from cache, if not in then load it to cache;
	// otherwise, read from blockfile
    bool read_block(Block b,int i);
    bool write_block(const Block b,int i);
	// write cache content to blockFile
    void flush();
};

#endif //__BLKFILE
