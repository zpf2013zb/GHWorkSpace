#ifndef __BTREE
#define __BTREE
// handle
#include "blk_file.h"

class BTDirNode;

class BTree
{
public:
    friend class BTDirNode;
    int UserField;
    int root;						// block # of root node
    BTDirNode  *root_ptr;           // root-node
    int num_of_inodes;	        	// # of stored directory pages
    CachedBlockFile *file;	  		// storage manager for harddisc blocks
    void load_root();            	// loads root_node into memory
    char *header;
protected:
    char *user_header;
	// btree header, num_of_inodes, root, UserField ,user_header
    void read_header(char *buffer);      // reads Rtree header
    void write_header(char *buffer);     // writes Rtree header
public:
	// construct a new btree
    BTree(char *fname,int _b_length, int cache_size);
	// load a new btree
    BTree(char *fname, int cache_size);
    virtual ~BTree();
};

class DirEntry
{
    friend class BTDirNode ;
    friend class BTree ;

public:
    BTree  *my_tree;        // pointer to my B-tree
    BTDirNode  *son_ptr;       // pointer to son if in main mem.
    int key;              // key value
    int son;                // block # of son
	// load the son block if son_ptr is null
    BTDirNode * get_son();  	// returns the son, loads son if necessary
	// read and write key&son from or to buffer
    void read_from_buffer(char *buffer);// reads data from buffer
    void write_to_buffer(char *buffer); // writes data to buffer
	// key&son
    static const int EntrySize=sizeof(int)+sizeof(int);
    // returns amount of needed buffer space
	// copy data from _d to current DirEntry, except my_tree
    virtual DirEntry  & operator = (DirEntry  &_d);
	// set my_tree and son_ptr
    DirEntry(BTree *bt = NULL);
    virtual ~DirEntry();
};

class BTDirNode
{
    friend class BTree ;
public:
    BTree  *my_tree;               		// pointer to B-tree
    int capacity;                       // max. # of entries
    int num_entries;                    // # of used entries
    bool dirty;                         // TRUE, if node has to be written
    int block;                          // disc block
    char level;                         // level of the node in the tree

    DirEntry  *entries;            			// array of entries
    // read(write) level, num_entries and entries from(to) buffer
	void read_from_buffer(char *buffer);	// reads data from buffer
    void write_to_buffer(char *buffer); 	// writes data to buffer
	// construct a new BTDirNode 
    BTDirNode(BTree  *rt);
	// load a BTDirNode from _block
    BTDirNode(BTree  *rt, int _block);
    virtual ~BTDirNode();
};



#endif // __BTREE
