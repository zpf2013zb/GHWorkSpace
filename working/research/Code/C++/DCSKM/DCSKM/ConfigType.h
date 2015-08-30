//
//  ConfigType.h
//  generateQueries
//
//  Created by 秦旭 on 14-6-28.
//  Copyright (c) 2014年 ZJU. All rights reserved.
//
// handle
#ifndef __RKSK__ConfigType__
#define __RKSK__ConfigType__

#include <iostream>
#include <map>
#include <string>

//-------------------------M--test for different platform--------------
// for linux
/*
#define PRO_HOME_DIR std::string("/Users/pengfei/Desktop/rksk_mac/")
#define CONFIG_PATH  PRO_HOME_DIR
#define DATA_PATH PRO_HOME_DIR+std::string("data/")
#define QUERY_FILE_PATH PRO_HOME_DIR+std::string("queryfiles/")
#define QUERY_RESULT_FILE_PATH PRO_HOME_DIR+std::string("queryresults/")
*/
std::string seperator = "\\";
// for windows
#define PRO_HOME_DIR std::string("F:\\experiment"+seperator)
#define CONFIG_PATH  PRO_HOME_DIR
#define DATA_PATH PRO_HOME_DIR+std::string("data\\")
#define INDEX_PATH PRO_HOME_DIR+std::string("index\\")
#define QUERY_PATH PRO_HOME_DIR+std::string("query\\")
//#define QUERY_FILE_PATH PRO_HOME_DIR+std::string("queryfiles/")
//#define QUERY_RESULT_FILE_PATH PRO_HOME_DIR+std::string("queryresults/")


class ConfigType
{
private:
	// used to record all the config information
    std::map<std::string,std::string> cr; 
public:
	// receive the information from file or cmd
    ConfigType(std::string& configFileName, int argc,char** argv);
    ~ConfigType();
	// list all the config information of cr
    void ListConfig();
	// get map file name with key "map"
    std::string getMapFileName();
	// get data file name
    std::string getDataFileName();
	// return the int value of "cachepages"
    int getParameterCachePages();
	// return the int value of k 
    int getParameterK();
	// return int value of "querykeywordsnumber"
    int getParameterQueryKeywordNumbers();
	// return int value of "numberOfQueryPoints"
    int getParameterNumberOfQueryPoints();
	// return int value of "OutlierDensity"
    int getParameterOutlierDensity();
	// return int value of "AvgKeywordsNumberOfOutliers"
    int getParameterAvgKeywordsNumberOfOutliers();
	// return the query file correspoding to data filename
    std::string getQueryFileName();
	// return the queryresult file correspoding to query filename
    std::string getQueryResultFileName();
	//---------------------M--add the new function-------------
	// data relevant 
	// index relevant
	std::string getIndexFileName();
	// query relevant
	// return the query distance constraint
	float getParameterDistanceConstraint();
	// return the query subspace dimension
	float getParameterSubspaceDimensions();



private:
	// read the config information from "config.prop"
    void AddConfigFromFile(std::string& configFileName);
	// remove the space such as {'\t','\n','\f','\r',' '};
    void TrimSpace(char* str);
	// // read the config information from cmd
    void AddConfigFromCmdLine(int argc,char** argv);
	// return the float value of key,transfrom from str to float
    float getConfigFloat(std::string& key);
	// return the int value of key,transfrom from str to int
    int getConfigInt(std::string& key);
	// return the str value of key, no transform
    std::string getConfigStr(std::string& key);

};

#endif /* defined(__generateQueries__ConfigType__) */
