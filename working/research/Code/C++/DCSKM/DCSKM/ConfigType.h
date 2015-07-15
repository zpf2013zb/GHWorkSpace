//
//  ConfigType.h
//  generateQueries
//
//  Created by 秦旭 on 14-6-28.
//  Copyright (c) 2014年 ZJU. All rights reserved.
//

#ifndef __RKSK__ConfigType__
#define __RKSK__ConfigType__

#include <iostream>
#include <map>
#include <string>

#define PRO_HOME_DIR std::string("/Users/qinxu/Desktop/rksk_mac/")
#define CONFIG_PATH  PRO_HOME_DIR
#define DATA_PATH PRO_HOME_DIR+std::string("data/")
#define QUERY_FILE_PATH PRO_HOME_DIR+std::string("queryfiles/")
#define QUERY_RESULT_FILE_PATH PRO_HOME_DIR+std::string("queryresults/")


class ConfigType
{
private:
    std::map<std::string,std::string> cr;
public:
    ConfigType(std::string& configFileName, int argc,char** argv);
    ~ConfigType();
    void ListConfig();
    std::string getMapFileName();
    std::string getDataFileName();
    int getParameterCachePages();
    int getParameterK();
    int getParameterQueryKeywordNumbers();
    int getParameterNumberOfQueryPoints();
    int getParameterOutlierDensity();
    int getParameterAvgKeywordsNumberOfOutliers();
    std::string getQueryFileName();
    std::string getQueryResultFileName();
private:
    void AddConfigFromFile(std::string& configFileName);
    void TrimSpace(char* str);
    void AddConfigFromCmdLine(int argc,char** argv);
    float getConfigFloat(std::string& key);
    int getConfigInt(std::string& key);
    std::string getConfigStr(std::string& key);



};

#endif /* defined(__generateQueries__ConfigType__) */
