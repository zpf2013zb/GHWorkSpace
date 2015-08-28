//
//  ConfigType.cpp
//  generateQueries
//
//  Created by 秦旭 on 14-6-28.
//  Copyright (c) 2014年 ZJU. All rights reserved.
//

#include "ConfigType.h"
#include <fstream>
#include <string>
#include <cstring>

ConfigType::ConfigType(std::string &configFileName,int argc,char** argv)
{
    AddConfigFromFile(configFileName);
    AddConfigFromCmdLine(argc, argv);
}

ConfigType::~ConfigType()
{

}

void ConfigType::AddConfigFromFile(std::string& configFileName)
{
    const int LINE_LEN=1024; // max length of a line
    char line[LINE_LEN],key[LINE_LEN],value[LINE_LEN];

    std::ifstream br(CONFIG_PATH + configFileName);

    if (!br.is_open())
    {
        std::cerr<<"Error opening file "<< CONFIG_PATH + "config.prop";
        exit (1);
    }

    while (br.getline(line,LINE_LEN))
    {
        if (strstr(line,"//")!=NULL || strstr(line,"/")!=NULL ||strstr(line, "*")!=NULL) continue; // remove comments
        char* chPos=strchr(line,'=');
        if (chPos!=NULL)
        {
            int pos=((int)(chPos-line))/sizeof(char);
            int keyLen=pos;
            unsigned long valueLen=strlen(line)-1-keyLen;
            memcpy(key,&line[0],keyLen);
            key[keyLen]='\0';
            memcpy(value,&line[pos+1],valueLen);
            value[valueLen]='\0';
            TrimSpace(key);
            TrimSpace(value);
            cr[key]=value;
        }
    }
    br.close();
}

void ConfigType::AddConfigFromCmdLine(int argc,char** argv)
{
    int i=0;
    while (i<argc)
    {
        while ((i<argc)&&(argv[i][0]!='-')) i++;
        if (i+1<argc)
        {
            char* key=&(argv[i][1]); //??取址
            char* value=argv[i+1];
            TrimSpace(key);
            TrimSpace(value);
            cr[key]=value;
            i+=2;
        }
        else
            return;
    }
}

void ConfigType::ListConfig()
{
    std::map<std::string, std::string>::iterator p=cr.begin();
    while (p!=cr.end())
    {
        printf("%s=%s\n",p->first.c_str(),p->second.c_str());
        p++;
    }
}

float ConfigType::getConfigFloat(std::string &key)
{
    float value = 0;
    if (cr.count(key))
        value=atof(cr[key].c_str());
    else
    {
        std::cerr<<"Config key not found :"<<key;
        exit(1);

    }
    return value;
}

int ConfigType::getConfigInt(std::string& key)
{
    int value = 0;
    if (cr.count(key))
        value=atoi(cr[key].c_str());
    else
    {

        std::cerr<<"Config key not found :"<<key;
        exit(1);

    }
    return value;
}

std::string ConfigType::getConfigStr(std::string& key)
{
    std::string value;
    if (cr.count(key))
        value=cr[key];
    else
    {
        std::cerr<<"Config key  not found :"<<key;
        exit(1);
    }
    return value;
}

void ConfigType::TrimSpace(char* str)
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

std::string ConfigType::getMapFileName()
{
    std::string key = "map";
    return DATA_PATH + getConfigStr(key);
}

std::string ConfigType::getDataFileName()
{
    return getMapFileName()+"_data_"+std::to_string(getParameterOutlierDensity())+"_"+std::to_string(getParameterAvgKeywordsNumberOfOutliers());
}

int ConfigType::getParameterCachePages()
{
    std::string key = "cachepages";
    return getConfigInt(key);
}

int ConfigType::getParameterK()
{
    std::string key = "k";
    return getConfigInt(key);
}
int ConfigType::getParameterQueryKeywordNumbers()
{
    std::string key = "querykeywordsnumber";
    return getConfigInt(key);

}
std::string ConfigType::getQueryFileName()
{
    auto pos = getDataFileName().rfind(seperator);
    std::string dataName = getDataFileName().substr(pos+1);
    std::string queryFileName = QUERY_PATH + seperator + "query"+dataName+"_queryPoints_"+ std::to_string(getParameterQueryKeywordNumbers())+"_"+std::to_string(getParameterK())+"_"+std::to_string(getParameterCachePages());
    return queryFileName;
}

int ConfigType::getParameterNumberOfQueryPoints()
{
    std::string key = "numberOfQueryPoints";
    return getConfigInt(key);
}

std::string ConfigType::getQueryResultFileName()
{
    auto pos = getQueryFileName().rfind(seperator);
    
    return QUERY_PATH + seperator + "result" + getQueryFileName().substr(pos+1)+std::string("_ResultFile");
}

std::string ConfigType::getIndexFileName()
{
	std::string key = "map";
	return INDEX_PATH + getConfigStr(key);
}

float ConfigType::getParameterDistanceConstraint()
{
	std::string key = "querydistanceconstraint";
	return getConfigInt(key);
}

float ConfigType::getParameterSubspaceDimensions()
{
	std::string key = "querysubspacedimensions";
	return getConfigInt(key);
}

int ConfigType::getParameterOutlierDensity()
{
    std::string key = "OutlierDensity";
    return getConfigInt(key);
}
int ConfigType::getParameterAvgKeywordsNumberOfOutliers()
{
    std::string key = "AvgKeywordsNumberOfOutliers";
    return getConfigInt(key);
}
