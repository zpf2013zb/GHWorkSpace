//
//  KeywordsGenerator.h
//  rksk_query
//
//  Created by 秦旭 on 14-7-1.
//  Copyright (c) 2014年 ZJU. All rights reserved.
//

#ifndef __rksk_query__KeywordsGenerator__
#define __rksk_query__KeywordsGenerator__

#include <iostream>
#include <random>

class KeywordsGenerator
{
public:
    static KeywordsGenerator& Instance()
    {
        static KeywordsGenerator theGenerator;
        return theGenerator;
    }
    
    static std::vector<unsigned long long> getKeywords(std::size_t totalNumberOfNodes,std::size_t avgKeywords);
    
    static std::vector<unsigned long long> getConstantKeywords(std::size_t totalNumberOfNodes, std::size_t number);
    
private:
    KeywordsGenerator();
    ~KeywordsGenerator();
    KeywordsGenerator(const KeywordsGenerator&);
    KeywordsGenerator& operator = (const KeywordsGenerator&);
    
};

#endif /* defined(__rksk_query__KeywordsGenerator__) */
