//
//  KeywordsGenerator.cpp
//  rksk_query
//
//  Created by 秦旭 on 14-7-1.
//  Copyright (c) 2014年 ZJU. All rights reserved.
//

#include "KeywordsGenerator.h"
#include <bitset>
#include <vector>

#define MAX_KEYWORDS 64

std::vector<unsigned long long> KeywordsGenerator::getKeywords(std::size_t totalNumberOfNodes, std::size_t avgKeywords)
{
    
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(avgKeywords,0.01);
    
    std::vector<unsigned long long> tmp;
    tmp.reserve(totalNumberOfNodes);
    std::bitset<MAX_KEYWORDS> keywordsSet;
    
    while (totalNumberOfNodes--) {
        keywordsSet.reset();
        int numberOfKeywords = distribution(generator);
        int percent = random()%1000;
        // precent used for what?
        if (percent < 50) {
            
            while(keywordsSet.count() < numberOfKeywords)
            {
                int index = random()%64;
                keywordsSet.set(index);
                
            }
            
        }
        else{
            
            while(keywordsSet.count() < numberOfKeywords)
            {
                int index;
                if (keywordsSet.count() < 8) {
                    
                    index = random()%8;
                }
                else
                    index = random()%20;
                
                keywordsSet.set(index+22);
                
            }
            
        }
        
        tmp.push_back(keywordsSet.to_ullong());

    }
    
    return tmp;
}

std::vector<unsigned long long> KeywordsGenerator::getConstantKeywords(std::size_t totalNumberOfNodes, std::size_t number)
{
    std::vector<unsigned long long> keys;
    keys.reserve(totalNumberOfNodes);
    std::bitset<MAX_KEYWORDS> keywords_set;
    while (totalNumberOfNodes --) {
        
        keywords_set.reset();
        
        while(keywords_set.count() < number )
        {
            int index;
            if (keywords_set.count() < 8) {
                index = random()%8;
            }
            else
                index = random()%20;
            
            keywords_set.set(index+22);
            
        }//Avoid generate vct#=0
        
        keys.push_back(keywords_set.to_ullong());
    }
    
    return keys;
    
}

KeywordsGenerator::KeywordsGenerator()
{
}

KeywordsGenerator::~KeywordsGenerator()
{
}


