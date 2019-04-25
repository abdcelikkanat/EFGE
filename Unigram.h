//
// Created by abdulkadir on 13/11/18.
//

#ifndef FAST_BERN_UNIGRAM_H
#define FAST_BERN_UNIGRAM_H

#include <math.h>
#include <random>
#include <iostream>

using namespace std;

const unsigned int UNIGRAM_TABLE_SIZE = 1000000;



class Unigram {
private:
    int table[UNIGRAM_TABLE_SIZE] = {0};

    default_random_engine generator;



public:
    Unigram();
    Unigram(int size, vector<int> freq, double expon);
    ~Unigram();

    void sample(int count, int samples[]);

};


#endif //FAST_BERN_UNIGRAM_H
