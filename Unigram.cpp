//
//
//

#include "Unigram.h"



Unigram::Unigram() {

}



Unigram::Unigram(int size, vector<int> freq, double expon) {

    cout << "--> The construction of unigram table has started." << endl;

    double norm = 0.0;

    for(int j=0; j<size; j++) {
        norm += pow(freq[j], expon);
    }

    double p=0;
    int i=0;
    for(int j=0; j<size; j++) {
        p += pow((double)freq[j], expon) / norm;
        while(i < UNIGRAM_TABLE_SIZE && (double)i/UNIGRAM_TABLE_SIZE < p ) {
            table[i] = j;
            i++;
        }
    }

    cout << "    + Completed!" << endl;

}

Unigram::~Unigram() {

}


void Unigram::sample(int count, int samples[]) {

    uniform_int_distribution<int> uni(0, UNIGRAM_TABLE_SIZE-1);


    for(int i=0; i<count; i++) {
        samples[i] = table[uni(generator)];
    }
}

