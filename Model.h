//
// Created by abdulkadir on 23/04/19.
//

#ifndef FAST_BERN_MODEL_H
#define FAST_BERN_MODEL_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <random>
#include <iomanip>
#include "math.h"
#include "Unigram.h"
#include "Vocabulary.h"
#include <algorithm>

using namespace std;

class Model {

private:
    string method_name;
    string file_path;
    int window_size, negative_sample_size, dim_size;
    unsigned long vocab_size;
    unordered_map <string, int> node2Id;
    vector <Node> vocab_items;
    int total_nodes;

    double starting_alpha, decay_rate, min_alpha;
    int num_of_iters;

    double **emb0, **emb1;

    default_random_engine generator;

    Unigram uni;

public:

    Model(string f_path, int w_size, int neg, double s_alpha, double d_rate, double m_alpha, int num_iters, int dim, string method_name);
    ~Model();
    vector <unordered_map <int, int>> getCoOccurenceCount();
    vector <unordered_map <int, double>> getRelativeFreq(vector <unordered_map <int, int>> cooccurences);
    vector <unordered_map <int, double>> getNormalizedFreq(vector <unordered_map <int, int>> cooccurences);
    double sigmoid(double z);
    void bernoulli_update(double alpha, vector <double> labels, int centerId, vector <int> contextIds);
    void bernoulli_update_v2(double alpha, vector <double> labels, int centerId, vector <int> contextIds);
    void poisson_update_v1(double alpha, vector <double> labels, int centerId, vector <int> contextIds);
    void poisson_update_v2(double alpha, vector <double> labels, int centerId, vector <int> contextIds);
    void exponential_update_v1(double alpha, vector <double> labels, int centerId, vector <int> contextIds);
    void gaussian_known_var(double alpha, vector <double> labels, int centerId, vector <int> contextIds);
    void gaussian_my_prior(double alpha, vector <double> labels, int centerId, vector <int> contextIds);
    void poisson_update_v3(double alpha, vector <double> labels, int centerId, vector <int> contextIds);
    void poisson_update_v4(double alpha, vector <double> labels, int centerId, vector <int> contextIds);
    void gaussian_kernel_update(double alpha, vector <double> labels, int centerId, vector <int> contextIds, double std);
    void run();
    void save_embeddings(string file_path);


};



#endif //FAST_BERN_MODEL_H
