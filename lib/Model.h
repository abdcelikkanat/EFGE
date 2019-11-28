//
//
//

#ifndef EFGE_MODEL_H
#define EFGE_MODEL_H

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
    string corpus_path;
    string method_name;

    int window_size, negative_sample_size, dim_size;
    double starting_alpha, min_alpha, decay_rate;
    int num_of_iters;
    double std_dev;

    int vocab_size;
    unordered_map <string, int> node2Id;
    vector <Node> vocab_items;
    int total_nodes;
    double **emb0, **emb1;

    default_random_engine generator;

    Unigram uni;

public:

    Model(string f_path, string method_name, double s_alpha, double m_alpha, double d_rate, int dim, int neg, int w_size, int num_iters, vector <double> &optionalParams);
    ~Model();
    double sigmoid(double z);
    void bernoulli_update(double alpha, vector <double> labels, int centerId, vector <int> contextIds);
    void poisson_update(double alpha, vector <double> labels, int centerId, vector <int> contextIds);
    void gaussian_known_var_exp(double alpha, vector <double> labels, int centerId, vector <int> contextIds);
    void run();
    void save_embeddings(string file_path);


};



#endif //EFGE_MODEL_H
