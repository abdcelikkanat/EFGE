//
// Created by abdulkadir on 13/11/18.
//

#ifndef FAST_BERN_MODEL_SON_H
#define FAST_BERN_MODEL_SON_H

#include "Graph.h"
#include <vector>
#include "math.h"
#include <random>
#include <sstream>
#include "Unigram.h"
#include "string.h"
#include <chrono>
#include <algorithm>

using namespace std;

class Model {
    Graph g;
    unsigned int num_of_nodes = 0;
    int dim_size = 0;
    vector <vector <int> > adj_list;
    vector <vector <int>> nb_list;
    default_random_engine generator;


    double **emb0, **emb1;


public:
    Model(Graph g, int dim);
    ~Model();

    double sigmoid(double z);
    void getNeighbors();
    void readGraph(string file_path, string filetype, bool directed);
    void run(double starting_alpha, double min_alpha, double decay_rate, int num_of_iters, int negative, int save_step, string save_file);
    void save_embeddings(string file_path);

    void computeTriangles(vector <vector <int>> &triangles, vector <vector <int>> &triangle_counts);

    void getNeighborNodes(vector <vector <int>> &nb_list, vector <vector <int>> &counts);
    void nodeseq2degreeseq(vector <vector <int>> &nb_list, vector <vector <int>> &degree_seq);
};


#endif //FAST_BERN_MODEL_SON_H

