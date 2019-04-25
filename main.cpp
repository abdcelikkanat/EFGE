#include <iostream>
#include "Vocabulary.h"
#include "Unigram.h"
#include "Model.h"
#include <string>
#include <sstream>

using namespace std;




int main() {

    stringstream input_path, embedding_file;
    string dataset = "citeseer_undirected";

    input_path << "/home/abdulkadir/Desktop/expon/walks/" << dataset << "_node2vec_p=1_q=1.corpus"; //_p=1_q=1
    embedding_file << "/home/abdulkadir/Desktop/expon/embeddings/" << dataset << "_test.embedding";


    int window_size = 10;
    int negative_sample_size = 5;
    int dim = 128;
    double starting_alpha = 0.025;
    double decay_rate = 1.0;
    double min_alpha = 0.0001;
    int num_iters = 1;


    Model model(input_path.str(), window_size, negative_sample_size, starting_alpha, decay_rate, min_alpha, num_iters, dim);
    model.run();
    model.save_embeddings(embedding_file.str());
    //model.getCoOccurenceCount();

    //Vocabulary(corpus_path.str());

    /*
    stringstream graph_path, embedding;
    string dataset = "citeseer_undirected";

    graph_path << "/home/abdulkadir/Desktop/datasets/" << dataset << ".edgelist";
    //graph_path << "./" << dataset << ".edgelist";

    Graph g;
    g.readGraph(graph_path.str(), "edgelist", false);

    cout << "Number of nodes: " << g.getNumOfNodes() << endl;
    cout << "Number of edges: " << g.getNumOfEdges() << endl;

    double starting_alpha=0.025;
    double min_alpha=0.0001;
    double decay_rate=0.001;
    int num_of_iters=1000;
    int negative=5;
    int save_step=10;
    string save_file="/home/abdulkadir/Desktop/expemb/embeddings/";

    Model model(g, 128);
    model.run(starting_alpha, min_alpha, decay_rate, num_of_iters, negative, save_step, save_file);
    embedding << save_file << "deneme" << "_last.embedding";
    //embedding << "./" << dataset << "_last.embedding";
    
    model.save_embeddings(embedding.str());
    */

    return 0;
}

