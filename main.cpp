#include <iostream>
#include "Vocabulary.h"
#include "Unigram.h"
#include "Model.h"
#include <string>
#include <sstream>

using namespace std;




int main(int argc, char** argv) {

    stringstream input_path, embedding_file, method_name;
    

    string dataset = "karate"; //"citeseer_undirected";

    input_path << "/home/abdulkadir/Desktop/expon/walks/" << dataset << "_node2vec_n=80_l=3.corpus";//"_n=10_l=80_p=1_q=1_method=n2v_sample=1.corpus";//"_node2vec_p=1_q=1_sample=1.corpus"; //"_n=80_l=10_method=dw_sample=1.corpus"; //"_node2vec_p=1_q=1.corpus"; //"_afaki.corpus"; //"_node2vec_p=1_q=1.corpus"; //_p=1_q=1
    embedding_file << "/home/abdulkadir/Desktop/expon/embeddings/" << dataset << "_method1_dim=2.embedding";
    method_name << "method1";

    /*
    if(argc == 4) {
        input_path << argv[1];
        embedding_file << argv[2];
        method_name << argv[3];
    } else {
        cout << "Format: ./run input.corpus output.embedding method_name" << endl;
        return 0;
    }
    */

    cout << input_path.str() << endl;
    cout << embedding_file.str() << endl;
    cout << method_name.str() << endl;

    int window_size = 10;
    int negative_sample_size = 5;
    int dim = 128;
    double starting_alpha = 0.025;
    double decay_rate = 1.0;
    double min_alpha = 0.0001;
    int num_iters = 1;


    Model model(input_path.str(), window_size, negative_sample_size, starting_alpha, decay_rate, min_alpha, num_iters, dim, method_name.str());
    model.run();
    model.save_embeddings(embedding_file.str());


    return 0;
}

