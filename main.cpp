#include <iostream>
#include "Vocabulary.h"
#include "Unigram.h"
#include "Model.h"
#include <string>
#include <sstream>

using namespace std;




int main(int argc, char** argv) {

    stringstream input_path, embedding_file, method_name;
    string format;
    double starting_alpha = 0.025;
    double min_alpha = 0.0001;
    double decay_rate = 1.0;
    int dim = 128;
    int negative_sample_size = 5;
    int window_size = 10;
    int num_iters = 1;


    if(argc >= 4) {
        input_path << argv[1];
        embedding_file << argv[2];
        method_name << argv[3];

        if(argc >= 5)
            starting_alpha = stod(argv[4]);
        if(argc >= 6)
            min_alpha = stod(argv[5]);
        if(argc >= 7)
            decay_rate = stod(argv[6]);
        if(argc >= 8)
            dim = stoi(argv[7]);
        if(argc >= 9)
            negative_sample_size = stoi(argv[8]);
        if(argc >= 10)
            window_size = stoi(argv[9]);
        if(argc == 11)
            num_iters = stoi(argv[10]);

    } else {
        format = "\nUsage: \n";
        format += "\t./expfamemb input_file.corpus output_file.embedding method_name[bern, pois, norm]\n";
        format +="\nOptional parameters:\n";
        format += "\tStarting alpha [Default: 0.025]\n\tMinimum alpha [Default: 0.0001]\n\tDecay rate: [Default: 1.0]\n";
        format += "\tDimension size [Default: 128]\n\tNegative sample size [Default: 5]\n\tWindow size [Default: 10]\n";
        format += "\tNumber of iterations [Default: 1]\n";
        cout << format << endl;
        return 0;
    }

    cout << "Input file: " << input_path.str() << endl;
    cout << "Output file: " << embedding_file.str() << endl;
    cout << "Method name: " << method_name.str() << endl;




    Model model(input_path.str(), method_name.str(), starting_alpha, min_alpha, decay_rate, dim, negative_sample_size, window_size, num_iters);
    model.run();
    model.save_embeddings(embedding_file.str());


    return 0;
}

