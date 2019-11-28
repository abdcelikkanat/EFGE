#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <vector>
#include "Utilities.h"
#include "Model.h"


using namespace std;



int main(int argc, char** argv) {

    // Set the default values
    string corpus_file, output_file, method_name;
    string format;
    int window_size = 10;
    int negative_sample_size = 5;
    int dim = 128;
    double starting_alpha = 0.025;
    double min_alpha = 0.0001;
    double decay_rate = 1.0;
    int num_iters = 1;
    vector <double> optionalParams; optionalParams.push_back(1.0); // std_dev

    int err_code =  parse_arguments(argc, argv, corpus_file, output_file, method_name,
                                    window_size, negative_sample_size, dim,
                                    starting_alpha, min_alpha, decay_rate, num_iters,
                                    optionalParams);

    if(err_code != 0) {
        if(err_code < 0)
            cout << "Error code: " << err_code << endl;
        return 0;
    }

    Model model(corpus_file, method_name, starting_alpha, min_alpha, decay_rate, dim, negative_sample_size, window_size, num_iters, optionalParams);
    model.run();
    model.save_embeddings(output_file);


    return 0;
}

