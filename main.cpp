#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <vector>
#include "Utilities.h"
#include "Model.h"


using namespace std;



//int argument_parser(int argc, char** argv,
//                     stringstream &input_path, stringstream &embedding_file, stringstream &method_name,
//                     double &starting_alpha, double &min_alpha, double &decay_rate, int &dim, int &negative_sample_size,
//                     int &window_size, int &num_iters) {
//
//
//    for(int i=1; i<argc; i=i+2) {
//
//
//        if(strcmp(argv[i], "--walks") == 0) {
//
//            input_path << argv[i+1];
//
//        } else if(strcmp(argv[i], "--output") == 0) {
//
//            embedding_file << argv[i+1];
//
//        } else if(strcmp(argv[i], "--method") == 0) {
//
//            method_name << argv[i+1];
//
//        } else if(strcmp(argv[i], "--starting_alpha") == 0) {
//
//            starting_alpha = stod(argv[i+1]);
//
//        }  else if(strcmp(argv[i], "--min_alpha") == 0) {
//
//            min_alpha = stod(argv[i+1]);
//
//        }  else if(strcmp(argv[i], "--decay_rate") == 0) {
//
//            decay_rate = stod(argv[i+1]);
//
//        }  else if(strcmp(argv[i], "--dim") == 0) {
//
//            dim = stoi(argv[i+1]);
//
//        }  else if(strcmp(argv[i], "--negative_sample_size") == 0) {
//
//            negative_sample_size = stoi(argv[i+1]);
//
//        }  else if(strcmp(argv[i], "--window_size") == 0) {
//
//            window_size = stoi(argv[i+1]);
//
//        }  else if(strcmp(argv[i], "--num_iters") == 0) {
//
//            num_iters = stoi(argv[i+1]);
//
//        } else if(strcmp(argv[i], "--help") == 0) {
//
//            string help_text;
//            help_text = "\nUsage: \n";
//            help_text += "\t./run --walks random_walks.corpus --output output_file.embedding --method method_name[bern, pois, norm]\n";
//            help_text +="\nOptional parameters:\n";
//            help_text += "\t--starting_alpha [Default: 0.025]\n\t --min_alpha [Default: 0.0001]\n\t--decay_rate: [Default: 1.0]\n";
//            help_text += "\t--dim [Default: 128]\n\t--negative_sample_size [Default: 5]\n\t--window_size [Default: 10]\n";
//            help_text += "\t--num_iters [Default: 1]\n";
//            cout << help_text << endl;
//
//            return -1;
//
//        } else {
//
//            cout << "You have typed an invalid argument name: " << argv[i] << endl;
//            cout << "For more details, you can use --help" << endl;
//
//            return -2;
//
//        }
//
//    }
//
//    // Check if the mandatory arguments are set
//    if(input_path.str().length() > 0 and embedding_file.str().length()>0 and method_name.str().length()>0 ) {
//        cout << "Input file: " << input_path.str() << endl;
//        cout << "Output file: " << embedding_file.str() << endl;
//        cout << "Method name: " << method_name.str() << endl;
//    } else {
//        cout << "Please enter the mandatory arguments: --walks CORPUS_FILE --output OUTPUT_FILE --method METHOD_NAME[bern, pois, norm]" << endl;
//
//        return -3;
//    }
//
//    return 0;
//
//}



int main(int argc, char** argv) {

    // Set the default values
    string corpus_file, output_file, method_name;
    stringstream deneme;
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

